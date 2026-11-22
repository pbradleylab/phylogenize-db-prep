"""Run in the matrix conda environment, or set PHYLOGENIZE_RSCRIPT and PHYLOGENIZE_MATRIX_PYTHON."""

import csv
import io
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
RSCRIPT = os.environ.get("PHYLOGENIZE_RSCRIPT") or shutil.which("Rscript")
MATRIX_PYTHON = os.environ.get("PHYLOGENIZE_MATRIX_PYTHON", sys.executable)


@unittest.skipUnless(RSCRIPT, "Rscript from the matrix environment is required")
class BinaryMatrixTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="binary matrices ")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.mapping = self.directory / "protein members.tsv"
        self.metadata = self.directory / "genome metadata.tsv"
        self.hits = self.directory / "combined hits.tsv"
        self.pairs = self.directory / "prepared" / "species families.csv"
        self.taxonomy = self.directory / "prepared" / "taxonomy.csv"
        self.binary = self.directory / "binary" / "profiles.rds"
        self.mapping.write_text(
            "c_first\tMGYG004_1\n"
            "c_shared\tMGYG002_1\n"
            "c_shared\tMGYG001.1_2\n"
            "c_second\tMGYG003_1\n"
            "c_second\tMGYG001.1_3\n"
            "c_shared\tMGYG001.1_4\n"
            "c_na\tMGYG005_1\n"
        )
        self.hits.write_text("query\ttarget\nc_first\t0999\nc_shared\t0002\nc_second\t0010\nc_na\tNA\n")
        self.metadata.write_text("Genome\tSpecies_rep\tLineage\n" + "".join(
            f"{genome}\t{species}\td__Bacteria;p__{phylum};c__C;o__O;f__F;g__G;s__S\n"
            for genome, species, phylum in [
                ("MGYG001.1", "MGYG001.1", "Firmicutes_A"),
                ("MGYG002.1", "MGYG001.1", "Firmicutes_A"),
                ("MGYG003.1", "MGYG003.1", "Firmicutes_A"),
                ("MGYG004.1", "MGYG004.1", ""),
                ("MGYG005.1", "MGYG005.1", "NA"),
            ]
        ))
        self.expected_pairs = {
            ("MGYG001.1", "0002"), ("MGYG001.1", "0010"),
            ("MGYG003.1", "0010"), ("MGYG004.1", "0999"), ("MGYG005.1", "NA"),
        }

    def run_command(self, arguments, error=None):
        result = subprocess.run(arguments, capture_output=True, text=True, timeout=90)
        details = result.stdout + result.stderr
        if error:
            self.assertNotEqual(result.returncode, 0, details)
            self.assertIn(error, details)
        else:
            self.assertEqual(result.returncode, 0, details)
        return result

    def prepare(self, *arguments, error=None):
        self.run_command([
            MATRIX_PYTHON, str(REPOSITORY / "workflow/scripts/reformat_taxonomy.py"),
            "--input", str(self.hits), "--mapping", str(self.mapping), "--tax", str(self.metadata),
            "--output", str(self.pairs), "--tax_output", str(self.taxonomy), *arguments,
        ], error=error)
        if error:
            self.assertFalse(self.pairs.exists())
            self.assertFalse(self.taxonomy.exists())

    def build(self, error=None):
        self.run_command([
            RSCRIPT, str(REPOSITORY / "workflow/scripts/make_binary.R"),
            str(self.pairs), str(self.taxonomy), str(self.binary),
        ], error=error)
        if error:
            self.assertFalse(self.binary.exists())

    def assert_pairs(self):
        with self.pairs.open() as handle:
            reader = csv.DictReader(handle)
            self.assertEqual(reader.fieldnames, ["cluster", "target"])
            rows = list(reader)
        self.assertEqual(len(rows), len(self.expected_pairs))
        self.assertEqual({(row["cluster"], row["target"]) for row in rows}, self.expected_pairs)

    def test_member_genomes_and_first_headerless_record_are_retained(self):
        self.prepare()
        self.assert_pairs()
        with self.taxonomy.open() as handle:
            taxonomy = {row["cluster"]: row["phylum"] for row in csv.DictReader(handle)}
        self.assertEqual(taxonomy, {
            "MGYG001.1": "Firmicutes_A", "MGYG003.1": "Firmicutes_A",
            "MGYG004.1": "unclassified", "MGYG005.1": "NA",
        })

    def test_binary_profiles_retain_singletons_and_correct_family_labels(self):
        self.prepare()
        self.build()
        result = self.run_command([RSCRIPT, "-e", '''
            library(Matrix)
            matrices <- readRDS(commandArgs(TRUE)[1])
            stopifnot(setequal(names(matrices), c("Firmicutes_A", "unclassified", "NA")))
            entries <- lapply(names(matrices), function(phylum) {
              m <- matrices[[phylum]]
              stopifnot(inherits(m, "lgCMatrix"), !anyNA(m@x), all(m@x),
                        !anyDuplicated(rownames(m)), !anyDuplicated(colnames(m)))
              x <- summary(m)
              data.frame(cluster=colnames(m)[x$j], target=rownames(m)[x$i])
            })
            stopifnot(!matrices[["Firmicutes_A"]]["0002", "MGYG003.1"])
            write.table(do.call(rbind, entries), stdout(), sep="\\t", quote=FALSE, row.names=FALSE)
        ''', str(self.binary)])
        self.assertEqual({
            (row["cluster"], row["target"])
            for row in csv.DictReader(io.StringIO(result.stdout), delimiter="\t")
        }, self.expected_pairs)

    def test_binary_continuous_and_tree_species_agree(self):
        self.prepare()
        self.build()
        continuous = self.directory / "continuous.rds"
        tree_input = self.directory / "source tree.nwk"
        tree_output = self.directory / "trees.rds"
        tree_input.write_text("((MGYG001:0.3,MGYG003:0.2):0.2,(MGYG004:0.4,MGYG005:0.6):0.1);\n")
        self.run_command([
            RSCRIPT, str(REPOSITORY / "workflow/scripts/sparse_continuous_pangenomes.R"),
            "--input", str(self.mapping), "--genome_metadata", str(self.metadata),
            "--combined_hits", str(self.hits), "--output_rds", str(continuous),
            "--output_csv", str(self.directory / "continuous.csv"), "--memory", "1", "--speedy", "TRUE",
        ])
        self.run_command([
            RSCRIPT, str(REPOSITORY / "workflow/scripts/make_tree.R"), "--input", str(tree_input),
            "--taxonomy", str(self.taxonomy), "--output", str(tree_output),
        ])
        self.run_command([RSCRIPT, "-e", '''
            library(Matrix)
            args <- commandArgs(TRUE)
            binary <- readRDS(args[1]); continuous <- readRDS(args[2]); trees <- readRDS(args[3])
            stopifnot(setequal(names(binary), names(continuous)), setequal(names(binary), names(trees)))
            for (phylum in names(binary)) {
              b <- binary[[phylum]]; c <- continuous[[phylum]]
              stopifnot(setequal(colnames(b), colnames(c)), setequal(rownames(b), rownames(c)),
                        setequal(colnames(b), trees[[phylum]]$tip.label))
              c <- c[rownames(b), colnames(b), drop=FALSE]
              stopifnot(Matrix::nnzero(b - (c > 0)) == 0)
            }
            stopifnot(continuous[["Firmicutes_A"]]["0010", "MGYG001.1"] == 0.5)
        ''', str(self.binary), str(continuous), str(tree_output)])

    def test_annotation_linker_accepts_species_family_pairs(self):
        self.prepare()
        annotations = self.directory / "annotations.tsv"
        annotations.write_text("gene_callers_id\tlinker_info\taccession\tfunction\n" + "".join(
            f"{index}\t{family}\tKO{index:05}\ttest function\n"
            for index, family in enumerate(["0002", "0010", "0999", "NA"])
        ))
        output = self.directory / "functions.csv"
        self.run_command([
            MATRIX_PYTHON, str(REPOSITORY / "workflow/scripts/sequence_to_annotation_linker.py"),
            str(self.hits), str(annotations), str(self.pairs), str(output),
        ])
        with output.open() as handle:
            self.assertEqual({(row["linker_info"], row["node_head"]) for row in csv.DictReader(handle)}, self.expected_pairs)

    def test_duplicate_hits_and_members_do_not_duplicate_pairs(self):
        self.mapping.write_text(self.mapping.read_text() + "c_shared\tMGYG001.1_2\n")
        self.hits.write_text(self.hits.read_text() + "c_shared\t0002\n")
        self.prepare()
        self.assert_pairs()

    def test_literal_custom_delimiter_is_supported(self):
        self.mapping.write_text("".join(
            centroid + "\t" + member.rsplit("_", 1)[0] + "|" + member.rsplit("_", 1)[1] + "\n"
            for centroid, member in (line.split("\t") for line in self.mapping.read_text().splitlines())
        ))
        self.prepare("--split_char", "|")
        self.assert_pairs()

    def test_unmatched_member_genome_is_rejected(self):
        self.mapping.write_text(self.mapping.read_text() + "c_shared\tMGYG999_1\n")
        self.prepare(error="Member genomes missing from metadata")

    def test_unmapped_centroid_is_rejected(self):
        self.hits.write_text(self.hits.read_text().replace("c_na\tNA\n", ""))
        self.prepare(error="Protein centroids missing from combined hits")

    def test_conflicting_family_assignments_are_rejected(self):
        self.hits.write_text(self.hits.read_text() + "c_shared\t0010\n")
        self.prepare(error="maps to more than one target family")

    def test_ambiguous_genome_versions_are_rejected(self):
        row = self.metadata.read_text().splitlines()[1].replace("MGYG001.1\t", "MGYG001.2\t", 1)
        self.metadata.write_text(self.metadata.read_text() + row + "\n")
        self.prepare(error="ambiguous after version normalization")

    def test_missing_source_species_is_rejected(self):
        self.mapping.write_text(self.mapping.read_text().replace("c_na\tMGYG005_1\n", ""))
        self.prepare(error="Source species missing from member mappings")

    def test_missing_representative_metadata_is_rejected(self):
        self.metadata.write_text(self.metadata.read_text().replace("MGYG003.1\tMGYG003.1", "MGYG003.1\tMGYG999.1"))
        self.prepare(error="Species representatives missing from metadata")

    def test_legacy_intermediate_is_rejected(self):
        self.prepare()
        self.pairs.write_text("query,target,other,cluster\nq,t,m,MGYG001.1\n")
        self.build(error="regenerate legacy intermediates")

    def test_missing_matrix_species_is_rejected(self):
        self.prepare()
        self.pairs.write_text(self.pairs.read_text().replace("MGYG005.1,NA\n", ""))
        self.build(error="exactly the same unique species IDs")

    def test_duplicate_taxonomy_species_is_rejected(self):
        self.prepare()
        self.taxonomy.write_text(self.taxonomy.read_text() + self.taxonomy.read_text().splitlines()[1] + "\n")
        self.build(error="exactly the same unique species IDs")

    def test_blank_family_id_is_rejected(self):
        self.prepare()
        self.pairs.write_text(self.pairs.read_text().replace(",0999\n", ",\n"))
        self.build(error="empty or missing IDs")


if __name__ == "__main__":
    unittest.main()
