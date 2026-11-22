"""CLI regression tests; set PHYLOGENIZE_RSCRIPT to an Rscript in the matrix environment."""

import csv
import io
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
RSCRIPT = os.environ.get("PHYLOGENIZE_RSCRIPT") or shutil.which("Rscript")


@unittest.skipUnless(RSCRIPT, "Rscript from the matrix environment is required")
class ContinuousPangenomeTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="continuous matrices ")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.mapping = self.directory / "protein members.tsv"
        self.metadata = self.directory / "genome metadata.tsv"
        self.hits = self.directory / "combined hits.tsv"
        self.output_rds = self.directory / "output" / "continuous.rds"
        self.output_csv = self.directory / "output" / "continuous.csv"
        lineages = {
            "Firmicutes_A": "d__Bacteria;p__Firmicutes_A;c__C;o__O;f__F;g__G;s__S",
            "blank": "d__Bacteria;p__;c__C;o__O;f__F;g__G;s__S",
            "NA": "d__Bacteria;p__NA;c__C;o__O;f__F;g__G;s__S",
        }
        self.metadata.write_text(
            "Genome\tSpecies_rep\tLineage\n"
            + "".join(
                f"{genome}\t{species}\t{lineages[phylum]}\n"
                for genome, species, phylum in [
                    ("MGYG003.1", "MGYG003.1", "Firmicutes_A"),
                    ("MGYG002.1", "MGYG001.1", "Firmicutes_A"),
                    ("MGYG001.1", "MGYG001.1", "Firmicutes_A"),
                    ("MGYG004.1", "MGYG004.1", "blank"),
                    ("MGYG005.1", "MGYG005.1", "NA"),
                ]
            )
        )
        self.mapping.write_text(
            "c_shared\tMGYG002_1\n"
            "c_shared\tMGYG001.1_2\n"
            "c_shared\tMGYG001.1_2\n"
            "c_second\tMGYG001.1_3\n"
            "c_second\tMGYG001.1_4\n"
            "c_other\tMGYG003_1\n"
            "c_blank\tMGYG004_1\n"
            "c_na\tMGYG005_1\n"
        )
        self.hits.write_text(
            "query\ttarget\n"
            "c_second\t0010\n"
            "c_shared\t0002\n"
            "c_other\t0002\n"
            "c_blank\t0999\n"
            "c_na\t0999\n"
        )

    def run_script(self, *arguments, error=None):
        result = subprocess.run(
            [
                RSCRIPT,
                str(REPOSITORY / "workflow/scripts/sparse_continuous_pangenomes.R"),
                "--input", str(self.mapping),
                "--genome_metadata", str(self.metadata),
                "--combined_hits", str(self.hits),
                "--output_rds", str(self.output_rds),
                "--output_csv", str(self.output_csv),
                "--memory", "1", "--threads", "1",
                *arguments,
            ],
            capture_output=True, text=True, timeout=120,
        )
        details = result.stdout + result.stderr
        if error:
            self.assertNotEqual(result.returncode, 0, details)
            self.assertIn(error, details)
            self.assertFalse(self.output_rds.exists())
        else:
            self.assertEqual(result.returncode, 0, details)

    def assert_outputs(self):
        with self.output_csv.open() as handle:
            rows = list(csv.DictReader(handle))
        observed = {
            (row["cluster"], row["column0"]): float(row["frac_observed"])
            for row in rows
        }
        self.assertEqual(observed, {
            ("MGYG001.1", "0002"): 1.0,
            ("MGYG001.1", "0010"): 0.5,
            ("MGYG003.1", "0002"): 1.0,
            ("MGYG004.1", "0999"): 1.0,
            ("MGYG005.1", "0999"): 1.0,
        })
        inspection = subprocess.run(
            [RSCRIPT, "-e", '''
                library(Matrix)
                matrices <- readRDS(commandArgs(TRUE)[1])
                stopifnot(setequal(names(matrices), c("Firmicutes_A", "unclassified", "NA")))
                entries <- lapply(names(matrices), function(phylum) {
                  m <- matrices[[phylum]]
                  stopifnot(inherits(m, "sparseMatrix"),
                            identical(rownames(m), sort(rownames(m))),
                            identical(colnames(m), sort(colnames(m))))
                  x <- summary(m)
                  data.frame(cluster=colnames(m)[x$j], family=rownames(m)[x$i], fraction=x$x)
                })
                stopifnot(matrices[["Firmicutes_A"]]["0010", "MGYG003.1"] == 0)
                write.table(do.call(rbind, entries), stdout(), sep="\\t", quote=FALSE, row.names=FALSE)
            ''', str(self.output_rds)],
            capture_output=True, text=True, timeout=30,
        )
        self.assertEqual(inspection.returncode, 0, inspection.stdout + inspection.stderr)
        matrix_rows = csv.DictReader(io.StringIO(inspection.stdout), delimiter="\t")
        self.assertEqual({
            (row["cluster"], row["family"]): float(row["fraction"])
            for row in matrix_rows
        }, observed)

    def test_fractional_presence_in_both_execution_modes(self):
        for speedy in ("FALSE", "TRUE"):
            with self.subTest(speedy=speedy):
                self.run_script("--speedy", speedy)
                self.assert_outputs()

    def test_configured_delimiter_is_literal(self):
        self.mapping.write_text(self.mapping.read_text().replace("_1\n", "|1\n").replace("_2\n", "|2\n").replace("_3\n", "|3\n").replace("_4\n", "|4\n"))
        self.run_script("--split_char", "|")
        self.assert_outputs()

    def test_unmapped_centroid_is_rejected(self):
        self.hits.write_text(self.hits.read_text().replace("c_other\t0002\n", ""))
        self.run_script(error="Protein centroids missing from combined hits")

    def test_conflicting_centroid_assignments_are_rejected(self):
        self.hits.write_text(self.hits.read_text() + "c_shared\t0010\n")
        self.run_script(error="more than one target family")

    def test_unknown_member_genome_is_rejected(self):
        self.mapping.write_text(self.mapping.read_text() + "c_shared\tMGYG999_1\n")
        self.run_script(error="Member genomes missing from metadata")

    def test_ambiguous_genome_versions_are_rejected(self):
        row = self.metadata.read_text().splitlines()[3].replace("MGYG001.1\t", "MGYG001.2\t", 1)
        self.metadata.write_text(self.metadata.read_text() + row + "\n")
        self.run_script(error="ambiguous IDs after version normalization")

    def test_species_without_members_is_rejected(self):
        row = self.metadata.read_text().splitlines()[1].replace("MGYG003.1", "MGYG999.1")
        self.metadata.write_text(self.metadata.read_text() + row + "\n")
        self.run_script(error="does not cover every source species")

    def test_denominator_includes_all_source_genomes(self):
        row = self.metadata.read_text().splitlines()[2].replace("MGYG002.1", "MGYG006.1")
        self.metadata.write_text(self.metadata.read_text() + row + "\n")
        self.run_script()
        with self.output_csv.open() as handle:
            rows = {row["column0"]: row for row in csv.DictReader(handle) if row["cluster"] == "MGYG001.1"}
        self.assertAlmostEqual(float(rows["0002"]["frac_observed"]), 2 / 3)
        self.assertAlmostEqual(float(rows["0010"]["frac_observed"]), 1 / 3)


if __name__ == "__main__":
    unittest.main()
