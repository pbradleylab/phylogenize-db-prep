"""Tree CLI regression tests; set PHYLOGENIZE_RSCRIPT to the matrix environment's Rscript."""

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
class MakeTreeTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="tree construction ")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.primary = self.directory / "bacterial tree.nwk"
        self.alternate = self.directory / "archaeal tree.nwk"
        self.taxonomy = self.directory / "taxonomy.csv"
        self.output = self.directory / "output" / "trees.rds"
        self.primary.write_text("(((B1:0.2,B2.1:0.3):0.4,B3:0.8):0.2,(B4:0.5,B_out:0.7):0.3);\n")
        self.alternate.write_text("((A1:0.1,A2.1:0.4):0.2,A3:0.8);\n")
        self.taxonomy.write_text(
            "cluster,phylum\n"
            "B1.1,Firmicutes_A\n"
            "B2.1,Firmicutes_A\n"
            "B3.1,NA\n"
            "B4.1,Single\n"
            "A1.1,Halobacteriota\n"
            "A2.1,Halobacteriota\n"
            "A3.1,Thermoplasmatota\n"
        )
        self.expected = {
            "Firmicutes_A": {"B1.1", "B2.1"},
            "NA": {"B3.1"},
            "Single": {"B4.1"},
            "Halobacteriota": {"A1.1", "A2.1"},
            "Thermoplasmatota": {"A3.1"},
        }

    def run_script(self, *arguments, alternate_flag="--alternate_input", error=None):
        command = [
            RSCRIPT, str(REPOSITORY / "workflow/scripts/make_tree.R"),
            "--input", str(self.primary), "--taxonomy", str(self.taxonomy),
            "--output", str(self.output),
        ]
        if alternate_flag:
            command.extend([alternate_flag, str(self.alternate)])
        result = subprocess.run([*command, *arguments], capture_output=True, text=True, timeout=60)
        details = result.stdout + result.stderr
        if error:
            self.assertNotEqual(result.returncode, 0, details)
            self.assertIn(error, details)
            self.assertFalse(self.output.exists())
        else:
            self.assertEqual(result.returncode, 0, details)
        return details

    def assert_output(self, expected=None):
        inspection = subprocess.run(
            [RSCRIPT, "-e", '''
                library(ape)
                trees <- readRDS(commandArgs(TRUE)[1])
                stopifnot(is.list(trees), !anyDuplicated(names(trees)))
                tips <- lapply(names(trees), function(phylum) {
                  tree <- trees[[phylum]]
                  stopifnot(inherits(tree, "phylo"), !anyDuplicated(tree$tip.label),
                            all(is.finite(tree$edge.length)), all(tree$edge.length >= 0))
                  if (length(tree$tip.label) > 1) stopifnot(is.ultrametric(tree))
                  data.frame(phylum=phylum, tip=tree$tip.label)
                })
                write.table(do.call(rbind, tips), stdout(), sep="\\t", quote=FALSE, row.names=FALSE)
            ''', str(self.output)],
            capture_output=True, text=True, timeout=30,
        )
        self.assertEqual(inspection.returncode, 0, inspection.stdout + inspection.stderr)
        actual = {}
        for row in csv.DictReader(io.StringIO(inspection.stdout), delimiter="\t"):
            actual.setdefault(row["phylum"], set()).add(row["tip"])
        self.assertEqual(actual, self.expected if expected is None else expected)

    def test_bacterial_and_archaeal_subtrees_use_canonical_ids(self):
        self.run_script()
        self.assert_output()

    def test_workflow_flags_keep_phylum_names_from_taxonomy(self):
        self.run_script("-n", "Archaea", alternate_flag="-a")
        self.assert_output()

    def test_legacy_long_alternate_flag_is_supported(self):
        self.run_script(alternate_flag="--altinput")
        self.assert_output()

    def test_primary_tree_can_be_used_alone(self):
        self.taxonomy.write_text("\n".join(
            line for line in self.taxonomy.read_text().splitlines() if not line.startswith("A")
        ) + "\n")
        self.run_script(alternate_flag=None)
        self.assert_output({name: tips for name, tips in self.expected.items() if not name.endswith("ota")})

    def test_omitted_archaeal_tree_is_rejected(self):
        self.run_script(alternate_flag=None, error="Unexpected species missing from input trees")

    def test_missing_species_is_rejected(self):
        self.primary.write_text(self.primary.read_text().replace("B3:", "unknown:"))
        self.run_script(error="B3.1")

    def test_ambiguous_taxonomy_versions_are_rejected(self):
        self.taxonomy.write_text(self.taxonomy.read_text() + "B1.2,Firmicutes_A\n")
        self.run_script(error="Taxonomy has missing or ambiguous species IDs")

    def test_ambiguous_normalized_tips_are_rejected(self):
        self.primary.write_text(self.primary.read_text().replace("B_out:", "B1.2:"))
        self.run_script(error="ambiguous after species-ID normalization")

    def test_overlapping_input_trees_are_rejected(self):
        self.alternate.write_text(self.primary.read_text())
        self.run_script(error="Species occur in more than one input tree")

    def test_phylum_split_between_input_trees_is_rejected(self):
        self.taxonomy.write_text(self.taxonomy.read_text().replace("B1.1,Firmicutes_A", "B1.1,Shared").replace("A1.1,Halobacteriota", "A1.1,Shared"))
        self.run_script(error="Phylum must be covered by exactly one input tree: Shared")

    def test_exact_audited_source_gaps_are_supported(self):
        self.primary.write_text(self.primary.read_text().replace("B3:", "unknown:"))
        exclusions = self.directory / "source exclusions.csv"
        exclusions.write_text("cluster,reason\nB3.1,absent_from_source_tree\n")
        details = self.run_script("--allow_missing_species", str(exclusions))
        self.assertIn("1 audited source species have no tree tips", details)
        self.assert_output({name: tips for name, tips in self.expected.items() if name != "NA"})

    def test_exclusion_list_must_match_actual_gaps(self):
        exclusions = self.directory / "source exclusions.csv"
        exclusions.write_text("cluster\nB3.1\n")
        self.run_script("--allow_missing_species", str(exclusions), error="exact audited --allow_missing_species list")

    def test_negative_branch_lengths_are_rejected(self):
        self.primary.write_text(self.primary.read_text().replace("B1:0.2", "B1:-0.2"))
        self.run_script(error="finite, nonnegative branch lengths")

    def test_multiple_trees_in_one_file_are_rejected(self):
        self.primary.write_text(self.primary.read_text() * 2)
        self.run_script(error="exactly one phylogenetic tree")


if __name__ == "__main__":
    unittest.main()
