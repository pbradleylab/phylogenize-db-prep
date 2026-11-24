"""Validate Linclust configuration and exercise the real rule on small proteins."""

import importlib.util
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("workflow_utils", REPOSITORY / "workflow/scripts/utils.py")
utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(utils)


class LinclustSettingsTests(unittest.TestCase):
    def settings(self, **options):
        return utils.get_linclust_settings({"mmseqs2": {"linclust": options}})

    def test_defaults_preserve_existing_clustering_criteria(self):
        self.assertEqual(self.settings(), {
            "min_seq_id": 0.0, "coverage": 0.8, "cov_mode": 0, "evalue": 0.001,
        })

    def test_existing_identity_key_and_new_options_are_honored(self):
        self.assertEqual(self.settings(seq_id_precent=0.5, coverage=0.9, cov_mode=1, evalue=1e-8), {
            "min_seq_id": 0.5, "coverage": 0.9, "cov_mode": 1, "evalue": 1e-8,
        })

    def test_invalid_thresholds_are_rejected(self):
        for key in ["seq_id_precent", "coverage", "evalue"]:
            invalid = [-0.1, float("inf"), float("nan"), True, "0.5", None]
            if key != "evalue":
                invalid.append(1.1)
            for value in invalid:
                with self.subTest(key=key, value=value), self.assertRaisesRegex(ValueError, key):
                    self.settings(**{key: value})

    def test_coverage_mode_requires_a_supported_integer(self):
        for value in [-1, 6, 0.5, 1.0, True, "0", None]:
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "cov_mode"):
                self.settings(cov_mode=value)


@unittest.skipUnless(importlib.util.find_spec("snakemake"), "Snakemake is required")
class LinclustRuleTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="linclust settings ")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        (self.directory / "workflow").symlink_to(REPOSITORY / "workflow", target_is_directory=True)
        (self.directory / "config").mkdir()
        self.tools = json.loads((REPOSITORY / "config/tools.json").read_text())
        (self.directory / "config/config.json").write_text(json.dumps({
            "database": "linclust-test", "mode": {"annotate_only": False},
            "files": {"fasta": {"catalogue": "unused.faa"}},
            "target_db": {"paths": {"reference": "unused-reference.faa"}},
        }))
        self.input = self.directory / "results/linclust-test/clustering/faSomeRecords/final_unmapped/catalogue_clustered.fa"
        self.input.parent.mkdir(parents=True)
        self.input.write_text(">protein\nM" + "A" * 100 + "\n")
        self.target = "results/linclust-test/clustering/mmseqs2_linclust/catalogue/unaligned_linclust_cluster.tsv"

    def run_rule(self, dry_run=True, **options):
        self.tools["mmseqs2"]["linclust"].update(options)
        self.tools["mmseqs2"]["linclust"]["tmp_dir"] = str(self.directory / "scratch with spaces")
        (self.directory / "config/tools.json").write_text(json.dumps(self.tools))
        result = subprocess.run([
            sys.executable, "-m", "snakemake", self.target,
            "--snakefile", "workflow/Snakefile", "--cores", "1", "--nolock",
            "--printshellcmds", "--allowed-rules", "mmseqs2_linclust",
            *(["--dry-run"] if dry_run else []),
        ], cwd=self.directory, env={
            **os.environ, "PYTHONDONTWRITEBYTECODE": "1", "XDG_CACHE_HOME": str(self.directory / "cache"),
        }, capture_output=True, text=True, timeout=90)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result.stdout + result.stderr

    def test_rule_passes_configured_criteria_to_mmseqs(self):
        output = self.run_rule(seq_id_precent=0.5, coverage=0.9, cov_mode=1, evalue=1e-8)
        for argument in ["--min-seq-id 0.5", "-c 0.9", "--cov-mode 1", "-e 1e-08"]:
            self.assertIn(argument, output)

    @unittest.skipUnless(shutil.which("mmseqs"), "MMseqs2 is required")
    def test_threshold_changes_affect_actual_families(self):
        generator = random.Random(37)
        alphabet = "ACDEFGHIKLMNPQRSTVWY"
        sequence = "".join(generator.choice(alphabet) for _ in range(300))
        variant = list(sequence)
        for index in generator.sample(range(300), 75):
            variant[index] = generator.choice(alphabet.replace(sequence[index], ""))
        self.input.write_text(f">full\n{sequence}\n>variant\n{''.join(variant)}\n")
        self.run_rule(dry_run=False, seq_id_precent=0.5, coverage=0.8, cov_mode=0)
        self.assertEqual(self.family_count(), 1)
        self.run_rule(dry_run=False, seq_id_precent=0.95)
        self.assertEqual(self.family_count(), 2)

        self.input.write_text(f">full\n{sequence}\n>fragment\n{sequence[:100]}\n")
        self.run_rule(dry_run=False, seq_id_precent=0.5, coverage=0.8)
        self.assertEqual(self.family_count(), 2)
        self.run_rule(dry_run=False, coverage=0.25)
        self.assertEqual(self.family_count(), 1)

    def family_count(self):
        rows = (self.directory / self.target).read_text().splitlines()
        self.assertEqual(len(rows), 2)
        return len({row.split("\t")[0] for row in rows})


if __name__ == "__main__":
    unittest.main()
