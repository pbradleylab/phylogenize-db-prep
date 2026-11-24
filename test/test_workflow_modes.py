"""Check mode flags, final targets, and Snakemake's annotation-only graph."""

import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
from types import SimpleNamespace
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("workflow_utils", REPOSITORY / "workflow/scripts/utils.py")
utils = importlib.util.module_from_spec(spec)
spec.loader.exec_module(utils)


class ModeFlagTests(unittest.TestCase):
    def test_native_and_legacy_booleans(self):
        for value in [True, "True", "true", " TRUE "]:
            self.assertTrue(utils.get_config_bool(value, "flag"))
        for value in [False, "False", "false", " FALSE "]:
            self.assertFalse(utils.get_config_bool(value, "flag"))

    def test_invalid_flags_are_rejected(self):
        for value in ["yes", "", 1, None]:
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "flag must be true or false"):
                utils.get_config_bool(value, "flag")

    def test_legacy_annotation_flag_is_supported(self):
        self.assertTrue(utils.get_annotate_only({"annotation": {"annotate_only": "True"}}))

    def test_modern_flag_takes_precedence(self):
        self.assertFalse(utils.get_annotate_only({
            "mode": {"annotate_only": "False"}, "annotation": {"annotate_only": "True"},
        }))

    def test_missing_flags_use_the_requested_default(self):
        self.assertFalse(utils.get_annotate_only({}))
        self.assertIsNone(utils.get_annotate_only({}, default=None))


class FinalTargetTests(unittest.TestCase):
    def targets(self, annotate_only=False, options_16s=None):
        patterns = {
            "get_binary": "results/{database}/binary/get_binary/{mapping_db}-binary.rds",
            "get_tree": "results/{database}/binary/get_tree/{mapping_db}-tree.rds",
            "run_annotation": "results/{database}/annotation/run_anntoation/{mapping_db}.done",
            "link_go_nodes": "results/{database}/annotation/link_go_nodes/{mapping_db}/node_annotations.tsv",
            "all_16S": "results/{database}/16S/.{mapping_db}.complete",
        }
        rules = SimpleNamespace(**{
            name: SimpleNamespace(output=[pattern]) for name, pattern in patterns.items()
        })
        rules.get_continuous = SimpleNamespace(output=SimpleNamespace(
            rds="results/{database}/binary/get_binary/get_continuous/{mapping_db}-continuous.rds"
        ))
        config = {
            "database": "mode-test", "files": {"fasta": {"gut": "gut.faa", "skin": "skin.faa"}},
            "target_db": {"paths": {"reference": "reference.faa"}},
        }
        if options_16s is not None:
            config["files"]["16S"] = options_16s
        namespace = {
            "config": config, "rules": rules, "annotate_only": annotate_only,
            "get_config_bool": utils.get_config_bool,
        }
        source = (REPOSITORY / "workflow/Snakefile").read_text()
        function = source[source.index("def get_final_output("):source.index("rule all:")]
        exec(function, namespace)
        return namespace["get_final_output"](None)

    def test_full_mode_requests_profiles_ko_and_go(self):
        outputs = self.targets()
        self.assertEqual(len(outputs), 10)
        for dataset in ["gut", "skin"]:
            self.assertIn(f"results/mode-test/binary/get_binary/{dataset}-binary.rds", outputs)
            self.assertIn(f"results/mode-test/binary/get_tree/{dataset}-tree.rds", outputs)
            self.assertIn(f"results/mode-test/binary/get_binary/get_continuous/{dataset}-continuous.rds", outputs)
            self.assertIn(f"results/mode-test/annotation/run_anntoation/{dataset}.done", outputs)
            self.assertIn(f"results/mode-test/annotation/link_go_nodes/{dataset}/node_annotations.tsv", outputs)

    def test_16s_targets_are_selected_per_dataset(self):
        outputs = self.targets(options_16s={
            "gut": {"generate_for_16S": True}, "skin": {"generate_for_16S": "False"},
        })
        self.assertIn("results/mode-test/16S/.gut.complete", outputs)
        self.assertNotIn("results/mode-test/16S/.skin.complete", outputs)

    def test_annotation_only_excludes_profiles_and_16s(self):
        outputs = self.targets(annotate_only=True, options_16s={"gut": {"generate_for_16S": "True"}})
        self.assertEqual(len(outputs), 4)
        self.assertTrue(all("/annotation/" in output for output in outputs))


@unittest.skipUnless(importlib.util.find_spec("snakemake"), "Snakemake is required")
class ModeGraphTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="workflow modes ")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        (self.directory / "workflow").symlink_to(REPOSITORY / "workflow", target_is_directory=True)
        (self.directory / "config").mkdir()
        (self.directory / "config/tools.json").write_text((REPOSITORY / "config/tools.json").read_text())
        self.inputs = self.directory / "inputs"
        self.inputs.mkdir()
        for name in ["catalogue.faa", "reference.faa"]:
            (self.inputs / name).write_text(">F1\nM" + "A" * 60 + "\n")
        for name in ["members.tsv", "metadata.tsv", "tree.nwk"]:
            (self.inputs / name).touch()
        (self.inputs / "genomes").mkdir()
        self.config = {
            "database": "mode-test", "indir": "", "mode": {"annotate_only": False},
            "files": {
                "fasta": {"catalogue": str(self.inputs / "catalogue.faa")},
                "mapping": {"catalogue": str(self.inputs / "members.tsv")},
                "taxonomy": {"catalogue": str(self.inputs / "metadata.tsv")},
                "tree": {"catalogue": str(self.inputs / "tree.nwk")},
                "16S": {"catalogue": {"generate_for_16S": False, "fasta_dir": str(self.inputs / "genomes")}},
            },
            "target_db": {"paths": {"reference": str(self.inputs / "reference.faa")}},
        }

    def cached_mappings(self):
        paths = {
            "results/mode-test/binary/combined_species_hits/catalogue_mode-test.tsv": "query\ttarget\nG1_1\tF1\nG1_2\tF2\n",
            "results/mode-test/binary/get_taxonomy/catalogue-binary-temp.csv": "cluster,target\nG1,F1\nG1,F2\n",
            "results/mode-test/clustering/mmseqs2_linclust/catalogue/unaligned_linclust_rep_seq.fasta": ">F2\nMAAAA\n",
            "results/mode-test/clustering/mmseqs2/convertalis_full_blast/reference_catalogue_convertlis.8": "G1_1\t60\tF1\t60\t0\t100\t60\t60\n",
        }
        for name, content in paths.items():
            path = self.directory / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content)

    def run_workflow(self, *options, error=None, target="all", dry_run=True, environment=None):
        (self.directory / "config/config.json").write_text(json.dumps(self.config))
        result = subprocess.run([
            sys.executable, "-m", "snakemake", target, "--snakefile", "workflow/Snakefile",
            "--directory", str(self.directory), "--nolock", "--cores", "1",
            *(["--dry-run"] if dry_run else []), *options,
        ], cwd=self.directory, env={
            **os.environ, "PYTHONDONTWRITEBYTECODE": "1", "XDG_CACHE_HOME": str(self.directory / "cache"),
            **(environment or {}),
        }, capture_output=True, text=True, timeout=90)
        output = result.stdout + result.stderr
        if error:
            self.assertNotEqual(result.returncode, 0, output)
            self.assertIn(error, output)
        else:
            self.assertEqual(result.returncode, 0, output)
        return output

    eggnog_files = {
        "eggnog_proteins.dmnd": "download_diamond",
        "eggnog.db": "download_eggnog_annotations",
        "eggnog.taxa.db": "download_eggnog_taxa",
        "eggnog.taxa.db.traverse.pkl": "download_eggnog_taxa",
    }
    emapper_target = "results/mode-test/annotation/emapper/catalogue/catalogue.part_001/go_annotations.emapper.annotations"

    def cached_annotation_chunks(self):
        self.cached_mappings()
        base = "results/mode-test/annotation"
        paths = {
            f"{base}/get_fasta_seqs_to_annotate/catalogue.fasta": ">F1\nMAAAA\n>F2\nMAAAA\n",
            f"{base}/get_fasta_seqs_to_annotate/catalogue.txt": "F1\n",
            f"{base}/anvio_reformat_fasta/catalogue.fasta": ">c_000000000001\nMAAAA\n>c_000000000002\nMAAAA\n",
            f"{base}/anvio_reformat_fasta/catalogue.names.tsv": "c_000000000001\tF1\nc_000000000002\tF2\n",
        }
        for index in [1, 2]:
            chunk = f"catalogue.part_{index:03}"
            paths[f"{base}/break_fasta_apart/catalogue/{chunk}.fasta"] = f">c_{index:012}\nMAAAA\n"
            paths[f"{base}/get_headers/catalogue/{chunk}_header.txt"] = f"gene_callers_id\tlinker_info\n0\tF{index}\n"
            paths[f"{base}/get_headers/catalogue/{chunk}.tmp"] = f"F{index}\n"
        for name, content in paths.items():
            path = self.directory / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(content)

    def test_emapper_downloads_databases_in_an_existing_empty_directory(self):
        self.cached_annotation_chunks()
        data_dir = self.directory / "resources/eggnog"
        data_dir.mkdir(parents=True)
        (data_dir / ".snakemake_timestamp").touch()
        (data_dir / "eggnog_proteins.dmnd.gz").touch()
        output = self.run_workflow("--printshellcmds", target=self.emapper_target)
        for rule in set(self.eggnog_files.values()):
            self.assertIn("rule " + rule + ":", output)
        self.assertIn("--data_dir resources/eggnog", output)
        self.assertNotIn("--data_dir resources/eggnog/eggnog_proteins.dmnd", output)

    def test_emapper_downloads_only_the_missing_database(self):
        self.cached_annotation_chunks()
        data_dir = self.directory / "resources/eggnog"
        data_dir.mkdir(parents=True)
        for filename in self.eggnog_files:
            (data_dir / filename).write_text("cached database")
        for filename, expected_rule in self.eggnog_files.items():
            with self.subTest(missing=filename):
                (data_dir / filename).unlink()
                output = self.run_workflow(target=self.emapper_target)
                for rule in set(self.eggnog_files.values()):
                    if rule == expected_rule:
                        self.assertIn("rule " + rule + ":", output)
                    else:
                        self.assertNotIn("rule " + rule + ":", output)
                (data_dir / filename).write_text("cached database")

    def test_emapper_reuses_complete_eggnog_databases(self):
        self.cached_annotation_chunks()
        data_dir = self.directory / "resources/eggnog"
        data_dir.mkdir(parents=True)
        for filename in self.eggnog_files:
            (data_dir / filename).write_text("cached database")
        output = self.run_workflow(target=self.emapper_target)
        self.assertIn("rule emapper:", output)
        for rule in set(self.eggnog_files.values()):
            self.assertNotIn("rule " + rule + ":", output)

    def test_emapper_runs_per_ko_chunk_and_combines_original_family_ids(self):
        self.config["mode"]["annotate_only"] = True
        self.cached_annotation_chunks()
        data_dir = self.directory / "resources/eggnog"
        data_dir.mkdir(parents=True)
        for filename in self.eggnog_files:
            (data_dir / filename).write_text("cached database")
        bin_dir = self.directory / "bin"
        bin_dir.mkdir()
        emapper = bin_dir / "emapper.py"
        emapper.write_text(f"#!{sys.executable}\n" + '''\
import argparse
from pathlib import Path

parser = argparse.ArgumentParser()
for option in ["-i", "-o", "--output_dir", "--cpu", "--go_evidence", "--data_dir"]:
    parser.add_argument(option, required=True)
parser.add_argument("--override", action="store_true", required=True)
args = parser.parse_args()
with open("emapper_calls.txt", "a") as calls:
    calls.write(args.i + "\\n")
header = ["#query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl", "COG_category", "Description", "Preferred_name", "GOs"]
with open(Path(args.output_dir) / (args.o + ".emapper.annotations"), "w") as output:
    output.write("\\t".join(header) + "\\n")
    for line in Path(args.i).read_text().splitlines():
        if line.startswith(">"):
            row = [line[1:].split()[0], "seed", "1e-20", "100", "OG", "2", "S", "description", "name", "GO:0008150"]
            output.write("\\t".join(row) + "\\n")
''')
        emapper.chmod(0o755)
        environment = {"PATH": str(bin_dir) + os.pathsep + os.environ["PATH"]}
        target = "results/mode-test/annotation/make_go_db/catalogue.3"
        output = self.run_workflow(target=target, dry_run=False, environment=environment)
        self.assertIn("rule combine_emapper:", output)
        expected_inputs = [
            f"results/mode-test/annotation/break_fasta_apart/catalogue/catalogue.part_{index:03}.fasta"
            for index in [1, 2]
        ]
        self.assertCountEqual((self.directory / "emapper_calls.txt").read_text().splitlines(), expected_inputs)
        self.assertEqual((self.directory / target).read_text(), "F1\tGO:0008150\nF2\tGO:0008150\n")
        combined = self.directory / "results/mode-test/annotation/combine_emapper/catalogue/go_annotations.emapper.annotations"
        rows = combined.read_text().splitlines()
        self.assertEqual(sum(row.startswith("#query\t") for row in rows), 1)
        self.assertFalse(any(row.startswith("c_") for row in rows))
        self.run_workflow(target=target, dry_run=False, environment=environment)
        self.assertEqual(len((self.directory / "emapper_calls.txt").read_text().splitlines()), 2)

    def test_full_mode_builds_profiles_and_both_annotation_branches(self):
        output = self.run_workflow()
        for name in ["get_binary", "get_tree", "get_continuous", "run_annotation", "link_go_nodes"]:
            self.assertIn("rule " + name + ":", output)
        self.assertNotIn("rule all_16S:", output)
        self.assertNotIn("Would remove temporary output results/mode-test/binary/get_taxonomy/catalogue-binary-temp.csv", output)

    def test_annotation_only_reuses_mappings_without_search_or_profiles(self):
        self.config["mode"]["annotate_only"] = "True"
        self.config["files"]["16S"]["catalogue"]["generate_for_16S"] = True
        self.cached_mappings()
        output = self.run_workflow()
        for name in ["map_query", "mmseqs2_linclust", "get_taxonomy", "get_binary", "get_tree", "get_continuous", "all_16S"]:
            self.assertNotIn("rule " + name + ":", output)
        self.assertIn("rule run_annotation:", output)
        self.assertIn("rule link_go_nodes:", output)

    def test_annotation_only_rejects_missing_intermediates(self):
        self.config["mode"]["annotate_only"] = True
        self.run_workflow(error="annotate_only requires an existing mapping intermediate")

    def test_legacy_command_line_flag_overrides_default_modern_flag(self):
        self.cached_mappings()
        override = self.directory / "legacy.json"
        override.write_text(json.dumps({"annotation": {"annotate_only": "True"}}))
        output = self.run_workflow("--configfile", str(override))
        self.assertIn("rule run_annotation:", output)
        self.assertNotIn("rule map_query:", output)
        self.assertNotIn("rule get_binary:", output)

    def test_modern_false_overrides_legacy_default(self):
        self.config.pop("mode")
        self.config["annotation"] = {"annotate_only": "True"}
        override = self.directory / "modern.json"
        override.write_text(json.dumps({"mode": {"annotate_only": False}}))
        output = self.run_workflow("--configfile", str(override))
        self.assertIn("rule map_query:", output)
        self.assertIn("rule get_binary:", output)

    def test_enabled_16s_schedules_final_alignment_and_tree(self):
        self.config["files"]["16S"]["catalogue"]["generate_for_16S"] = "True"
        output = self.run_workflow()
        for name in ["all_16S", "fix_aln_file", "prune_16S_tree"]:
            self.assertIn("rule " + name + ":", output)
        self.assertIn("results/mode-test/16S/.catalogue.complete", output)


if __name__ == "__main__":
    unittest.main()
