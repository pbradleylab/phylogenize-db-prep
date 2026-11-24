"""Check real annotation handoffs without downloading or searching databases."""

import csv
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
SCRIPTS = REPOSITORY / "workflow/scripts"
sys.path.insert(0, str(SCRIPTS))
from annotation_io import fasta_records
from combine_ko_annotations import combine_annotations
from make_annotation_headers import make_headers
from make_go_db import make_go_db
from merge_annotations import KO_COLUMNS, merge_annotations
from select_annotation_sequences import select_sequences

GO_PYTHON = os.environ.get("PHYLOGENIZE_GO_PYTHON") or (sys.executable if importlib.util.find_spec("goatools") else None)
MATRIX_PYTHON = os.environ.get("PHYLOGENIZE_MATRIX_PYTHON") or (sys.executable if importlib.util.find_spec("polars") else None)
OBO = '''format-version: 1.2
data-version: test-release

[Term]
id: GO:0000001
name: example process
namespace: biological_process
alt_id: GO:0000004
def: "An example process." []

[Term]
id: GO:0000002
name: example activity
namespace: molecular_function

[Term]
id: GO:0000003
name: example compartment
namespace: cellular_component

[Term]
id: GO:0000005
name: obsolete example
namespace: biological_process
is_obsolete: true
'''
EMAPPER_HEADER = "#query\tGOs\tseed_ortholog\tevalue\tscore\n"


class AnnotationHandoffTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="annotation handoffs '")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.hits = self.write("hits.tsv", "query\ttarget\nq1\t0002\nq2\tdenovo-family\nq3\tNA\nq4\t0002\n")
        self.pairs = self.write("pairs.csv", "cluster,target\nS1.1,0002\nS2.1,0002\nS2.1,denovo-family\nS3.1,NA\n")
        self.obo = self.write("go.obo", OBO)

    def write(self, name, content):
        path = self.directory / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
        return path

    def rows(self, path, delimiter="\t"):
        with path.open(newline="") as handle:
            return list(csv.DictReader(handle, delimiter=delimiter))

    def run_script(self, interpreter, script, *arguments, error=None):
        result = subprocess.run([interpreter, str(SCRIPTS / script), *map(str, arguments)],
                                capture_output=True, text=True, timeout=60,
                                env={**os.environ, "PYTHONDONTWRITEBYTECODE": "1"})
        if error:
            self.assertNotEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn(error, result.stderr)
        else:
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        return result

    def test_representatives_are_complete_unique_and_retry_safe(self):
        reference = self.write("reference.faa", ">0002 description\nMAAA\n>NA\nMCCC\n>unused\nMDDD\n")
        duplicate = self.write("other.faa", ">0002\nMAAA\n")
        denovo = self.write("denovo.faa", ">denovo-family\nMEEE\n")
        output, ids = self.directory / "selected.faa", self.directory / "ids.txt"
        for _ in range(2):
            select_sequences(self.hits, [reference, duplicate], denovo, output, ids)
            self.assertEqual([identifier for identifier, _ in fasta_records(output)], ["0002", "NA", "denovo-family"])
            self.assertEqual(ids.read_text(), "0002\nNA\ndenovo-family\n")

    def test_conflicting_or_missing_representatives_are_rejected(self):
        reference = self.write("reference.faa", ">0002\nMAAA\n>NA\nMCCC\n")
        denovo = self.write("denovo.faa", ">0002\nMEEE\n")
        with self.assertRaisesRegex(ValueError, "Conflicting"):
            select_sequences(self.hits, [reference], denovo, self.directory / "out", self.directory / "ids")
        denovo.write_text("")
        with self.assertRaisesRegex(ValueError, "Missing representatives"):
            select_sequences(self.hits, [reference], denovo, self.directory / "out", self.directory / "ids")

    def test_chunk_headers_use_report_ids_and_local_gene_numbers(self):
        fasta = self.write("chunk.faa", ">arbitrary_B\nMAAA\n>arbitrary_A\nMCCC\n")
        report = self.write("report.tsv", "arbitrary_A\tNA description\narbitrary_B\t0002 description\n")
        output = self.directory / "headers.tsv"
        make_headers(fasta, report, output, self.directory / "ids")
        self.assertEqual(self.rows(output), [
            {"gene_callers_id": "0", "linker_info": "0002"},
            {"gene_callers_id": "1", "linker_info": "NA"},
        ])
        report.write_text("arbitrary_A\tNA\n")
        with self.assertRaisesRegex(ValueError, "missing from rename report"):
            make_headers(fasta, report, output, self.directory / "ids")

    def prepare_ko(self):
        fasta = self.write("chunk.faa", ">first\nMAAA\n>second\nMCCC\n")
        calls = self.directory / "calls.tsv"
        result = subprocess.run(["bash", str(SCRIPTS / "make_external_gene_calls.sh"), str(fasta), str(calls)],
                                capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertEqual([row["gene_callers_id"] for row in self.rows(calls)], ["0", "1"])
        linker = self.write("headers.tsv", "gene_callers_id\tlinker_info\n0\t0002\n1\tNA\n")
        functions = self.write("functions.tsv", "gene_callers_id\tsource\taccession\tfunction\te_value\n"
                               "0\tKOfam\tK00001\tfirst, function\t1e-20\n"
                               "0\tKOfam\tK00002\tsecond function\t1e-15\n"
                               "1\tOther\tother\tother function\t0\n"
                               "1\tKOfam\tK00003\tthird function\t1e-10\n")
        return functions, calls, linker

    def test_ko_merge_preserves_multiple_annotations_and_literal_ids(self):
        functions, calls, linker = self.prepare_ko()
        output = self.directory / "merged.tsv"
        merge_annotations(functions, calls, linker, output)
        self.assertEqual([(row["linker_info"], row["accession"]) for row in self.rows(output)],
                         [("0002", "K00001"), ("0002", "K00002"), ("NA", "K00003")])

    def test_ko_merge_rejects_missing_gene_links(self):
        functions, calls, linker = self.prepare_ko()
        linker.write_text("gene_callers_id\tlinker_info\n0\t0002\n")
        with self.assertRaisesRegex(ValueError, "different gene IDs"):
            merge_annotations(functions, calls, linker, self.directory / "out")

    def test_ko_chunks_have_one_header_and_allow_empty_annotations(self):
        first = self.write("first.tsv", "\t".join(KO_COLUMNS) + "\n0\t0002\tK00001\tfirst\n")
        empty = self.write("empty.tsv", "\t".join(KO_COLUMNS) + "\n")
        last = self.write("last.tsv", "\t".join(KO_COLUMNS) + "\n0\tNA\tK00002\tlast\n")
        output = self.directory / "combined.tsv"
        combine_annotations([first, empty, last], output)
        self.assertEqual(len(self.rows(output)), 2)
        self.assertEqual(output.read_text().count("gene_callers_id"), 1)
        with self.assertRaisesRegex(ValueError, "more than one KO chunk"):
            combine_annotations([first, first], output)
        combine_annotations([empty], output)
        self.assertEqual(self.rows(output), [])

    def test_go_lookup_uses_named_column(self):
        annotations = self.write("emapper.tsv", EMAPPER_HEADER + "0002\tGO:0000001\tseed\t1e-20\t100\n"
                                 "family-with-hyphens\tGO:0000002\tseed\t1e-10\t50\nNA\t-\tseed\t1e-10\t50\n")
        output = self.directory / "lookup"
        make_go_db(annotations, output)
        self.assertEqual(output.read_text(), "0002\tGO:0000001\nfamily-with-hyphens\tGO:0000002\n")

    @unittest.skipUnless(GO_PYTHON, "Set PHYLOGENIZE_GO_PYTHON to an environment with GOATOOLS")
    def test_go_triples_keep_query_scores_and_ontology_status(self):
        triples = self.write("triples.tsv", "0002\tGO:0000004\t0.700\nNA\tGO:9999999\t0.200\n"
                             "denovo-family\tGO:0000005\t0.100\n")
        output = self.directory / "descriptions.tsv"
        self.run_script(GO_PYTHON, "annotate_go_terms.py", "--obo", self.obo, "--go_ids", triples, "--out", output)
        rows = self.rows(output)
        self.assertEqual([(row["query"], row["score"], row["status"]) for row in rows],
                         [("0002", "0.700", "active"), ("NA", "0.200", "unknown"), ("denovo-family", "0.100", "obsolete")])
        self.assertEqual(rows[0]["go_id"], "GO:0000001")
        self.assertEqual(rows[0]["original_go_id"], "GO:0000004")
        triples.write_text("")
        self.run_script(GO_PYTHON, "annotate_go_terms.py", "--obo", self.obo, "--go_ids", triples, "--out", output)
        self.assertEqual(self.rows(output), [])

    @unittest.skipUnless(GO_PYTHON, "Set PHYLOGENIZE_GO_PYTHON to an environment with GOATOOLS")
    def test_direct_go_export_includes_denovo_and_separates_invalid_terms(self):
        annotations = self.write("emapper.tsv", EMAPPER_HEADER +
                                 "0002\tGO:0000001,GO:0000002,GO:0000003,GO:0000004\tseed1\t1e-20\t100\n"
                                 "denovo-family\tGO:0000002\tseed2\t1e-15\t80\n"
                                 "NA\tGO:0000005,GO:9999999\tseed3\t1e-10\t50\n")
        families, evidence = self.directory / "families.tsv", self.directory / "evidence.tsv"
        self.run_script(GO_PYTHON, "export_go_annotations.py", "--annotations", annotations, "--hits", self.hits,
                        "--obo", self.obo, "--families", families, "--evidence", evidence, "--go-evidence", "non-electronic")
        rows = self.rows(families)
        self.assertEqual(len(rows), 4)
        self.assertEqual({row["aspect"] for row in rows}, {"P", "F", "C"})
        self.assertIn(("denovo-family", "GO:0000002"), {(row["node_head"], row["accession"]) for row in rows})
        self.assertEqual({row["status"] for row in self.rows(evidence)}, {"active", "obsolete", "unknown"})
        self.assertEqual(self.rows(evidence)[0]["seed_score"], "100")

    @unittest.skipUnless(MATRIX_PYTHON, "Set PHYLOGENIZE_MATRIX_PYTHON to an environment with Polars")
    def test_species_links_match_matrix_pairs_without_duplicate_centroid_expansion(self):
        annotations = self.write("families.tsv", "node_head\taccession\tfunction\n0002\tGO:0000001\tfirst, function\n"
                                 "denovo-family\tGO:0000002\tsecond function\nNA\tGO:0000003\tthird function\n")
        output = self.directory / "species.csv"
        self.run_script(MATRIX_PYTHON, "sequence_to_annotation_linker.py", self.hits, annotations, self.pairs, output)
        self.assertEqual({(row["linker_info"], row["node_head"]) for row in self.rows(output, ",")},
                         {("S1.1", "0002"), ("S2.1", "0002"), ("S2.1", "denovo-family"), ("S3.1", "NA")})
        annotations.write_text("node_head\taccession\tfunction\n")
        self.run_script(MATRIX_PYTHON, "sequence_to_annotation_linker.py", self.hits, annotations, self.pairs, output)
        self.assertEqual(self.rows(output, ","), [])
        annotations.write_text("node_head\taccession\tfunction\nUNKNOWN\tGO:0000001\tfunction\n")
        self.run_script(MATRIX_PYTHON, "sequence_to_annotation_linker.py", self.hits, annotations, self.pairs, output,
                        error="families absent from combined hits")


@unittest.skipUnless(importlib.util.find_spec("snakemake") and GO_PYTHON and MATRIX_PYTHON,
                     "Snakemake, GOATOOLS, and Polars interpreters are required")
class AnnotationWorkflowTests(unittest.TestCase):
    write = AnnotationHandoffTests.write
    rows = AnnotationHandoffTests.rows

    def setUp(self):
        AnnotationHandoffTests.setUp(self)
        (self.directory / "workflow").symlink_to(REPOSITORY / "workflow", target_is_directory=True)
        reference = self.write("reference.faa", ">0002\nMAAA\n>NA\nMCCC\n")
        self.write("config/config.json", json.dumps({
            "database": "handoff-test", "mode": {"annotate_only": True},
            "files": {"fasta": {"catalogue": "unused.faa"}},
            "target_db": {"paths": {"reference": str(reference)}},
        }))
        self.write("config/tools.json", (REPOSITORY / "config/tools.json").read_text())
        self.write("results/handoff-test/binary/combined_species_hits/catalogue_handoff-test.tsv", self.hits.read_text())
        self.write("results/handoff-test/binary/get_taxonomy/catalogue-binary-temp.csv", self.pairs.read_text())
        self.write("results/handoff-test/clustering/mmseqs2_linclust/catalogue/unaligned_linclust_rep_seq.fasta",
                   ">denovo-family\nMEEE\n")
        self.write("results/handoff-test/clustering/mmseqs2/convertalis_full_blast/reference_catalogue_convertlis.8",
                   "q1\t4\t0002\t4\t1e-20\t100\t4\t4\n")
        self.base = "results/handoff-test/annotation"
        self.write(f"{self.base}/get_fasta_seqs_to_annotate/catalogue.fasta", ">0002\nMAAA\n>NA\nMCCC\n>denovo-family\nMEEE\n")
        self.write(f"{self.base}/get_fasta_seqs_to_annotate/catalogue.txt", "0002\nNA\ndenovo-family\n")
        self.write(f"{self.base}/anvio_reformat_fasta/catalogue.fasta",
                   ">c_000000000001\nMAAA\n>c_000000000002\nMCCC\n>c_000000000003\nMEEE\n")
        self.write(f"{self.base}/anvio_reformat_fasta/catalogue.names.tsv",
                   "c_000000000001\t0002\nc_000000000002\tNA\nc_000000000003\tdenovo-family\n")
        self.write("resources/anvio/kofams/.snakemake_timestamp", "")
        for index, fasta in [(1, ">c_000000000001\nMAAA\n>c_000000000002\nMCCC\n"),
                             (2, ">c_000000000003\nMEEE\n")]:
            chunk = f"catalogue.part_{index:03}"
            chunk_fasta = self.write(f"{self.base}/break_fasta_apart/catalogue/{chunk}.fasta", fasta)
            calls = self.directory / f"{self.base}/make_gene_calls/catalogue/{chunk}.tsv"
            calls.parent.mkdir(parents=True, exist_ok=True)
            subprocess.run(["bash", str(SCRIPTS / "make_external_gene_calls.sh"), str(chunk_fasta), str(calls)], check=True)
            self.write(f"{self.base}/anvio_contigs_db/catalogue{chunk}.db", "cached engine output")
            self.write(f"results/temp/handoff-test/annotation/anvio_kegg_annotation/catalogue/{chunk}_anvio_run_kegg_kofams.0", "")
            functions = "gene_callers_id\tsource\taccession\tfunction\te_value\n"
            functions += "0\tKOfam\tK00001\tfirst function\t1e-20\n"
            if index == 1:
                functions += "1\tKOfam\tK00002\tsecond function\t1e-15\n"
            self.write(f"{self.base}/anvio_export_functions/catalogue_{chunk}.tsv", functions)
            annotations = EMAPPER_HEADER
            for line in fasta.splitlines():
                if line.startswith(">"):
                    annotations += f"{line[1:]}\tGO:0000001,GO:0000002\tseed\t1e-20\t100\n"
            self.write(f"{self.base}/emapper/catalogue/{chunk}/go_annotations.emapper.annotations", annotations)
        self.write(f"{self.base}/break_fasta_apart/catalogue/.snakemake_timestamp", "")
        self.write("resources/go_ontology/obo/go.obo", OBO)
        for name in ["eggnog_proteins.dmnd", "eggnog.db", "eggnog.taxa.db", "eggnog.taxa.db.traverse.pkl"]:
            self.write("resources/eggnog/" + name, "cached test resource")
        wrapper = self.write("bin/python", f"#!{sys.executable}\n" +
                             "import os, sys\nfrom pathlib import Path\n" +
                             f"interpreter = {sys.executable!r}\n" +
                             "if Path(sys.argv[1]).name in ['annotate_go_terms.py', 'export_go_annotations.py']:\n" +
                             f"    interpreter = {GO_PYTHON!r}\n" +
                             "elif Path(sys.argv[1]).name == 'sequence_to_annotation_linker.py':\n" +
                             f"    interpreter = {MATRIX_PYTHON!r}\n" +
                             "os.execv(interpreter, [interpreter, *sys.argv[1:]])\n")
        wrapper.chmod(0o755)

    def test_both_handoffs_produce_matching_species_family_ids(self):
        marker = f"{self.base}/run_anntoation/catalogue.done"
        legacy = f"{self.base}/reannotate_blast_results/reference_catalogue.tsv"
        allowed = ["break_fasta_apart", "get_headers", "make_gene_calls", "combine_emapper", "make_go_db",
                   "link_information", "combine_link_info", "export_go_annotations", "link_nodes",
                   "link_go_nodes", "run_annotation", "parse_blast", "reannotate_blast_results"]
        result = subprocess.run([
            sys.executable, "-m", "snakemake", marker, legacy, "--snakefile", "workflow/Snakefile",
            "--cores", "1", "--nolock", "--allowed-rules", *allowed,
        ], cwd=self.directory, env={
            **os.environ, "PATH": str(self.directory / "bin") + os.pathsep + os.environ["PATH"],
            "PYTHONDONTWRITEBYTECODE": "1", "XDG_CACHE_HOME": str(self.directory / "cache"),
        }, capture_output=True, text=True, timeout=120)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        expected = {("S1.1", "0002"), ("S2.1", "0002"), ("S2.1", "denovo-family"), ("S3.1", "NA")}
        for rule in ["link_nodes", "link_go_nodes"]:
            output = self.directory / f"{self.base}/{rule}/catalogue/node_annotations.tsv"
            rows = self.rows(output, ",")
            self.assertEqual({(row["linker_info"], row["node_head"]) for row in rows}, expected)
            if rule == "link_go_nodes":
                self.assertEqual(len(rows), 8)
        self.assertTrue((self.directory / marker).is_file())
        self.assertEqual(self.rows(self.directory / legacy)[0]["query"], "q1")
        self.assertEqual(self.rows(self.directory / legacy)[0]["status"], "active")
        self.assertNotIn("mmseqs search", result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
