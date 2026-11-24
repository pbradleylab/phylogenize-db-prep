"""Check that merging eggNOG chunks preserves family IDs and annotation columns."""

import importlib.util
from pathlib import Path
import tempfile
import unittest


REPOSITORY = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    "combine_emapper", REPOSITORY / "workflow/scripts/combine_emapper_annotations.py"
)
combine_emapper = importlib.util.module_from_spec(spec)
spec.loader.exec_module(combine_emapper)

HEADER = "#query\tseed_ortholog\tevalue\tscore\teggNOG_OGs\tmax_annot_lvl\tCOG_category\tDescription\tPreferred_name\tGOs\n"


class CombineEmapperTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="emapper chunks ")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        self.output = self.directory / "combined.annotations"
        self.fastas = []
        self.linkers = []
        self.annotations = []
        for index, family in enumerate(["0002", "family-with-hyphens", "NA"], start=1):
            query = f"c_{index:012}"
            fasta = self.directory / f"chunk{index}.fasta"
            fasta.write_text(f">{query}\nMAAAA\n")
            linker = self.directory / f"chunk{index}.tsv"
            linker.write_text(f"gene_callers_id\tlinker_info\n0\t{family}\n")
            annotations = self.directory / f"chunk{index}.annotations"
            annotations.write_text(
                "## emapper metadata\n" + HEADER +
                f"{query}\tseed\t1e-20\t100\tOG\t2\tS\tdescription\tname\tGO:0008150,GO:0003674\n" +
                "## 1 queries scanned\n"
            )
            self.fastas.append(fasta)
            self.linkers.append(linker)
            self.annotations.append(annotations)

    def combine(self):
        combine_emapper.combine_annotations(self.annotations, self.fastas, self.linkers, self.output)

    def test_restores_literal_family_ids_and_keeps_one_header(self):
        self.combine()
        lines = self.output.read_text().splitlines()
        self.assertEqual(sum(line.startswith("#query\t") for line in lines), 1)
        rows = [line.split("\t") for line in lines if not line.startswith("#")]
        self.assertEqual([row[0] for row in rows], ["0002", "family-with-hyphens", "NA"])
        self.assertTrue(all(row[9] == "GO:0008150,GO:0003674" for row in rows))

    def test_accepts_a_chunk_without_hits(self):
        self.annotations[1].write_text(HEADER)
        self.combine()
        rows = [line.split("\t")[0] for line in self.output.read_text().splitlines()
                if not line.startswith("#")]
        self.assertEqual(rows, ["0002", "NA"])

    def test_accepts_all_chunks_without_hits(self):
        for annotations in self.annotations:
            annotations.write_text(HEADER)
        self.combine()
        self.assertEqual(self.output.read_text().splitlines()[-1], HEADER.rstrip("\n"))

    def test_maps_within_chunk_order_instead_of_global_anvio_number(self):
        self.fastas[1].write_text(">c_000000010001\nMAAAA\n>c_000000010002\nMAAAA\n")
        self.linkers[1].write_text("gene_callers_id\tlinker_info\n1\tsecond\n0\tfirst\n")
        self.annotations[1].write_text(HEADER + "c_000000010002\tseed\t1e-20\t100\tOG\t2\tS\tdescription\tname\t-\n")
        self.combine()
        self.assertIn("second\tseed\t", self.output.read_text())

    def test_rejects_inconsistent_mapping(self):
        self.linkers[1].write_text("gene_callers_id\tlinker_info\n1\twrong-index\n")
        with self.assertRaisesRegex(ValueError, "FASTA and header mapping disagree"):
            self.combine()

    def test_rejects_unknown_queries(self):
        self.annotations[1].write_text(HEADER + "unknown\tseed\t1e-20\t100\tOG\t2\tS\tdescription\tname\tGO:0008150\n")
        with self.assertRaisesRegex(ValueError, "Unknown query ID"):
            self.combine()

    def test_rejects_duplicate_family_annotations(self):
        self.linkers[1].write_text("gene_callers_id\tlinker_info\n0\t0002\n")
        with self.assertRaisesRegex(ValueError, "Duplicate annotated family ID"):
            self.combine()

    def test_rejects_mismatched_annotation_columns(self):
        self.annotations[1].write_text(HEADER.replace("GOs", "other_column"))
        with self.assertRaisesRegex(ValueError, "Inconsistent eggNOG columns"):
            self.combine()


if __name__ == "__main__":
    unittest.main()
