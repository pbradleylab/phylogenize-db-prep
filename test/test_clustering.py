"""Check the clustering file contracts with MMseqs2 and faSomeRecords on PATH."""

from pathlib import Path
import random
import shutil
import subprocess
import tempfile
import unittest


MMSEQS = shutil.which("mmseqs")
FA_SOME_RECORDS = shutil.which("faSomeRecords")


@unittest.skipUnless(MMSEQS, "MMseqs2 from the clustering environment is required")
class ClusteringTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="clustering-")
        self.addCleanup(self.temporary.cleanup)
        self.directory = Path(self.temporary.name)
        generator = random.Random(37)
        sequence = "".join(generator.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(150))
        other = "".join(generator.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(150))
        self.fasta = self.directory / "proteins.faa"
        self.fasta.write_text(f">0001\n{sequence}\n>NA\n{sequence}\n>002\n{other}\n")

    def run_command(self, arguments):
        result = subprocess.run(arguments, capture_output=True, text=True, timeout=120)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)

    def test_database_creation_keeps_the_expected_index(self):
        self.run_command([MMSEQS, "createdb", str(self.fasta), str(self.directory / "queries")])
        self.assertTrue((self.directory / "queries.index").is_file())
        self.assertTrue((self.directory / "queries.dbtype").is_file())

    def test_linclust_members_and_representative_fasta_agree(self):
        self.run_command([
            MMSEQS, "easy-linclust", str(self.fasta), str(self.directory / "clusters"),
            str(self.directory / "tmp"), "--threads", "1",
        ])
        rows = [line.split("\t") for line in (self.directory / "clusters_cluster.tsv").read_text().splitlines()]
        self.assertTrue(all(len(row) == 2 for row in rows), rows)
        self.assertEqual(len(rows), 3)
        representatives = {member: representative for representative, member in rows}
        self.assertEqual(set(representatives), {"0001", "NA", "002"})
        self.assertEqual(representatives["0001"], representatives["NA"])
        fasta_ids = {
            line[1:].split()[0]
            for line in (self.directory / "clusters_rep_seq.fasta").read_text().splitlines()
            if line.startswith(">")
        }
        self.assertEqual(set(representatives.values()), fasta_ids)

    @unittest.skipUnless(FA_SOME_RECORDS, "faSomeRecords from the clustering environment is required")
    def test_fasta_selection_preserves_literal_ids(self):
        identifiers = self.directory / "ids.txt"
        identifiers.write_text("0001\nNA\n")
        output = self.directory / "selected.faa"
        self.run_command([FA_SOME_RECORDS, str(self.fasta), str(identifiers), str(output)])
        selected = [line[1:].split()[0] for line in output.read_text().splitlines() if line.startswith(">")]
        self.assertEqual(len(selected), 2)
        self.assertEqual(set(selected), {"0001", "NA"})


if __name__ == "__main__":
    unittest.main()
