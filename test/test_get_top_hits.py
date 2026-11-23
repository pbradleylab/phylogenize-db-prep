"""CLI regression tests for exact sequence IDs and highest-identity matching."""

import csv
import importlib.util
from pathlib import Path
import random
import subprocess
import sys
import tempfile
import unittest


SCRIPT = Path(__file__).resolve().parents[1] / "workflow/scripts/get_top_hits.py"


@unittest.skipUnless(importlib.util.find_spec("pandas"), "pandas is required")
class GetTopHitsTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="top hits ")
        self.addCleanup(self.temporary.cleanup)
        self.input = Path(self.temporary.name) / "input.tsv"
        self.output = Path(self.temporary.name) / "output.tsv"

    def select(self, rows):
        with self.input.open("w") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["query", "target", "pident"])
            writer.writerows(rows)
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--input", str(self.input), "--output", str(self.output)],
            capture_output=True, text=True, timeout=30,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        with self.output.open() as handle:
            selected = list(csv.reader(handle, delimiter="\t"))
        self.assertTrue(all(len(row) == 3 for row in selected))
        self.assertEqual(len(selected), len({row[0] for row in selected}))
        return {query: (target, float(pident)) for query, target, pident in selected}

    def test_queries_with_leading_zeros_remain_distinct(self):
        self.assertEqual(self.select([
            ("001", "family_A", 90), ("1", "family_B", 99),
            ("001", "family_C", 95),
        ]), {"001": ("family_C", 95), "1": ("family_B", 99)})

    def test_target_ids_keep_leading_zeros(self):
        self.assertEqual(self.select([
            ("query_A", "0001", 99), ("query_A", "1", 90),
            ("query_B", "0010", 95),
        ]), {"query_A": ("0001", 99), "query_B": ("0010", 95)})

    def test_decimal_and_scientific_notation_ids_are_literal(self):
        identifiers = ["1", "1.0", "01", "1e3", "1000"]
        self.assertEqual(self.select([
            (identifier, "family_" + identifier, 90) for identifier in identifiers
        ]), {identifier: ("family_" + identifier, 90) for identifier in identifiers})

    def test_large_numeric_ids_keep_their_precision(self):
        identifiers = ["9007199254740992", "9007199254740993", "1.5"]
        self.assertEqual(self.select([
            (identifier, "family_" + identifier, 90) for identifier in identifiers
        ]), {identifier: ("family_" + identifier, 90) for identifier in identifiers})

    def test_missing_value_marker_names_are_valid_ids(self):
        identifiers = ["NA", "NaN", "NULL", "None", "N/A"]
        self.assertEqual(self.select([
            (identifier, identifier, 90) for identifier in identifiers
        ]), {identifier: (identifier, 90) for identifier in identifiers})

    def test_highest_identity_is_compared_numerically(self):
        self.assertEqual(self.select([
            ("query_A", "family_9", "9"),
            ("query_A", "family_99", "99.9"),
            ("query_A", "family_100", "100"),
        ]), {"query_A": ("family_100", 100)})

    def test_ties_keep_the_first_maximum_hit(self):
        self.assertEqual(self.select([
            ("query_A", "lower", 90), ("query_A", "first", 99),
            ("query_A", "second", 99),
        ]), {"query_A": ("first", 99)})

    def test_normal_catalogue_matching_matches_expected_first_maximum(self):
        rng = random.Random(7281)
        rows = [
            (f"MGYG{query:09}_001", f"UniRef50_family_{hit}", rng.choice([51, 88.5, 99, 100]))
            for query in range(50) for hit in range(20)
        ]
        rng.shuffle(rows)
        expected = {}
        for query, target, score in rows:
            if query not in expected or score > expected[query][1]:
                expected[query] = (target, score)
        self.assertEqual(self.select(rows), expected)

    def test_header_only_input_still_produces_an_empty_headerless_file(self):
        self.assertEqual(self.select([]), {})
        self.assertEqual(self.output.read_text(), "")


if __name__ == "__main__":
    unittest.main()
