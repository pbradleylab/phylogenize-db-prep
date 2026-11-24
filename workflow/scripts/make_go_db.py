"""Extract the named GOs column while preserving literal family identifiers."""

import argparse
import csv

from annotation_io import emapper_rows


def make_go_db(annotations, output):
    seen = set()
    with open(output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        for row in emapper_rows(annotations):
            family = row["#query"]
            if family in seen:
                raise ValueError(f"Repeated annotated family {family!r}")
            seen.add(family)
            if row["GOs"] not in ["", "-"]:
                writer.writerow([family, row["GOs"]])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("annotations")
    parser.add_argument("output")
    args = parser.parse_args()
    make_go_db(args.annotations, args.output)
