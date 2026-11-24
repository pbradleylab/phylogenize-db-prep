"""Combine chunk annotations with one header and validate family ownership."""

import argparse
import csv

from annotation_io import read_table
from merge_annotations import KO_COLUMNS


def combine_annotations(paths, output):
    families = set()
    with open(output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(KO_COLUMNS)
        for path in paths:
            chunk_families = set()
            for row in read_table(path, KO_COLUMNS):
                family = row["linker_info"]
                if family in families:
                    raise ValueError(f"Family {family!r} occurs in more than one KO chunk")
                chunk_families.add(family)
                writer.writerow([row[column] for column in KO_COLUMNS])
            families.update(chunk_families)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True)
    parser.add_argument("chunks", nargs="*")
    args = parser.parse_args()
    combine_annotations(args.chunks, args.output)
