"""Join exported KOfam annotations to a chunk's literal family identifiers."""

import argparse
import csv

from annotation_io import read_table


KO_COLUMNS = ["gene_callers_id", "linker_info", "accession", "function"]


def merge_annotations(functions, gene_calls, linker, output):
    calls = [row["gene_callers_id"] for row in read_table(gene_calls, ["gene_callers_id"])]
    families = {}
    for row in read_table(linker, ["gene_callers_id", "linker_info"]):
        if row["gene_callers_id"] in families:
            raise ValueError("Duplicate gene ID in family mapping")
        families[row["gene_callers_id"]] = row["linker_info"]
    if len(set(calls)) != len(calls) or set(calls) != families.keys():
        raise ValueError("Gene calls and family mapping have different gene IDs")
    with open(output, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(KO_COLUMNS)
        seen = set()
        for row in read_table(functions, ["gene_callers_id", "source", "accession", "function"]):
            if row["gene_callers_id"] not in families:
                raise ValueError(f"Annotation refers to unknown gene ID {row['gene_callers_id']!r}")
            if row["source"] != "KOfam":
                continue
            record = (row["gene_callers_id"], families[row["gene_callers_id"]], row["accession"], row["function"])
            if record not in seen:
                writer.writerow(record)
                seen.add(record)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["functions", "gene_calls", "linker", "output"]:
        parser.add_argument(name)
    args = parser.parse_args()
    merge_annotations(args.functions, args.gene_calls, args.linker, args.output)
