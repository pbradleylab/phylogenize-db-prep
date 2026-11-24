"""Restore family IDs for a chunk using Anvi'o's explicit FASTA rename report."""

import argparse
import csv

from annotation_io import fasta_records


def make_headers(fasta, report, output, identifiers):
    queries = [identifier for identifier, _ in fasta_records(fasta)]
    if len(set(queries)) != len(queries):
        raise ValueError(f"Duplicate sequence IDs in {fasta}")
    wanted = set(queries)
    families = {}
    with open(report) as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t", 1)
            if len(fields) != 2 or not fields[1].split():
                raise ValueError(f"Malformed FASTA rename report: {report}")
            if fields[0] in wanted:
                if fields[0] in families:
                    raise ValueError(f"Repeated sequence ID in rename report: {fields[0]}")
                families[fields[0]] = fields[1].split()[0]
                if len(families) == len(wanted):
                    break
    if wanted != families.keys():
        raise ValueError(f"FASTA chunk IDs missing from rename report: {fasta}")
    if len(set(families.values())) != len(families):
        raise ValueError(f"Duplicate family IDs in {fasta}")
    with open(output, "w", newline="") as handle, open(identifiers, "w") as ids:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_callers_id", "linker_info"])
        for index, query in enumerate(queries):
            writer.writerow([index, families[query]])
            ids.write(families[query] + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ["fasta", "report", "output", "identifiers"]:
        parser.add_argument("--" + option, required=True)
    args = parser.parse_args()
    make_headers(args.fasta, args.report, args.output, args.identifiers)
