"""Select exactly one representative per final family, including de novo families."""

import argparse
from contextlib import closing
import hashlib
import shutil
import subprocess

from annotation_io import fasta_records, read_table


def selected_records(path, identifiers):
    selector = shutil.which("faSomeRecords")
    if selector is None:
        yield from fasta_records(path)
        return
    command = [selector, str(path), str(identifiers), "stdout"]
    process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)
    try:
        yield from fasta_records(process.stdout)
    finally:
        process.stdout.close()
        returncode = process.wait()
    if returncode:
        raise subprocess.CalledProcessError(returncode, command)


def select_sequences(hits, references, denovo, output, identifiers):
    wanted = {row["target"] for row in read_table(hits, ["query", "target"])}
    seen = {}
    # Retain the existing C-based FASTA selector for large reference databases.
    with open(identifiers, "w") as mapping:
        mapping.writelines(family + "\n" for family in sorted(wanted))
    with open(output, "w") as fasta:
        for path in [*references, denovo]:
            with closing(selected_records(path, identifiers)) as records:
                for family, sequence in records:
                    if family not in wanted:
                        continue
                    checksum = hashlib.sha256(sequence.encode("ascii")).digest()
                    if family in seen:
                        if seen[family] != checksum:
                            raise ValueError(f"Conflicting representative sequences for family {family!r}")
                        continue
                    seen[family] = checksum
                    fasta.write(f">{family}\n{sequence}\n")
    with open(identifiers, "w") as mapping:
        mapping.writelines(family + "\n" for family in seen)
    missing = wanted - seen.keys()
    if missing:
        raise ValueError("Missing representatives for final families: " + ", ".join(sorted(missing)[:20]))
    if not wanted:
        raise ValueError("Combined hits contain no families to annotate")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hits", required=True)
    parser.add_argument("--references", nargs="+", required=True)
    parser.add_argument("--denovo", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--identifiers", required=True)
    args = parser.parse_args()
    select_sequences(args.hits, args.references, args.denovo, args.output, args.identifiers)
