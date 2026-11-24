"""Small streaming readers shared by the annotation handoff scripts."""

import csv
from contextlib import nullcontext


def read_table(path, required, delimiter="\t"):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not set(required).issubset(reader.fieldnames or []):
            raise ValueError(f"Missing required columns {required} in {path}")
        for number, row in enumerate(reader, start=2):
            if None in row or any(row[key] is None or not row[key].strip() for key in required):
                raise ValueError(f"Malformed or incomplete row {number} in {path}")
            yield row


def fasta_records(path):
    identifier = None
    sequence = []
    with (nullcontext(path) if hasattr(path, "read") else open(path)) as handle:
        for number, line in enumerate(handle, start=1):
            if line.startswith(">"):
                if identifier is not None:
                    if not sequence:
                        raise ValueError(f"Empty protein {identifier!r} in {path}")
                    yield identifier, "".join(sequence)
                fields = line[1:].split()
                if not fields:
                    raise ValueError(f"Empty FASTA identifier at line {number} in {path}")
                identifier, sequence = fields[0], []
            elif line.strip():
                if identifier is None:
                    raise ValueError(f"Sequence before FASTA header in {path}")
                sequence.append("".join(line.split()).upper())
    if identifier is not None:
        if not sequence:
            raise ValueError(f"Empty protein {identifier!r} in {path}")
        yield identifier, "".join(sequence)


def emapper_rows(path):
    with open(path, newline="") as handle:
        header = None
        for line in handle:
            if line.startswith("#query\t"):
                if header is not None:
                    raise ValueError(f"Repeated eggNOG column header in {path}")
                header = line.rstrip("\r\n").split("\t")
                if "GOs" not in header:
                    raise ValueError(f"Missing GOs column in {path}")
            elif not line.strip() or line.startswith("#"):
                continue
            else:
                fields = line.rstrip("\r\n").split("\t")
                if header is None or len(fields) != len(header) or not fields[0]:
                    raise ValueError(f"Malformed eggNOG annotation in {path}")
                yield dict(zip(header, fields))
        if header is None:
            raise ValueError(f"Missing eggNOG column header in {path}")
