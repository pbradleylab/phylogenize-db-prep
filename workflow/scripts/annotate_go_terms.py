"""Describe GO IDs or protein/GO/score triples without losing their linkage."""

import argparse
import csv
import math
import re

from goatools.obo_parser import GODag


def load_ontology(path):
    return GODag(str(path), optional_attrs={"def"}, load_obsolete=True, prt=None)


def term_metadata(go_id, ontology):
    if not re.fullmatch(r"GO:[0-9]{7}", go_id):
        raise ValueError(f"Malformed GO identifier: {go_id!r}")
    term = ontology.get(go_id)
    if term is None:
        return {"go_id": go_id, "name": "NOT FOUND", "aspect": "?", "hint": "", "status": "unknown"}
    aspect = {"biological_process": "P", "molecular_function": "F", "cellular_component": "C"}.get(term.namespace, "?")
    definition = getattr(term, "defn", "").strip()
    return {
        "go_id": term.id, "name": term.name, "aspect": aspect,
        "hint": definition.split(".")[0].strip('"'),
        "status": "obsolete" if term.is_obsolete else "active",
    }


def annotate_terms(obo, go_ids, output):
    ontology = load_ontology(obo)
    columns = ["query", "go_id", "score", "original_go_id", "name", "aspect", "hint", "status"]
    with open(go_ids) as source, open(output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for line in source:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) == 1:
                query, go_id, score = "", fields[0], ""
            elif len(fields) == 3:
                query, go_id, score = fields
                if not query or not math.isfinite(float(score)):
                    raise ValueError("Invalid protein identifier or GO support score")
            else:
                raise ValueError("Expected one GO ID or a protein/GO/score triple per line")
            writer.writerow({"query": query, "score": score, "original_go_id": go_id,
                             **term_metadata(go_id, ontology)})


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--obo", required=True)
    parser.add_argument("--go_ids", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    annotate_terms(args.obo, args.go_ids, args.out)
