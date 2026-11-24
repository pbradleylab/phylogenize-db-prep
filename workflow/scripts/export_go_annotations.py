"""Export direct eggNOG family/GO annotations and their supporting evidence."""

import argparse
import csv

from annotation_io import emapper_rows, read_table
from annotate_go_terms import load_ontology, term_metadata


def export_annotations(annotations, hits, obo, families_output, evidence_output, go_evidence):
    ontology = load_ontology(obo)
    wanted = {row["target"] for row in read_table(hits, ["query", "target"])}
    columns = ["node_head", "accession", "function", "aspect"]
    evidence_columns = [*columns, "original_go_id", "status", "source", "go_evidence",
                        "seed_ortholog", "seed_evalue", "seed_score", "ontology_version"]
    seen = set()
    with open(families_output, "w", newline="") as families, open(evidence_output, "w", newline="") as evidence:
        family_writer = csv.DictWriter(families, fieldnames=columns, delimiter="\t", lineterminator="\n")
        evidence_writer = csv.DictWriter(evidence, fieldnames=evidence_columns, delimiter="\t", lineterminator="\n")
        family_writer.writeheader()
        evidence_writer.writeheader()
        for row in emapper_rows(annotations):
            family = row["#query"]
            if family not in wanted or family in seen:
                raise ValueError(f"Unknown or repeated annotated family {family!r}")
            seen.add(family)
            if row["GOs"] in ["", "-"]:
                continue
            written = set()
            for go_id in sorted(set(row["GOs"].split(","))):
                metadata = term_metadata(go_id.strip(), ontology)
                record = {"node_head": family, "accession": metadata["go_id"],
                          "function": metadata["name"], "aspect": metadata["aspect"]}
                evidence_writer.writerow({
                    **record, "original_go_id": go_id, "status": metadata["status"],
                    "source": "eggNOG-mapper", "go_evidence": go_evidence,
                    "seed_ortholog": row.get("seed_ortholog", ""),
                    "seed_evalue": row.get("evalue", ""), "seed_score": row.get("score", ""),
                    "ontology_version": ontology.data_version or ontology.version,
                })
                if metadata["status"] == "active" and metadata["go_id"] not in written:
                    family_writer.writerow(record)
                    written.add(metadata["go_id"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for option in ["annotations", "hits", "obo", "families", "evidence", "go-evidence"]:
        parser.add_argument("--" + option, required=True)
    args = parser.parse_args()
    export_annotations(args.annotations, args.hits, args.obo, args.families, args.evidence, args.go_evidence)
