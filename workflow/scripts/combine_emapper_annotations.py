"""Combine eggNOG chunks and restore the family IDs used by the KO branch."""

import csv


def read_query_mapping(fasta_path, linker_path):
    with open(fasta_path) as handle:
        queries = [line[1:].split()[0] for line in handle if line.startswith(">")]
    with open(linker_path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["gene_callers_id", "linker_info"]:
            raise ValueError(f"Unexpected header in {linker_path}")
        families = {}
        for row in reader:
            index = row["gene_callers_id"]
            family = row["linker_info"]
            if index in families or not family:
                raise ValueError(f"Duplicate index or missing family ID in {linker_path}")
            families[index] = family
    if len(set(queries)) != len(queries):
        raise ValueError(f"Duplicate FASTA IDs in {fasta_path}")
    if set(families) != {str(index) for index in range(len(queries))}:
        raise ValueError(f"FASTA and header mapping disagree: {fasta_path}, {linker_path}")
    return {query: families[str(index)] for index, query in enumerate(queries)}


def combine_annotations(annotation_paths, fasta_paths, linker_paths, output_path):
    if not annotation_paths or not (len(annotation_paths) == len(fasta_paths) == len(linker_paths)):
        raise ValueError("Expected one annotation file and header mapping per FASTA chunk")
    expected_header = None
    seen_families = set()
    with open(output_path, "w") as output:
        output.write("## Combined eggNOG-mapper chunks; original family IDs restored\n")
        for annotations, fasta, linker in zip(annotation_paths, fasta_paths, linker_paths):
            query_mapping = read_query_mapping(fasta, linker)
            header = None
            with open(annotations) as handle:
                for line in handle:
                    if line.startswith("#query\t"):
                        header = line.rstrip("\r\n").split("\t")
                        if expected_header is None:
                            expected_header = header
                            output.write("\t".join(header) + "\n")
                        elif header != expected_header:
                            raise ValueError(f"Inconsistent eggNOG columns in {annotations}")
                    elif line.startswith("#") or not line.strip():
                        continue
                    else:
                        fields = line.rstrip("\r\n").split("\t")
                        if header is None or len(fields) != len(header):
                            raise ValueError(f"Missing header or malformed annotation in {annotations}")
                        query = fields[0]
                        if query not in query_mapping:
                            raise ValueError(f"Unknown query ID {query!r} in {annotations}")
                        family = query_mapping[query]
                        if family in seen_families:
                            raise ValueError(f"Duplicate annotated family ID {family!r}")
                        seen_families.add(family)
                        fields[0] = family
                        output.write("\t".join(fields) + "\n")
            if header is None:
                raise ValueError(f"Missing eggNOG column header in {annotations}")


if __name__ == "__main__":
    combine_annotations(
        snakemake.input.annotations,
        snakemake.input.fastas,
        snakemake.input.linkers,
        snakemake.output[0],
    )
