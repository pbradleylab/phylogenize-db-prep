"""Link family annotations directly to the species/family pairs used by matrices."""

import argparse

import polars as pl


def read_table(path, separator, columns):
    frame = pl.read_csv(path, separator=separator, infer_schema=False).select(columns)
    if frame.filter(pl.any_horizontal(pl.all().is_null() | (pl.all().str.strip_chars() == ""))).height:
        raise ValueError(f"Missing identifiers or annotations in {path}")
    return frame.unique(maintain_order=True)


def link_annotations(hits_path, annotations_path, pairs_path, output, separator=","):
    hits = read_table(hits_path, "\t", ["query", "target"])
    if hits["query"].is_duplicated().any():
        raise ValueError("A protein centroid maps to more than one target family")
    pairs = read_table(pairs_path, ",", ["cluster", "target"])
    annotations = pl.read_csv(annotations_path, separator="\t", infer_schema=False)
    family_column = "node_head" if "node_head" in annotations.columns else "linker_info"
    annotations = annotations.select(
        pl.col(family_column).alias("node_head"), "accession", "function"
    ).unique(maintain_order=True)
    if annotations.filter(pl.any_horizontal(pl.all().is_null() | (pl.all().str.strip_chars() == ""))).height:
        raise ValueError("Family annotations contain missing identifiers or functions")
    valid_families = hits.select("target").unique()
    for frame, column in [(pairs, "target"), (annotations, "node_head")]:
        if frame.join(valid_families, left_on=column, right_on="target", how="anti").height:
            raise ValueError("Species pairs or annotations contain families absent from combined hits")
    merged = pairs.join(annotations, left_on="target", right_on="node_head", how="inner")
    merged = merged.select(
        pl.col("cluster").alias("linker_info"), pl.col("target").alias("node_head"),
        "accession", "function",
    ).unique().sort(["linker_info", "node_head", "accession", "function"])
    merged.write_csv(output, separator=separator, include_header=True)
    print(f"Linked {annotations.height} family annotations to {merged.height} species/family annotations")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ["hits", "annotations", "pairs", "output"]:
        parser.add_argument(name)
    parser.add_argument("--separator", choices=[",", "tab"], default=",",
                        help="Output delimiter; comma preserves the existing Phylogenize format")
    args = parser.parse_args()
    link_annotations(args.hits, args.annotations, args.pairs, args.output,
                     separator="\t" if args.separator == "tab" else ",")
