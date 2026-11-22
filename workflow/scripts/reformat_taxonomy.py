"""Map protein-family members to canonical species and their taxonomy."""

import argparse
from pathlib import Path

import polars as pl


def genome_key(column):
    return pl.col(column).str.replace(r"\.[0-9]+$", "")


def load_metadata(path):
    metadata = pl.read_csv(path, separator="\t", infer_schema=False).select(
        "Genome", "Species_rep", "Lineage"
    )
    if metadata.is_empty() or metadata.filter(
        pl.any_horizontal(pl.all().is_null())
        | pl.any_horizontal(pl.all().str.strip_chars() == "")
    ).height:
        raise ValueError("Genome metadata contains missing IDs or lineages")
    metadata = metadata.with_columns(genome_key("Genome").alias("genome_key"))
    collisions = metadata.filter(pl.col("genome_key").is_duplicated())
    if collisions.height:
        raise ValueError(
            "Genome IDs are ambiguous after version normalization: "
            + ", ".join(collisions["Genome"].head(20).to_list())
        )
    return metadata


def make_taxonomy(metadata):
    ranks = ["domain", "phylum", "class", "order", "family", "genus", "species"]
    representatives = metadata.filter(pl.col("Genome") == pl.col("Species_rep"))
    missing = set(metadata["Species_rep"]) - set(representatives["Species_rep"])
    if missing:
        raise ValueError("Species representatives missing from metadata: " + ", ".join(sorted(missing)[:20]))
    if representatives.filter(pl.col("Lineage").str.count_matches(";") != 6).height:
        raise ValueError("Representative lineages must contain seven taxonomy ranks")
    return (
        representatives.with_columns(
            pl.col("Lineage").str.split_exact(";", 6).struct.rename_fields(ranks)
        )
        .unnest("Lineage")
        .with_columns(pl.col(rank).str.replace(r"^[a-z]__", "") for rank in ranks)
        .with_columns(pl.when(pl.col("phylum").str.strip_chars() == "").then(pl.lit("unclassified")).otherwise(pl.col("phylum")).alias("phylum"))
        .select(*ranks, pl.col("Species_rep").alias("cluster"))
        .sort("cluster")
    )


def make_pairs(args, metadata):
    hits = pl.read_csv(args.input, separator="\t", infer_schema=False).select("query", "target").unique()
    if hits.is_empty() or hits.filter(
        pl.any_horizontal(pl.all().is_null())
        | pl.any_horizontal(pl.all().str.strip_chars() == "")
    ).height:
        raise ValueError("Combined hits are empty or contain missing protein/family IDs")
    if hits["query"].is_duplicated().any():
        raise ValueError("A protein centroid maps to more than one target family")
    members = pl.scan_csv(
        args.mapping,
        separator="\t",
        has_header=False,
        new_columns=["query", "other"],
        schema_overrides=[pl.String, pl.String],
    )
    pairs = (
        members.join(hits.lazy(), on="query", how="left")
        .with_columns(pl.col("other").str.split(args.split_char).list.get(0).alias("member_genome"))
        .with_columns(genome_key("member_genome").alias("genome_key"))
        .join(metadata.lazy().select("genome_key", pl.col("Species_rep").alias("cluster")), on="genome_key", how="left")
        .select(
            "cluster", "target",
            pl.when(pl.col("cluster").is_null()).then(pl.col("member_genome")).alias("unmatched_genome"),
            pl.when(pl.col("target").is_null()).then(pl.col("query")).alias("unmapped_centroid"),
        )
        .unique()
        .collect(engine="streaming")
    )
    unmatched = pairs["unmatched_genome"].drop_nulls().unique().to_list()
    unmapped = pairs["unmapped_centroid"].drop_nulls().unique().to_list()
    if unmatched:
        raise ValueError("Member genomes missing from metadata: " + ", ".join(sorted(unmatched)[:20]))
    if unmapped:
        raise ValueError("Protein centroids missing from combined hits: " + ", ".join(sorted(unmapped)[:20]))
    pairs = pairs.select("cluster", "target")
    if pairs.is_empty() or pairs.filter(pl.any_horizontal(pl.all().is_null())).height:
        raise ValueError("Member mappings produced empty or incomplete species/family pairs")
    missing = set(metadata["Species_rep"]) - set(pairs["cluster"])
    if missing:
        raise ValueError("Source species missing from member mappings: " + ", ".join(sorted(missing)[:20]))
    return pairs


def main(args):
    if not args.split_char:
        raise ValueError("Protein-to-genome delimiter must not be empty")
    metadata = load_metadata(args.tax)
    taxonomy = make_taxonomy(metadata)
    pairs = make_pairs(args, metadata)
    for output in [args.output, args.tax_output]:
        Path(output).parent.mkdir(parents=True, exist_ok=True)
    pairs.write_csv(args.output)
    taxonomy.write_csv(args.tax_output)
    print(f"Wrote {pairs.height} species/family pairs for all {taxonomy.height} source species")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", "-o", required=True, help="Species/family CSV to write")
    parser.add_argument("--input", "-i", required=True, help="Combined centroid-to-family hits TSV")
    parser.add_argument("--split_char", "-s", default="_", help="Protein-to-genome delimiter")
    parser.add_argument("--mapping", "-m", required=True, help="Headerless centroid/member protein TSV")
    parser.add_argument("--tax", "-t", required=True, help="Source genome metadata TSV")
    parser.add_argument("--tax_output", "-to", required=True, help="Canonical species taxonomy CSV to write")
    main(parser.parse_args())
