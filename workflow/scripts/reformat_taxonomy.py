""" This script is for combining species into a larger matrix using polars. The 
matrix  has protein names on the Y axis and the species accession on the top. Values
are a binary for if the centroid is present for the different predefined species.

In our case, this is used for the nucleotide clustered pangenomes for midas' GTDB 
database. 
"""
import argparse
import polars as pl

def transform(df, mapping):
    df=df.join(mapping, on="query", how="left")
    assert(True == False)
    df=df.with_columns([
        pl.col("query").str.split(args.split_char).list.get(0).alias("Genome")])
    prefix_ids=(df.select("Genome").unique().with_columns([(pl.arange(100001, 100001 + pl.len()).alias("cluster"))]))
    df=df.join(prefix_ids, on="Genome", how="left").drop("query")
    return(df)

def translate(df, tax):
    merged = df.join(tax["Genome", "Species_rep", "Lineage"], on="Genome", how="left")
    merged = merged.with_columns(
        pl.col("Lineage").str.replace_all(r"\b[a-z]__+", "").str.replace_all(r"_", " ").str.split(";"))
    merged = merged.with_columns([
        pl.col("Lineage").list.get(i).alias(name)
        for i, name in enumerate(["domain", "phylum", "class", "order", "family", "genus", "species"])]).drop("Lineage")
    return(merged)

def write_tax(tax, args):
    final_tax=tax["domain", "phylum", "class", "order", "family", "genus", "species","cluster"]
    final_tax.write_csv(args.tax_output, separator=",")

def main(args):
    df=pl.read_csv(args.input, separator="\t", truncate_ragged_lines=True)
    tax=pl.read_csv(args.tax, separator="\t"),
    mapping=pl.read_csv(args.tax, separator="\t", has_header=False, new_columns=["query","target"])
    df=translate(transform(df, mapping), tax)

    write_tax(df, args)
    df.write_csv(args.output, separator=",")
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output","-o",
            help = "The matrix file to write to")
    parser.add_argument("--input","-i",
            help = "Directory containing all the files to combine")
    parser.add_argument("--split_char","-s",
            help = "Which charcter string to split on")
    parser.add_argument("--mapping","-t",
            help = "Mapping file with all centroids mapping to their nodes")
    parser.add_argument("--tax_output","-to",
            help = "Taxonomy file to be made for phylogenize database")
    args=parser.parse_args()
    main(args)
