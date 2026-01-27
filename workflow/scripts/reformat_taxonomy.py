""" This script is for combining species into a larger matrix using polars. The 
matrix  has protein names on the Y axis and the species accession on the top. Values
are a binary for if the centroid is present for the different predefined species.

In our case, this is used for the nucleotide clustered pangenomes for midas' GTDB 
database. 
"""
import argparse
import polars as pl
import sys

# Makes the taxonomy file.
def transform(df, mapping, tax, args):
    tmp=df.join(mapping, on="query", how="left")
    tmp=tmp.with_columns([
        pl.col("query").str.split(args.split_char).list.get(0).alias("query")])
   
    print(tmp)
    sys.exit()
    tax=tax["Genome","Species_rep"]
    tax.columns=["query","cluster"]

    tax=tax.with_columns(pl.col("query").cast(pl.Utf8))
    tmp=tmp.with_columns(pl.col("query").cast(pl.Utf8))
    tmp=tmp.with_columns(pl.col("target").cast(pl.Utf8))
    

    print(tax)
    print(tmp)
    sys.exit()
    tmp=tmp.join(tax, on="query", how="left")
   
    out=tmp.with_columns([
        pl.col("other").str.split(args.split_char).list.get(0).alias("other")])
    
    return(out)

def translate(df, tax):
    tax=tax.with_columns([
            pl.col("Lineage").str.replace_all(r"\b[a-z]__+", "").str.replace_all(r" ", " ").str.split(";"),
            pl.col("Genome").str.split(args.split_char).list.get(0).alias("query")])
    tax=tax.with_columns([
        pl.col("Lineage").list.get(i).alias(name)
        for i, name in enumerate(["domain", "phylum", "class", "order", "family", "genus", "species"])]).drop("Lineage")
    
    merged = df.join(tax["query", "domain", "phylum", "class", "order", "family", "genus", "species"], on="query", how="left").drop(["other","target","query"])
    
    return(merged)

def write_tax(tax, args):
    final_tax=tax["domain", "phylum", "class", "order", "family", "genus", "species","cluster"].unique()
    final_tax.write_csv(args.tax_output, separator=",")

def main(args):
    df=pl.read_csv(args.input, separator="\t", truncate_ragged_lines=True)
    tax=pl.read_csv(args.tax, separator="\t")
    mapping=pl.read_csv(args.mapping, separator="\t", has_header=True, new_columns=["query","other"])
    out_binary=transform(df, mapping, tax, args)
    
    del(df); del(mapping)

    print("Transformed the binary")
    out_tax=translate(out_binary, tax)
    print("Translated the taxonomy")
    write_tax(out_tax, args)
    out_binary.write_csv(args.output, separator=",")
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output","-o",
            help = "The matrix file to write to")
    parser.add_argument("--input","-i",
            help = "Directory containing all the files to combine")
    parser.add_argument("--split_char","-s",
            help = "Which charcter string to split on")
    parser.add_argument("--mapping","-m",
            help = "Mapping file with all centroids mapping to their nodes")
    parser.add_argument("--tax","-t",
            help = "input taxonomy file")
    parser.add_argument("--tax_output","-to",
            help = "Taxonomy file to be made for phylogenize database")
    args=parser.parse_args()
    main(args)
