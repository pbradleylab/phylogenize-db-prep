#!/usr/bin/env python
# Script to make a single 16S file, named appropriately, from a MGnify Genomes database directory

import click
import os
import gzip
import logging
import polars as pl
from Bio import Seq, SeqIO, SeqRecord
from gff3 import Gff3

logger = logging.getLogger(__name__)

@click.command()
@click.option('--input_path', '-i', default=".", help="Path to traverse for GFF files")
@click.option('--log_file', '-l', default='extract_ssu.log', help="Log file for output")
@click.option('--metadata_file', '-m', default=None, help="Tab-separated file of genome metadata")
@click.option('--output_file', '-o', default="ssu_output.fa", help="Output FASTA file")
# Main script logic. Weird comments are to disable complaints from linters that don't know about click options
def run(input_path, log_file, metadata_file, output_file):
    logging.basicConfig(filename=log_file, level=logging.INFO) # noqa: F821  # pyright: ignore
    all_parsed = parse_all_gffs(input_path) 
    if metadata_file: # noqa: F821 # pyright: ignore
        # if provided, rename outputs 
        all_parsed = rename_16s(all_parsed, metadata_file) # noqa: F821 # pyright: ignore
    SeqIO.write(iter_nested(all_parsed), output_file, "fasta")  # noqa: F821 # pyright: ignore

# Rename 16S sequences given a genome metadata file
def rename_16s(all_parsed, metadata_file):
    nz_parsed = [n for n in all_parsed.keys() if len(all_parsed[n]) > 0]
    metadata = pl.read_csv(metadata_file, separator="\t")
    md_name_dict = dict(
            metadata.filter(
            pl.col("Genome").is_in(nz_parsed)
        ).with_columns(
            Lineage=pl.col("Lineage").str.split_exact(";", 6).
                struct.rename_fields([
                    "domain",
                    "phylum",
                    "class",
                    "order",
                    "family",
                    "genus",
                    "species"
                ])
        ).unnest("Lineage").with_columns(
            output_name = pl.concat_str(
                ["genus", "Species_rep"],
                separator=";;"
            )
        ).select(["Genome", "output_name"]).iter_rows()
    )
    for genome in nz_parsed:
        for gene in all_parsed[genome]:
            seqr = all_parsed[genome][gene]
            # Be robust to genome being missing for some reason
            if genome in md_name_dict.keys():
                seqr.id = ";;".join([gene, md_name_dict[genome]])
                seqr.name = seqr.id
                all_parsed[genome][gene] = seqr
            else:
                logging.warning(f"Genome {genome} not found in metadata")
                seqr.id = ";;".join([gene, "unknown", genome])
    return(all_parsed)

# SeqIO.write wants an iterable, but our sequences are buried in a nested dict. This one also automatically skips empty entries
def iter_nested(d):
    for k in d.keys():
        if (len(d[k]) == 0):
            continue
        for v in d[k].values():
            yield v

# Traverse a path looking for GFFs and get all 16S/SSU sequences
def parse_all_gffs(input_path="."):
    genomes = dict()
    for r, ds, fs in os.walk(input_path, topdown=False):
        for f in fs:
            f_path = os.path.join(r, f)
            if f.endswith(".gff.gz"):
                f_name = re.sub(".gff.gz$", "", f)
                with gzip.open(f_path, 'rt') as fh:
                    logging.info(f"GFF file found at {f_path}\n")
                    genomes[f_name] = get_16s(fh)
            if f.endswith(".gff"):
                f_name = re.sub(".gff$", "", f)
                with open(f_path, 'r') as fh:
                    logging.info(f"GFF file found at {f_path}\n")
                    genomes[f_name] = get_16s(fh)
    return(genomes)

# Parse a single GFF file and return any 16S sequences we find that were generated with barrnap or INFERNAL
def get_16s(gff_handle):
    parsed = Gff3(logger=logger)
    logging.disable(logging.CRITICAL) # temporarily turn off logging since this is very verbose
    parsed.parse(gff_handle)
    logging.disable(logging.NOTSET)
    ks = parsed.features.keys()
    seqs_16S = dict()
    for k in ks:
        for feature in parsed.features[k]:
            is_16S = False
            if "source" in feature.keys():
                if feature["source"].startswith("barrnap"):
                    if feature["attributes"]["product"].startswith("16S ribosomal RNA"):
                        is_16S = True
                if feature["source"].startswith("INFERNAL"):
                    # check if rfam HMM is for bacterial/archaeal SSU RNA
                    if feature["attributes"]["rfam"] in ["RF00177", "RF01959"]:
                        is_16S = True
            if is_16S:
                name = feature["attributes"]["ID"]
                logging.info(f"Found a 16S sequence called {name}")
                seqs_16S[name] = extract_seq_from_contig(parsed, feature)
                is_16S = False
    return(seqs_16S)

# Helper function to get a specific feature from the contigs at the end of the FASTA file
def extract_seq_from_contig(gff, feature):
    contig_name = feature["seqid"]
    contig = Seq.Seq(gff.fasta_embedded[contig_name]["seq"])
    gseq = contig[(feature["start"]-1):feature["end"]]
    if feature["strand"] == "-":
        gseq = gseq.reverse_complement()
    rec = SeqRecord.SeqRecord(gseq,
        id=feature["attributes"]["ID"])
    return(rec)
    
# Run the script if appropriate
if __name__ == '__main__':
    run()
