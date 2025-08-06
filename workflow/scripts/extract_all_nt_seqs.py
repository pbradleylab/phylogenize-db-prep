#!/usr/bin/env python
# Script to make a single 16S file, named appropriately, from a MGnify Genomes database directory

import click
import os
import re
import gzip
import logging
import csv
from multiprocessing import Pool
from functools import reduce
from functools import partial
from Bio import Seq, SeqIO, SeqRecord
from gff3 import Gff3

logger = logging.getLogger(__name__)

@click.command()
@click.option('--input_path', '-i', default=".", help="Path to traverse for GFF files for given species")
@click.option("--species", '-s', help="Species representative name (e.g. MGYG000000010)")
@click.option('--log_file', '-l', default="/dev/null", help="Optional log file")
@click.option('--output_file', '-o', default="ssu_output.fa", help="Output FASTA file")
@click.option('--n_processes', '-n', default=1, help="Launch this many processes at a time")
# Main script logic. Weird comments are to disable complaints from linters that don't know about click options
def run(input_path, species, log_file, output_file, n_processes):
    logging.basicConfig(filename=log_file, level=logging.INFO) # noqa: F821  # pyright: ignore
    logging.info("Finding/parsing GFFs...\n")
    all_parsed = parse_all_gffs(seqs, input_path, n_processes=n_processes)  # noqa: F821  # pyright: ignore
    logging.info("Writing FASTA output...\n")
    SeqIO.write(iter_nested(all_parsed), output_file, "fasta")  # noqa: F821 # pyright: ignore

# SeqIO.write wants an iterable, but our sequences are buried in a nested dict. This one also automatically skips empty entries
def iter_nested(d):
    for k in d.keys():
        if (len(d[k]) == 0):
            continue
        for v in d[k].values():
            yield v


# get all matching sequences from traversing a path - multithreaded
def parse_all_gffs(seq_list, input_path=".", n_processes=4):
    genomes = dict()
    with Pool(n_processes) as p:
        subdicts = p.map(
            partial(read_and_parse_gff, seq_list=seq_list),
            os.walk(input_path, topdown=False))
    genomes = reduce(lambda a, b: a | b, subdicts, {})
    return(genomes)

# Look for and parse gffs in a given directory given a list of sequences to match
def read_and_parse_gff(x):
    (r, ds, fs) = x
    subdict = {}
    for f in fs:
        f_path = os.path.join(r, f)
        if f.endswith(".gff.gz"):
            f_name = re.sub(".gff.gz$", "", f)
            with gzip.open(f_path, 'rt') as fh:
                logging.info(f"GFF file found at {f_path}\n")
                subdict[f_name] = extract_seqs_in_list(fh)
        if f.endswith(".gff"):
            f_name = re.sub(".gff$", "", f)
            with open(f_path, 'r') as fh:
                logging.info(f"GFF file found at {f_path}\n")
                subdict[f_name] = extract_seqs_in_list(fh)
    return(subdict)


# Parse a single GFF file and return any nucl. sequences for proteins that match a given list
def extract_seqs_in_list(gff_handle):
    parsed = Gff3(logger=logger)
    logging.disable(logging.CRITICAL) # temporarily turn off logging since this is very verbose
    parsed.parse(gff_handle)
    logging.disable(logging.NOTSET)
    ks = parsed.features.keys()
    seqs = dict()
    for k in ks:
        for feature in parsed.features[k]:
            if "type" in feature.keys():
                if feature["type"] == "CDS":
                    name = feature["attributes"]["ID"]
                    logging.info(f"Found a CDS that matched called {name}")
                    seqs[name] = extract_seq_from_contig(parsed, feature)
    return(seqs)

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
