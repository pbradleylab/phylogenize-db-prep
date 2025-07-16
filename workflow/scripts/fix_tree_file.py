#!/usr/bin/env python
# Script to make a single 16S file, named appropriately, from a MGnify Genomes database directory

import click
import re

@click.command()
@click.option('--input_path', '-i', default=".", help="Path to traverse for GFF files")
@click.option('--output_file', '-o', default="ssu_output.fa", help="Output FASTA file")
# Take only first line, remove quote characters, and replace ;; with ____
def run(input_path, output_file):
    with open(input_path, 'r') as fh:
        first_line = fh.readline()
        dequote = re.sub("'", "", first_line)
        redelim = re.sub(";;", "____", dequote)
    with open(output_file, 'w') as fh:
        fh.write(redelim)

# Run the script if appropriate
if __name__ == '__main__':
    run()
