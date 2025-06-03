import polars as pl
import argparse


# Read in the arguments
parser = argparse.ArgumentParser(description="Merge anvio function file from merge_annotations.py into a linker file")
parser.add_argument('functions', help='functions file')
parser.add_argument('linker', help='finished matrix from main pipeline with combined cluster information')
parser.add_argument('output', help='final linker file')
args = parser.parse_args()

# merge the files together by the gene calling id
functions=pl.read_csv(args.functions, separator='\t')
linker=pl.read_csv(args.linker, separator='\t', has_header=True)
functions.columns=["query","linker_info"]

merged=functions.join(linker, on="linker_info")
merged.columns=["linker_info","node_head","gene_id","accession","function"]

merged["linker_info","node_head","accession","function"].write_csv(args.output, include_header=True)
