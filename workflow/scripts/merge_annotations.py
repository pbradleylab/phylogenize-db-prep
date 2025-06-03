import pandas as pd
import argparse


# Read in the arguments
parser = argparse.ArgumentParser(description="Merge anvio function and gene calls together to make an annotated fasta")
parser.add_argument('functions', help='functions file')
parser.add_argument('gene_calls', help='gene calls')
parser.add_argument('linker', help='linker')
parser.add_argument('output', help='merged file')
args = parser.parse_args()

# merge the files together by the gene calling id
functions=pd.read_csv(args.functions, sep='\t')
genecalls=pd.read_csv(args.gene_calls, sep='\t')
linker=pd.read_csv(args.linker, sep='\t')

merged=pd.merge(functions, genecalls, on="gene_callers_id")
merged=pd.merge(merged, linker, on="gene_callers_id")
filt=merged[merged["source_x"]=="KOfam"]
filt=filt[["gene_callers_id","linker_info", "accession","function"]]

filt.to_csv(args.output, index=False, sep="\t")
