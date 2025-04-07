import polars as pl
import argparse


# Read in the arguments
parser = argparse.ArgumentParser(description="Linked the ")
parser.add_argument('node', help='node list from mmseq2')
parser.add_argument('linker', help='finished matrix from main pipeline with combined cluster information')
parser.add_argument('output', help='final linker file')
args = parser.parse_args()

# merge the files together by the gene calling id
node=pl.read_csv(args.node, separator='\t', has_header=False)
node.columns=["linker_info","node_head"]
linker=pl.read_csv(args.linker, separator='\t')

# Run until one of the DataFrames is out of lines
batchn=900000
cnt=0
while not (node.is_empty() or linker.is_empty()):
    # Perform the inner join and the anti join
    merged=node.join(linker.head(batchn), on=["linker_info","node_head"], how="inner")
    node=node.join(linker.head(batchn), on=["linker_info","node_head"], how="anti")
    linker=linker.tail(linker.shape[0] - batchn)

    # Check if linker still has rows, if not break the loop
    if linker.is_empty():
        break
    # Since all genes in a single gene cluster are not annotated the same in reality, cases with 
    # multiple functions associated with a gene cluster are decided by choosing the most frequent 
    # annotation or the one that the greatest number of genes in the gene cluster matches. Ties 
    # are broken arbitrarily following what the published tool anvio does 
    # https://merenlab.org/2016/11/08/pangenomics-v2/
    grouped=merged.groupby(["node_head","accession"]).agg([pl.count().alias("count")])
    grouped=grouped.sort(by=["node_head","count"], descending=[False, True])
    filtered=grouped.groupby("node_head").agg([
        pl.col("accession").first().alias("accession"),
        pl.col("count").max().alias("max_count")
    ])
    result=filtered.join(merged, on=["node_head", "accession"])
    # Write out the dataframe to a checkpoint file
    result.write_csv(args.output+str(cnt), separator='\t', include_header=True)
    cnt=cnt+1
