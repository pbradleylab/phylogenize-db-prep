""" Script to parse the tab delimited output of hmmsearch from hmmer3. This is run via python3.
Do. Not. Use. Python. 2.
"""
import platform
import sys
import glob
import os
from Bio import SearchIO
import argparse
from collections import defaultdict
import polars
from progress.bar import Bar

# Exit if being used with python2
def version_check():
	if platform.python_version().startswith('2'):
		print("You have tried running with python2. Only use python3")
		sys.exit(0) 

# Read directory of files with specific extension
def read_file_dir(in_dir, ext):
	return glob.glob(os.path.join(os.getcwd(), in_dir, "*"+str(ext)))

# Collect user input
def inparser():
	parser = argparse.ArgumentParser()
	parser.add_argument("--input","-i",
			required=True,
			help = "Input directory with the hmmsearch output files")
	parser.add_argument("--extension","-e",
			required=True,
			help = "Extension of the files to retrieve")
	parser.add_argument("--output","-o",
			required=True,
			help = "Output file for the specififed hits after collecting the best hit")
	args = parser.parse_args()
	return args

# Parse the hmmsearch output with Biopython
def hmm_parser(files):
	out = []
	with Bar("Processing", max = len(files)) as bar:
		for f in files:
			all_hits=defaultdict(list)
			with open(f) as hmm:
				for result in SearchIO.parse(hmm, "hmmer3-tab"):
					hits = result.hits
				
					# Gather ine info needed for all the hits	
					if len(hits) > 0:
						for i in range(0,len(hits)): 
							query = hits[i].query_id # hit decription 
							best = hits[i].evalue # evalue of hit	
							target = result[i].id # target name of hit
						
							all_hits["eval"].append(best)
							all_hits["query"].append(query)
							all_hits["target"].append(target)

			df = polars.from_dict(all_hits)
			if len(out) == 0:
				out=df
			else:
				if not len(df) == 0:
					out = polars.concat([out, df])		
			bar.next()
	return(out)


def main():
	args = inparser()
	dirs = args.input.split(',')
	
	files_lst=[]
	for folder in dirs:
		files_lst=files_lst+read_file_dir(folder, args.extension)
	df = hmm_parser(files_lst)
	df = df.group_by(["target","query"]).agg(polars.col("eval").min().alias("best_score"))
	
	highest = df.group_by(["target","query"]).first()
	highest.write_csv(args.output)
		
	

if __name__ == "__main__":
	main()
