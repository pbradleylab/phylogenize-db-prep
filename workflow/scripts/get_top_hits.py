import pandas as pd
import argparse

def main(args):
        frame = pd.read_csv(args.input, sep='\t')
        print(frame)
        frame["pident"] = pd.to_numeric(frame["pident"])
        frame = frame.loc[frame.groupby(["query"])["pident"].idxmax()]
        frame.to_csv(args.output, sep='\t', index=None, header=None)

if __name__ == "__main__":
	parser = argparse.ArgumentParser()
	parser.add_argument("--input","-i",
                help = "The input file with a numerical statistic to take the top hit from")
	parser.add_argument("--output","-o",
                help = "The file to write to")
	args = parser.parse_args()
	main(args)
