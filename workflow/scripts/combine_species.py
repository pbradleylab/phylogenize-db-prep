""" This script is for combining species into a larger matrix using pandas. The 
matrix  has centroids on the Y axis and the species accession on the top. Values
are a binary for if the centroid is present for the different predefined species.

In our case, this is used for the nucleotide clustered pangenomes for midas' GTDB 
database. 
"""
import os, re
import pandas as pd
import argparse
import numpy as np


def main(args):
        centroids_dict={}
        species = [found for found in os.listdir(args.dir) if os.path.isfile(os.path.join(args.dir, found)) and found.endswith(args.ext)]
        # Get pair and unpair files in separate lists
        for specie in species:
            frame = pd.read_csv(os.path.join(args.dir, specie), delimiter = '\t')
            accessions = list(set(["".join(re.split("(\.\d*_)", x)[0:2])[:-1] for x in list(frame.iloc[:, 0])]))
            centroids = [re.split("\.\d*_", x)[1].split('_')[0] for x in list(frame.iloc[:, 0])]
            for key in centroids:
                if not key in centroids_dict.keys():
                    centroids_dict[key] = []
                else:
                    species_lst = centroids_dict[key]
                    species_lst = species_lst + accessions
                    centroids_dict[key]=list(set(species_lst))
        
        accessions_frame = pd.DataFrame({"accessions": list(centroids_dict.values())})
        accessions_matrix = accessions_frame['accessions'].apply(pd.value_counts).fillna(0).astype(int)
        ids = list(centroids_dict.keys())
        accessions_matrix.index = ids
        out = accessions_matrix.reindex(sorted(accessions_matrix.columns), axis=1)	
        print(out)

if __name__ == "__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("--output","-i",
                help = "The matrix file to write to")
        parser.add_argument("--dir","-d",
                help = "Directory containing all the files to combine")
        parser.add_argument("--ext","-e",
                help = "Extension of the files to combine")
        args = parser.parse_args()

        main(args)
