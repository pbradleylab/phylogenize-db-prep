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
        dups=open(args.duplicates, "w")
        centroids_dict={}
        species = [found for found in os.listdir(args.dir) if os.path.isfile(os.path.join(args.dir, found)) and found.endswith(args.ext)]
        # Get pair and unpair files in separate lists
        frame = pd.read_csv(os.path.join(args.dir, species[0]), delimiter = '\t')
        accessions = [x.split("_")[0] for x in list(frame.iloc[:, 0])]
        centroids = list(frame.iloc[:, 1])
        
        # Make the dictionary for the genomes each protein was aligned or clustered to.
        for i in range(0,len(centroids)):
            key=centroids[i]
            if not key in centroids_dict.keys():
                centroids_dict[key] = [accessions[i]]
            else:
                species_lst = centroids_dict[key]
                species_lst.append(accessions[i])
                centroids_dict[key]=list(set(species_lst))
                # Write out duplicate hits found in target databases.
                if len(species_lst) >1:
                    dups.write(key+"\n")
                
        # Convert the dictionary into a binary matrix
        accessions_frame = pd.DataFrame({"accessions": list(centroids_dict.values())})
        accessions_matrix = accessions_frame['accessions'].apply(pd.value_counts).fillna(0).astype(int)
        ids = list(centroids_dict.keys())
        accessions_matrix.index = ids
        out = accessions_matrix.reindex(sorted(accessions_matrix.columns), axis=1)	
        
        # Test that the input number matches the output number of entries in the matrix.  
        sums = pd.DataFrame(out.sum(axis=0))
        len(centroids) == sum(sums[0])
        
        out.to_csv(args.output, sep=",")
        dups.close()

if __name__ == "__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("--output","-o",
                help = "The matrix file to write to")
        parser.add_argument("--dir","-d",
                help = "Directory containing all the files to combine")
        parser.add_argument("--ext","-e",
                help = "Extension of the files to combine")
        parser.add_argument("--duplicates","-p",
                help = "The file that contains duplicate entries found")
        args = parser.parse_args()
        
        main(args)
