""" This script is for combining species into a larger matrix using pandas. The 
matrix  has centroids on the Y axis and the species accession on the top. Values
are a binary for if the centroid is present for the different predefined species.

In our case, this is used for the nucleotide clustered pangenomes for midas' GTDB 
database. 
"""
import os, re, sys
import pandas as pd
import argparse
import numpy as np
from dask import dataframe as dd 
import dask

def clean_proteins(lst):
        centroids=[]
        for x in lst:
            if "GCA" in x or "GCF" in x:
                x=x.split("_")
                x="".join(x[1:])
            if x.startswith("UniRef50"):
                x=x.replace("UniRef50_", "")
            centroids.append(x)
        return centroids


def main(args):
        centroids_dict={}
        binary=[]
        
        # Duplicate file
        dups=open(args.duplicates, "w")
        # Get pair and unpair files in separate lists
        frame = pd.read_csv(args.input, delimiter='\t', header=None, names=["species", "protein"])
        # Read in species and proteins lists
        #proteins = pd.read_csv(args.protein, header=None, names=["protein"], sep='\t')
        # Collect the species that are relavant for the binary if a subset has been selected, this works as well. 
        frame["species"] = [x.split("_")[0] for x in list(frame["species"])]
        frame["protein"] = clean_proteins(list(frame["protein"]))
        #print("Made Centroid list")
        frame["presence"] = 1
        print(frame)
         
        # Collect species we are making the binary for
        #proteins["protein"] = clean_proteins(list(proteins["protein"]))
        #frame=pd.merge(frame, proteins, on="protein", how="outer", indicator=True)
        # Make a template frame for all the proteins not found for any given species
        # and then concat them to all those that are shared
        #template=frame[frame["_merge"] == "right_only"]
        #template["_merge"] = False
        #frame=frame[frame["_merge"] == "both"]
        #frame["_merge"] = True
        #print(frame)
        #for species in list(set(frame["species"])):
        #     template["species"]=species
        #     frame=pd.concat([frame, template])
        #     print("species: "+species+" added.")
        #frame = frame.rename({"_merge":"presence"})
        #frame = frame.drop_duplicates()
        #frame["species"]=frame["species"].astype('int8')
        #frame["protein"]=frame["protein"].astype('string[pyarrow]')
        frame.reset_index(drop=True).to_feather(args.output) 

if __name__ == "__main__":
        parser = argparse.ArgumentParser()
        parser.add_argument("--output","-o",
                help = "The matrix file to write to")
        parser.add_argument("--input","-i",
                help = "Directory containing all the files to combine")
        parser.add_argument("--duplicates","-d",
                help = "The file that contains duplicate entries found")
        args = parser.parse_args()
        
        main(args)
