""" Add rules to this section that are related to quality control and post processing.
"""
from scripts.utils import *
configfile: "config/config.json"

include: "resources.smk"
include: "mapping.smk"
include: "database_management.smk"

def get_pangenomes(wildcards):
    pangenomes = get_subsample_attributes(wildcards.pangenome, "pangenomes", pep)
    return pangenomes


# Translate nucleotides per genome to peptide sequences per genome.
# Please check the config to set if stop codons shouldn't convert
# from the default '*' character to an 'X' representing any animo acid.
rule transeq:
    input:get_pangenomes
    output: "results/{database}/transeq/{pangenome}.ffn"
    params:
        clean=config["transeq"]["convert_missing_to_x"]
    log: "logs/{database}/transeq/{pangenome}.log"
    conda: "../envs/transformation.yml"
    shell:
        """
        transeq {input} {output} -clean {params.clean} 2> {log}
        """

# Get all the sequences that are not labeled as unaligned in the .sam
# file. The output is a bam as to retain headers as samtools doesn't
# keep header information which results in many downstream errors for
# programs.
rule samtools_aligned_fasta:
     input: rules.mmseqs2_convertalis_sam.output,
     output: "results/{database}/samtools/unaligned/{database}_unaligned.bam",
     log: "logs/{database}/samtools/mapping/{database}_map.log"
     conda: "../envs/transformation.yml"
     shell:
         """
         samtools view -b {input} -o {output} 2> {log}
         """

# Convert the unaligned samples to a fasta to build a new database.
# We consider these to be temporary files as we retain the bam for QC
# afterwards since reads may have failed to align due to other reasons
# and not solely based on species sequence homology.
rule samtools_fasta:
     input: rules.samtools_aligned_fasta.output
     output: "results/{database}/samtools/fasta/{database}.fasta"
     conda: "../envs/transformation.yml"
     shell:
         """
         samtools fasta {input} > {output}
         """

# Combines species with a 90% or greater identity match to the target database, 
# and the unmapped regions to a list of species specific vectors by their centroid.
rule combine_species_hits:
    input:
        unmapped=rules.mmseqs2_convertalis_unmapped_blast.output,
        identity_90=rules.get_top_90_evals.output
    output: 
         txt="results/{database}/mmseqs2/combined_species_hits/{database}.txt",
         out_dir=directory("results/{database}/mmseqs2/combined_species_hits/")
    conda: "../envs/transformation.yml"
    log: "logs/{database}/mmseqs2/hits_90/mmseqs2_hits_90.log"
    shell:
        """
        samtools view {input.unmapped} | cut -f1 > {output.txt}
        cat {input.identity_90} | cut -f1 >> {output.txt}
        """

rule create_species_matrix:
    input: rules.combine_species_hits.output.out_dir
    output: "results/{database}/final/species_matrix/{database}.txt"
    conda: "../envs/transformation.yml"
    log: "logs/{database}/mmseqs2/hits_90/mmseqs2_hits_90.log"
    shell:
        """
        python workflow/scripts/combine_species.py --output {output} --dir {input} --ext ".txt"
        """
