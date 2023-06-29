""" Add rules to this section that are related to quality control and post processing.
"""
from scripts.utils import *
configfile: "config/config.json"

include: "resources.smk"
include: "mapping.smk"

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

# Converts the database's mappings to a sam format. The unmapped (unaligned)
# sequences are then taken to generate a new database in rule:
# `samtools_get_aligned`.
rule mmseqs2_convertalis_sam:
     input:
         query=rules.create_mmseqs2_query_db.output.query_path,
         target=rules.create_mmseqs2_target_db.output.uniprot90_path,
         mapped=rules.mmseqs2_map.output.out_dir
     output: "results/{database}/mmseqs2/convertalis/{database}_convertlis.sam"
     params:
         prefix=rules.mmseqs2_map.params.prefix,
         query_prefix=rules.create_mmseqs2_query_db.params.query_prefix,
         target_prefix=rules.create_mmseqs2_target_db.params.uniprot90_prefix
     threads: config["mmseqs2"]["convertalis"]["threads"]
     conda: "../envs/transformation.yml"
     shell:
         """
         mmseqs convertalis \
             {input.query}/{params.query_prefix} \
             {input.target}/{params.target_prefix} \
             {input.mapped}/{params.prefix} \
             {output} --format-mode 1
         """


rule mmseqs2_convertalis_blast:
     input:
         query=rules.create_mmseqs2_query_db.output.query_path,
         target=rules.create_mmseqs2_target_db.output.uniprot90_path,
         map=rules.mmseqs2_map.output.out_dir
     output: "results/{database}/mmseqs2/convertalis/{database}_convertlis.8"
     params:
         prefix=rules.mmseqs2_map.params.prefix,
         query_prefix=rules.create_mmseqs2_query_db.params.query_prefix,
         target_prefix=rules.create_mmseqs2_target_db.params.uniprot90_prefix
     threads: config["mmseqs2"]["convertalis"]["threads"]
     conda: "../envs/transformation.yml"
     shell:
         """
         mmseqs convertalis {input.query}/{params.query_prefix} {input.target}/{params.target_prefix} {input.map}/{params.prefix} {output} --format-mode 4 --format-output query,target,pident
         """

rule get_top_90_evals:
     input: rules.mmseqs2_convertalis_blast.output
     output: "results/{database}/mmseqs2/top_90/{database}_convertlis.tsv"
     conda: "../envs/transformation.yml"
     shell:
         """
         awk '$3>90 {{print}}' {input} > {output}
         """

# Get all the sequences that are not labeled as unaligned in the .sam
# file. The output is a bam as to retain headers as samtools doesn't
# keep header information which results in many downstream errors for
# programs.
rule samtools_get_unaligned:
     input: rules.mmseqs2_convertalis_sam.output,
     output: "results/{database}/samtools/unaligned/{database}_unaligned.bam",
     log: "logs/{database}/samtools/mapping/{database}_map.log"
     conda: "../envs/transformation.yml"
     shell:
         """
         samtools view -b -f 4 {input} -o {output} 2> {log}
         """

# Convert the unaligned samples to a fasta to build a new database.
# We consider these to be temporary files as we retain the bam for QC
# afterwards since reads may have failed to align due to other reasons
# and not solely based on species sequence homology.
rule samtools_fasta:
     input: rules.samtools_get_unaligned.output
     output: "results/{database}/samtools/fasta/{database}.fasta"
     conda: "../envs/transformation.yml"
     shell:
         """
         samtools fasta {input} > {output}
         """

# Create a new database that is declared as temporary. This database
# holds the unaligned peptide sequences that are assumed as potential
# species specific alignments.
rule create_mmseqs2_unaligned_db:
    input: rules.samtools_fasta.output
    output:
        out_dir=directory("resources/{database}/mmseqs2/unmapped/"),
        index="resources/{database}/mmseqs2/unmapped/unmapped.index",
    params:
        unaligned_prefix="unmapped"
    conda: "../envs/transformation.yml"
    log: "logs/{database}/mmseqs2/create_mmseqs2_unaligned/mmseqs2_create_mmseqs2_unaligned.log"
    shell:
        """
        mkdir -p {output.out_dir}
        mmseqs createdb {input} {output.out_dir}/{params.unaligned_prefix} --dbtype 1 2> {log}
        mmseqs createindex {output.out_dir}/{params.unaligned_prefix} \
            /tmp 2> {log}
        """

# Combines species with a 90% or greater identity match to the target database, 
# and the unmapped regions to a list of species specific vectors by their centroid.
rule combine_species_hits:
    input:
        unmapped=rules.samtools_get_unaligned.output,
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
