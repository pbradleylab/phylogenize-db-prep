""" Add rules to this section that are related to quality control and post processing.
"""
from scripts.utils import *
configfile: "config/config.json"

def get_pangenomes(wildcards):
    pangenomes = get_subsample_attributes(wildcards.pangenome, "reads", pep)
    return pangenomes

def get_mmseqs2_input(wildcards):
    outputLST = []
    for subsample in pep.subsample_table.subsample.tolist():
        project = get_subsample_attributes(subsample, "project", pep)
        # Always run rules on the outside
        outputLST.append(rules.transeq.output[0].format(database=project, pangenome=subsample))
    return outputLST

# Translate nucleotides per genome to peptide sequences per genome.
rule transeq:
    input:get_pangenomes
    output: "results/{database}/transeq/{pangenome}.ffn"
    conda: "../envs/transformation.yml"
    shell:
        """
        transeq {input} {output}
        """

rule create_mmseqs2_query_db:
    input: get_mmseqs2_input
    output:
        index="resources/{database}/custom/custom.index",
        db_prefix=directory("resources/{database}/custom")
    params:
        query_prefix="custom"
    conda: "../envs/transformation.yml"
    log: "logs/{database}/custom/{database}.log"
    shell:
        """
        mmseqs createdb {input} {output.db_prefix}/{params.query_prefix} --dbtype 1 2> {log}
        mmseqs createindex {output.db_prefix}/{params.query_prefix} /tmp 2> {log}
        """

rule create_mmseqs2_target_db:
    output: 
         fasta="resources/{database}/uniprot90/tmp/latest/uniref90.fasta.gz",
         db_prefix="resources/{database}/uniprot90"
    params:
        target_prefix="UniRef90"
    conda: "../envs/transformation.yml"
    threads: 32
    shell:
        """
        mmseqs databases UniRef90 {params.target_prefix} {output.db_prefix}/{params.target_prefix} --threads {threads}
        """

# Map the amino acid sequences by similarity in the UniProt 90 database
# The internal prefilter module is called which is high sensitivity to 
# detect high scores and ungapped alignment. This could be exchanged for
# the `mmseqs search` command for lower sensitivity.
# 
# Note: We call params here from previous rules. This method is continued
#    throughout all subsequent methods for continuety. These can be abstracted
#    out to a config however it looses some of the automation that way.
rule mmseqs2_map:
     input:
         query=rules.create_mmseqs2_query_db.output.db_prefix,
         target=rules.create_mmseqs2_target_db.output.db_prefix
     output:
         out_dir=directory("results/{database}/mmseqs2/mapping/"),
         index="results/{database}/mmseqs2/mapping/{database}_map.index"
     log: "logs/{database}/mmseqs2/mapping/mmseqs2_map.log"
     params:
         prefix="{database}_map",
         query_prefix=rules.create_mmseqs2_query_db.params.query_prefix,
         target_prefix=rules.create_mmseqs2_target_db.params.target_prefix
     threads: config["mmseqs2"]["threads"]
     conda: "../envs/transformation.yml"
     shell:
         """
         mkdir -p {output.out_dir}
         mmseqs map --threads {threads} {input.query}/{params.query_prefix} {input.target}/{params.target_prefix} {output.out_dir}/{params.prefix} /tmp/tmp -a 2> {log}
         """

# Converts the database's mappings to a sam format. The unmapped (unaligned)
# sequences are then taken to generate a new database in rule:
# `samtools_get_aligned`.
rule mmseqs2_convertalis:
     input:
         query=rules.create_mmseqs2_query_db.output.db_prefix,
         target=rules.create_mmseqs2_target_db.output.db_prefix,
         mapped=rules.mmseqs2_map.output.out_dir
     output: "results/{database}/mmseqs2/convertalis/{database}_convertlis.sam"
     params:
         prefix=rules.mmseqs2_map.params.prefix,
         query_prefix=rules.create_mmseqs2_query_db.params.query_prefix,
         target_prefix=rules.create_mmseqs2_target_db.params.target_prefix
     threads: config["mmseqs2"]["threads"]
     conda: "../envs/transformation.yml"
     shell:
         """
         mmseqs convertalis {input.query}/{params.query_prefix} {input.target}/{params.target_prefix} {input.mapped}/{params.prefix} {output} --format-mode 1
         """

# Get all the sequences that are not labeled as unaligned in the .sam 
# file. The output is a bam as to retain headers as samtools doesn't 
# keep header information which results in many downstream errors for
# programs.
rule samtools_get_unaligned:
     input:
         unmapped=rules.mmseqs2_convertalis.output
     output: "results/{database}/samtools/unaligned/{database}_aligned.bam"
     log: "logs/{database}/samtools/mapping/{database}_map.log"
     conda: "../envs/transformation.yml"
     shell:
         """
         samtools view -b -f 4 {input} -o {output}
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
        mmseqs createindex {output.out_dir}/{params.unaligned_prefix} /tmp 2> {log}
        """

# Cluster the unaligned protein sequence database from 
# mmseqs' search command.
rule mmseqs2_linclust:
     input: rules.create_mmseqs2_unaligned_db.output.out_dir
     output: 
         database="results/{database}/mmseqs2/linclust/unaligned_linclust.index",
         out_dir=directory("results/{database}/mmseqs2/linclust/")
     params:
         unaligned_prefix=rules.create_mmseqs2_unaligned_db.params.unaligned_prefix,
         prefix="unaligned_linclust",
         seq_id_precent=config["mmseqs2"]["seq_id_precent"],
         tmp_dir=config["mmseqs2"]["tmp_dir"]
     conda: "../envs/transformation.yml"
     log: "logs/{database}/mmseqs2/linclust/mmseqs2_linclust.log"
     threads: config["mmseqs2"]["threads"]
     shell:
         """
         mmseqs linclust {input}/{params.unaligned_prefix} {output.out_dir}/{params.prefix} {params.tmp_dir} --min-seq-id {params.seq_id_precent} --threads {threads} 2> {log}
         """

# Add taxonomy to the database. Mmseqs2 only uses uniprot internally,
# therefore since ids may be from an assortment of databases we assume
# the user can suply a mapping file as explained in the readme prior 
# to running this workflow. Uniprot90 ids are used by default to match
# mmseqs2.
#
# Generate a taxonomy database from the clustered database of unmapped 
# sequences representing a single species.

# Generate a maxtrix for all of the databases made with the taxon id
# on the top and the protein family on the y axis.

# rule peptide_matrix_generation:
#     input:rules.build_pangenome_database.output
#     output: 
#     params:
#     log: 
#     resources: 
#     conda: "../envs/transformation.yml"
#     shell:
#         """
#         """
