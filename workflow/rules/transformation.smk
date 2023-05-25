""" Add rules to this section that are related to quality control and post processing.
"""
from scripts.utils import *
configfile: "config/config.json"

def get_pangenomes(wildcards):
    pangenomes = get_subsample_attributes(wildcards.pangenome, "reads", pep)
    return pangenomes

def get_mmseq2_input(wildcards):
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

# Map the protein sequences to 
rule create_mmseq2_query_db:
    input: get_mmseq2_input
    output: "resources/{database}.tsv"
    params:
        db_name="resources/{database}"
    conda: "../envs/transformation.yml"
    shell:
        """
        mmseqs createdb {input} {params.db_name} --dbtype 1
        mmseqs createtsv {params.db_name} {params.db_name} {output}
        """

rule download_UniRef90:
     input: get_mmseq2_input
    output: "resources/{database}.tsv"
    params:
        db_name="resources/{database}"
    conda: "../envs/transformation.yml"
    shell:
        """
        mmseqs databases UniRef90/Swiss-Prot resources/UniRef90/swissprot tmp
        """   

rule create_mmseq2_target_db:
     input: rules.download_UniRef90.output
    output: "resources/UniRef90.tsv"
    params:
        db_name="resources/UniRef90"
    conda: "../envs/transformation.yml"
    shell:
        """
        mmseqs createdb {input} {params.db_name} --dbtype 1
        mmseqs createtsv {params.db_name} {params.db_name} {output}
        """

rule create_mmseq2_target_index:
     input: rules.download_UniRef90.output
    output: "resources/UniRef90.idx"
    params:
        db_name="resources/UniRef90"
    conda: "../envs/transformation.yml"
    shell:
        """
        mmseqs createindex {params.db_name} tmp
        """

rule mmseq2_query:
    input: 
        rules.create_mmseq2_query_db.output,
        rules.create_mmseq2_target_index.output
    output: "results/{database}/mmseq2/{database}.m8"
    params:
        query="resources/{database}",
        target="resources/UniRef90",
        results_db="results/mmseq2/db/{database}"
    shell:
        """
        mmseqs search {params.query} {params.target} {params.results_db} tmp
        """

# rule mmseq2_convert_blast:
#     input: rules.mmseq2_query.output
#     output: "results/{database}/mmseq2/{database}.tsv"
#     params:
#         query="resources/{database}",
#         target="resources/UniRef90",
#         results_db="results/mmseq2/db/{database}"
#     shell:
#         """
#         mmseqs convertalis {params.query} {params.target} {params.results_db} {input}
#         """

# rule mmseq2_cluster:
#     input: rules.mmseq2_query.output
#     output: 
#         database="results/{database}/mmseq2/{database}.dbtype"
#     params:
#         db_name="resources/{database}",
#         out_dir="results/{database}/mmseq2/{database}",
#         seq_id_precent=config["mmseq2"]["seq_id_precent"],
#         tmp_dir=config["mmseq2"]["tmp_dir"]
#     conda: "../envs/transformation.yml"
#     threads: config["mmseq2"]["threads"]
#     shell:
#         """
#         mmseqs taxonomy {params.db_name}
#         mmseqs cluster {params.db_name} {params.out_dir} {params.tmp_dir} --min-seq-id {params.seq_id_precent} --threads {threads}
#         """

# rule midas:
#     input:rules.mmseq2.output
#     output: 
#     conda: "../envs/transformation.yml"
#     shell:
#         """
#         """

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