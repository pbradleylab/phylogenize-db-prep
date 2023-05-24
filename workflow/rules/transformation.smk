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

rule transeq:
    input:get_pangenomes
    output: "results/{database}/transeq/{pangenome}.ffn"
    conda: "../envs/transformation.yml"
    shell:
        """
        transeq {input} {output}
        """

rule create_mmseq2_db:
    input: get_mmseq2_input
    output: "resources/{database}.index"
    params:
        db_name="resources/{database}"    
    conda: "../envs/transformation.yml"
    shell:
        """
        mmseqs createdb {input} {params.db_name} --dbtype 1
        """

rule mmseq2:
    input: 
        get_mmseq2_input,
        rules.create_mmseq2_db.output
    output: 
        database="results/{database}/mmseq2/{database}.dbtype"
    params:
        db_name="resources/{database}",
        out_dir="results/{database}/mmseq2/{database}",
        seq_id_precent=0.9,
        tmp_dir="tmp90"
    conda: "../envs/transformation.yml"
    threads: 4
    shell:
        """
        mmseqs cluster {params.db_name} {params.out_dir} {params.tmp_dir} --min-seq-id {params.seq_id_precent} --threads {threads}
        """

# rule build_pangenome_database:
#     input:rules.90_clustering.output
#     output: 
#     params:
#     log: 
#     resources: 
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