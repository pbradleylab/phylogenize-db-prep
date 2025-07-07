""" Resources retrieval and any rule relating to the transformation of resource files    
should be placed here. 
"""
from scripts.utils import *
configfile: "config/config.json"
include: "translation.smk"

def get_mmseqs2_input(wildcards):
    outputLST = []
    for subsample in pep.subsample_table.subsample.tolist():
        project = get_subsample_attributes(subsample, "project", pep)
        # Always run rules on the outside
        outputLST.append(rules.transeq.output[0].format(database=project, pangenome=subsample))
    return outputLST


# Creates a query database, query being the database containing the
# pangenomes that is being made into a final species level protein
# binary for Phylogenize.
rule create_mmseqs2_query_db:
    input: rules.reduce_complexity.output
    output:
        index="resources/{database}/custom/custom.index",
        query_path=directory("resources/{database}/custom")
    params:
        query_prefix="custom"
    conda: "../envs/query_db.yml"
    log: "logs/{database}/custom/{database}.log"
    threads: config["mmseqs2"]["createdb"]["threads"]
    shell:
        """
        mmseqs createdb {input} {output.query_path}/{params.query_prefix} --dbtype 1 2> {log}
        mmseqs createindex {output.query_path}/{params.query_prefix} /tmp 2> {log}
        """
