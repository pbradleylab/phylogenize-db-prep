""" Add rules to this section that are related to quality control and post processing.
"""
from scripts.utils import *
configfile: "config/config.json"

def get_pangenomes(wildwards):
    return get_pangenome_attributes(wildcards.sample, "sample_name", pep)

rule convert_to_neucleotides:
    input:get_pangenomes()
    output: 
    params:
    log: 
    resources: 
    conda: "../envs/transformation.yml"
    shell:
        """
        """

rule 90_clustering:
    input:rules.convert_to_neucleotides.output
    output: 
    params:
    log: 
    resources: 
    conda: "../envs/transformation.yml"
    shell:
        """
        """

rule build_pangenom_database:
    input:rules.90_clustering.output
    output: 
    params:
    log: 
    resources: 
    conda: "../envs/transformation.yml"
    shell:
        """
        """

rule peptide_matrix_generation:
    input:rules.build_pangenom_database.output
    output: 
    params:
    log: 
    resources: 
    conda: "../envs/transformation.yml"
    shell:
        """
        """