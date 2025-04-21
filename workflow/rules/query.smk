""" Resources retrieval and any rule relating to the transformation of resource files    
should be placed here. 
"""
from scripts.utils import *
configfile: "config/config.json"
include: "translation.smk"

def get_queries(wildcards):
    query = config["files"]["fasta"].get(wildcards.mapping_db)
    if query:
        return query
    else:
        raise ValueError(f"Query not found for database {wildcards.mapping_db}")

def get_previous_target_db(wildcards):
    """Get the previous target database in the sequence"""
    try:
        idx = TARGET_DBS.index(wildcards.target_db)
        if idx > 0:
            return TARGET_DBS[idx-1]
        return None
    except ValueError:
        raise ValueError(f"Target database {wildcards.target_db} not found in config")

def get_cumulative_unmapped_queries(wildcards):
    """
    Generate the appropriate input based on the target database position in the sequence
    """
    # Get the initial query sequences
    query = config["files"]["fasta"].get(wildcards.mapping_db)
    current_target_index = TARGET_DBS.index(wildcards.target_db)
    if current_target_index == 0:
        return query

    # Otherwise, use the unmapped sequences from the previous database
    previous_target = TARGET_DBS[current_target_index - 1]
    return(rules.faSomeRecords.output.cumulative_unmapped.format(database=wildcards.database, mapping_db=wildcards.mapping_db, target_db=previous_target))


# Process input sequences for current iteration
rule prepare_current_iteration_input:
    input: get_cumulative_unmapped_queries
    output:"results/{database}/query/prepare_current_iteration_input/{target_db}/{mapping_db}_input.fa"
    shell:
        """
        case {input} in
            *.fa.gz|*.fasta.gz|*.faa.gz) gunzip -c {input} > {output} ;;
            *) cp {input} {output} ;;
        esac
        """

# Creates a query database, query being the database containing the
# pangenomes that is being made into a final species level protein
# binary for Phylogenize.
# Create databases from queries
rule make_query_databases:
    input: get_queries
    output:
        index="results/{database}/query/make_query_databases/{mapping_db}/{mapping_db}.index",
        query_path=directory("results/{database}/query/make_query_databases/{mapping_db}/")
    conda: "../envs/query.yml"
    log: "logs/{database}/query/make_query_databases/{mapping_db}.log"
    threads: config["mmseqs2"]["createdb"]["threads"]
    shell:
        """
        mmseqs createdb {input} {output.query_path}/{wildcards.mapping_db} --dbtype 1 2> {log}
        mmseqs createindex {output.query_path}/{wildcards.mapping_db} /tmp 2> {log}
        """

# Ensures that for each target the mapping database is fully run 
checkpoint database_processing_checkpoint:
    input:
        db_results=lambda wildcards: [
            rules.faSomeRecords.output.cumulative_unmapped.format(
                database=wildcards.database,
                mapping_db=mapping_db,
                target_db=wildcards.target_db)
            for mapping_db in config["files"]["fasta"].keys()]
    output: touch("results/{database}/checkpoints/database_processing_checkpoint/{mapping_db}_{target_db}_processed.done")
    shell:
        """
        echo "Processing of target database {wildcards.target_db} and query database {wildcards.mapping_db} is complete."
        """

# Create mmseqs database from input for current iteration
rule make_current_query_databases:
    input: rules.prepare_current_iteration_input.output
    output:
        index="results/{database}/query/make_current_query_databases/{target_db}/{mapping_db}/{mapping_db}.index",
        query_path=directory("results/{database}/query/make_current_query_databases/{target_db}/{mapping_db}/")
    conda: "../envs/query.yml"
    log: "logs/{database}/query/make_current_query_databases/{target_db}_{mapping_db}.log"
    threads: config["mmseqs2"]["createdb"]["threads"]
    shell:
        """
        mmseqs createdb {input} {output.query_path}/{wildcards.mapping_db} --dbtype 1 2> {log}
        mmseqs createindex {output.query_path}/{wildcards.mapping_db} /tmp 2> {log}
        """
