""" Resources retrieval and any rule relating to the transformation of resource files    
should be placed here. 
"""
from scripts.utils import *
configfile: "config/config.json"
include: "translation.smk"

def get_queries(wildcards):
    query = config["annotation"]["mapping_databases"].get(wildcards.mapping_db)
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
    query = config["annotation"]["mapping_databases"].get(wildcards.mapping_db)
    if not query:
        raise ValueError(f"Query not found for database {wildcards.mapping_db}")

    # Find the current target database's index
    current_target_index = TARGET_DBS.index(wildcards.target_db)

    # If it's the first target database, use the initial full query
    if current_target_index == 0:
        return query

    # Otherwise, use the unmapped sequences from the previous database
    previous_target = TARGET_DBS[current_target_index - 1]
    return f"results/{wildcards.database}/clustering/faSomeRecords/cumulative_unmapped/{wildcards.mapping_db}_after_{previous_target}.fa"


# Process input sequences for current iteration
rule prepare_current_iteration_input:
    input: get_cumulative_unmapped_queries
    output:
        current_input="results/{database}/annotations/current_input/{target_db}/{mapping_db}_input.fa"
    shell:
        """
        case {input} in
            *.fa.gz|*.fasta.gz|*.faa.gz) gunzip -c {input} > {output.current_input} ;;
            *) cp {input} {output.current_input} ;;
        esac
        """

# Creates a query database, query being the database containing the
# pangenomes that is being made into a final species level protein
# binary for Phylogenize.
# Create databases from queries
rule make_query_databases:
    input: get_queries
    output:
        index="results/{database}/make_query_databases/{mapping_db}/{mapping_db}.index",
        query_path=directory("results/{database}/make_query_databases/{mapping_db}/")
    conda: "../envs/query.yml"
    log: "logs/{database}/make_query_databases/{mapping_db}.log"
    threads: config["mmseqs2"]["createdb"]["threads"]
    shell:
        """
        mmseqs createdb {input} {output.query_path}/{wildcards.mapping_db} --dbtype 1 2> {log}
        mmseqs createindex {output.query_path}/{wildcards.mapping_db} /tmp 2> {log}
        """

# Ensures that for each target the mapping database is fully run 
checkpoint database_processing_checkpoint:
    input:
        db_results = lambda wildcards: expand(
            "results/{database}/clustering/faSomeRecords/cumulative_unmapped/{mapping_db}_after_{target_db}.fa",
            database=wildcards.database,
            mapping_db=config["annotation"]["mapping_databases"].keys(),
            target_db=wildcards.target_db
        )
    output:touch("results/{database}/checkpoints/{target_db}_processed.done")
    shell:
        """
        echo "Processing of target database {wildcards.target_db} is complete."
        """

# Create mmseqs database from input for current iteration
rule make_current_query_databases:
    input: rules.prepare_current_iteration_input.output.current_input
    output:
        index="results/{database}/make_current_query_databases/{target_db}/{mapping_db}/{mapping_db}.index",
        query_path=directory("results/{database}/make_current_query_databases/{target_db}/{mapping_db}/")
    conda: "../envs/query.yml"
    log: "logs/{database}/make_current_query_databases/{target_db}_{mapping_db}.log"
    threads: config["mmseqs2"]["createdb"]["threads"]
    shell:
        """
        mmseqs createdb {input} {output.query_path}/{wildcards.mapping_db} --dbtype 1 2> {log}
        mmseqs createindex {output.query_path}/{wildcards.mapping_db} /tmp 2> {log}
        """
