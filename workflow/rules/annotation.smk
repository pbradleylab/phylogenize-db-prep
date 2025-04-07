def get_initial_queries(wildcards):
    """
    Return the initial full set of query sequences
    """
    query = config["annotation"]["mapping_databases"].get(wildcards.mapping_db)
    if query:
        return query
    else:
        raise ValueError(f"Query not found for database {wildcards.mapping_db}")

def get_cumulative_unmapped_queries(wildcards):
    """
    Generate the cumulative unmapped sequences across all previous target databases
    """
    # Track the sequence of target databases from configuration
    target_dbs = list(config["target_db"]["paths"].keys())
    
    # Find the current target database's index
    current_target_index = target_dbs.index(wildcards.target_db)
    
    # If it's the first target database, use the initial full query
    if current_target_index == 0:
        return get_initial_queries(wildcards)
    
    # Construct the path to the cumulative unmapped sequences
    previous_targets = target_dbs[:current_target_index]
    
    # Generate a list of unmapped sequence files from previous iterations
    unmapped_files = expand(
        "results/{{database}}/annotation/faSomeRecords/cumulative_unmapped/{mapping_db}_after_{target_db}.fa",
        database=wildcards.database,
        mapping_db=wildcards.mapping_db,
        target_db=previous_targets
    )
    
    return unmapped_files

def get_targets(wildcards):
    target = config["target_db"]["paths"].get(wildcards.target_db)
    if target:
        return target
    else:
        raise ValueError(f"Target not found for database {wildcards.target_db}")

def get_cumulative_unmapped_queries(wildcards):
    """
    Generate the cumulative unmapped sequences across all previous target databases
    """
    # Track the sequence of target databases from configuration
    target_dbs = list(config["target_db"]["paths"].keys())
    
    # Get the initial query sequences
    query = config["annotation"]["mapping_databases"].get(wildcards.mapping_db)
    if not query:
        raise ValueError(f"Query not found for database {wildcards.mapping_db}")
    
    # If no target_db is specified (like during database creation), return initial query
    if not hasattr(wildcards, 'target_db'):
        return query
    
    # Find the current target database's index
    current_target_index = target_dbs.index(wildcards.target_db)
    
    # If it's the first target database, use the initial full query
    if current_target_index == 0:
        return query
    
    # Construct the path to the cumulative unmapped sequences
    previous_targets = target_dbs[:current_target_index]
    
    # Generate a list of unmapped sequence files from previous iterations
    unmapped_files = expand(
        "results/{{database}}/annotation/faSomeRecords/cumulative_unmapped/{mapping_db}_after_{target_db}.fa",
        database=wildcards.database,
        mapping_db=wildcards.mapping_db,
        target_db=previous_targets
    )
    
    return unmapped_files


rule make_query_databases:
    input: get_initial_queries
    output:
        index="results/{database}/annotations/make_query_databases/{mapping_db}/{mapping_db}.index",
        query_path=directory("results/{database}/annotations/make_query_databases/{mapping_db}/")
    conda: "../envs/clustering.yml"
    log: "logs/{database}/annotations/make_query_databases/{mapping_db}.log"
    threads: config["mmseqs2"]["createdb"]["threads"]
    shell:
        """
        mmseqs createdb {input} {output.query_path}/{wildcards.mapping_db} --dbtype 1 2> {log}
        mmseqs createindex {output.query_path}/{wildcards.mapping_db} /tmp 2> {log}
        """

rule make_target_databases:
    input: get_targets
    output:
        index="results/{database}/annotations/make_target_databases/{target_db}/{target_db}.index",
        target_path=directory("results/{database}/annotations/make_target_databases/{target_db}/")
    conda: "../envs/clustering.yml"
    log: "logs/{database}/annotations/make_target_databases/{target_db}.log"
    threads: config["mmseqs2"]["createdb"]["threads"]
    shell:
        """
        mmseqs createdb {input} {output.target_path}/{wildcards.target_db} --dbtype 1 2> {log}
        mmseqs createindex {output.target_path}/{wildcards.target_db} /tmp 2> {log}
        """

rule make_cumulative_unmapped_queries:
    input: get_cumulative_unmapped_queries
    output:
        cumulative_queries="results/{database}/annotations/make_query_databases/{mapping_db}/{mapping_db}_cumulative.fa"
    shell:
        """
        # If multiple files are input, concatenate them
        if [ $(echo "{input}" | wc -w) -gt 1 ]; then
            cat {input} > {output.cumulative_queries}
        else
            cp {input} {output.cumulative_queries}
        fi
        """

rule make_cumulative_query_databases:
    input: rules.make_cumulative_unmapped_queries.output.cumulative_queries
    output:
        index="results/{database}/annotations/make_cumulative_query_databases/{mapping_db}/{mapping_db}_cumulative.index",
        query_path=directory("results/{database}/annotations/make_cumulative_query_databases/{mapping_db}/")
    conda: "../envs/clustering.yml"
    log: "logs/{database}/annotations/make_cumulative_query_databases/{mapping_db}.log"
    threads: config["mmseqs2"]["createdb"]["threads"]
    shell:
        """
        mmseqs createdb {input} {output.query_path}/{wildcards.mapping_db}_cumulative --dbtype 1 2> {log}
        mmseqs createindex {output.query_path}/{wildcards.mapping_db}_cumulative /tmp 2> {log}
        """

rule map_query:
    input:
       target=rules.make_target_databases.output.target_path,
       query=rules.make_cumulative_query_databases.output.query_path
    output:
        outdir=directory("results/{database}/annotation/map_query/{target_db}/{mapping_db}/"),
        #index="results/{database}/annotation/map_query/{target_db}/{mapping_db}/{mapping_db}_map.index"
    params:
        sensitivity=config["mmseqs2"]["map"]["sensitivity"]
    conda: "../envs/clustering.yml"
    log: "logs/{database}/annotations/map_query/{target_db}_{mapping_db}.log"
    threads: config["mmseqs2"]["map"]["threads"]
    shell:
        """
        mmseqs search --threads {threads} {input.query}/{wildcards.mapping_db}_cumulative \
            {input.target}/{wildcards.target_db} {output.outdir}/{wildcards.mapping_db} \
            results/tmp50 --min-seq-id 0.50 {params.sensitivity} 2> {log}
        """

rule mmseqs2_convertalis:
     input:
         query=rules.make_cumulative_query_databases.output.query_path,
         target=rules.make_target_databases.output.target_path,
         map=rules.map_query.output.outdir
     output:
         blast="results/{database}/annotation/mmseqs2/convertalis/{target_db}_{mapping_db}_convertlis.8",
         list="results/{database}/annotation/mmseqs2/convertalis/{target_db}_{mapping_db}_convertlis.list"
     conda: "../envs/clustering.yml"
     log: "logs/{database}/annotation/mmseqs2_convertalis_blast/{target_db}_{mapping_db}.log"
     threads: config["mmseqs2"]["convertalis"]["threads"]
     shell:
         """
         mmseqs convertalis {input.query}/{wildcards.mapping_db}_cumulative \
             {input.target}/{wildcards.target_db} \
             {input.map}/{wildcards.mapping_db} {output.blast} --format-mode 4 \
             --format-output query,target,pident 2> {log}
         cut -f1 {output.blast} | sed '1d' > {output.list}
         """

rule get_top_50_evals:
     input: rules.mmseqs2_convertalis.output.blast
     output:
         unfiltered="results/{database}/annotation/mmseqs2/top_50/{target_db}_{mapping_db}_convertlis.tsv",
         tophits="results/{database}/annotation/mmseqs2/top_50/{target_db}_{mapping_db}_convertlis_tophits.tsv"
     params:
         outdir="results/{database}/annotation/mmseqs2/top_50/"
     shell:
         """
         touch {output.unfiltered}
         awk '$3>50 {{print}}' {input} > {output.unfiltered}
         python workflow/scripts/get_top_hits.py -i {output.unfiltered} -o {output.tophits}
         """

rule get_unaligned_sequences:
    input:
        aligned=rules.get_top_50_evals.output.unfiltered,
        all_sequences=rules.make_cumulative_unmapped_queries.output.cumulative_queries
    output:
        unmapped="results/{database}/annotation/faSomeRecords/unmapped/{target_db}_{mapping_db}.fa",
        cumulative_unmapped="results/{database}/annotation/faSomeRecords/cumulative_unmapped/{mapping_db}_after_{target_db}.fa"
    conda: "../envs/clustering.yml"
    shell:
        """
        cut -f1 {input.aligned} | sed '1d' | uniq > /tmp/tmp.aligned
        grep '>' {input.all_sequences} | sed "s/>//g" > /tmp/tmp.all
        grep -F -v -x -f /tmp/tmp.aligned /tmp/tmp.all > /tmp/tmp.unaligned
        
        # Current iteration's unmapped
        faSomeRecords {input.all_sequences} /tmp/tmp.unaligned {output.unmapped}
        
        # Cumulative unmapped for next iterations
        cp {output.unmapped} {output.cumulative_unmapped}
        """

