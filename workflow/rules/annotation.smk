def get_queries(wildcards):
    query=config["annotation"]["mapping_databases"].get(wildcards.mapping_db)
    if query:
        return query
    else:
        raise ValueError(f"Query not found for database {wildcards.database}")

def get_targets(wildcards):
    target=config["target_db"]["paths"].get(wildcards.target_db)
    if target:
        return target
    else:
        raise ValueError(f"Target not found for database {wildcards.database}")

def get_next_queries(wildcards):
    # Wait for the checkpoint to finish
    checkpoint_output = checkpoints.get_unaligned_sequences.get(**wildcards).output.unaligned
    
    # If the unaligned file is empty, stop iterating
    if os.path.exists(checkpoint_output) and os.stat(checkpoint_output).st_size > 0:
        return checkpoint_output
    else:
        raise WorkflowError(f"No more unaligned sequences for {wildcards.database}")


rule make_query_databases:
    input: lambda wildcards: get_next_queries(wildcards)
    output:
        index="results/{database}/annotations/make_query_databases/{mapping_db}/{mapping_db}.index",
        query_path=directory("results/{database}/annotations/make_query_databases/{mapping_db}/")
    conda: "../envs/annotations.yml"
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
        query_path=directory("results/{database}/annotations/make_target_databases/{target_db}/")
    conda: "../envs/annotations.yml"
    log: "logs/{database}/annotations/make_target_databases/{target_db}.log"
    threads: config["mmseqs2"]["createdb"]["threads"]
    shell:
        """
        mmseqs createdb {input} {output.query_path}/{wildcards.target_db} --dbtype 1 2> {log}
        mmseqs createindex {output.query_path}/{wildcards.target_db} /tmp 2> {log}
        """

rule map_query:
    input:
       target=rules.make_target_databases.output,
       query=rules.make_query_databases.output
    output:
        outdir=directory("results/{database}/annotation/map_query/{target_db}/{mapping_db}/"),
        index="results/{database}/annotation/map_query/{target_db}/{mapping_db}/{mapping_db}_map.index"
    params:
        sensitivity=config["mmseqs2"]["map"]["sensitivity"]
    conda: "../envs/annotations.yml"
    log: "logs/{database}/annotations/map_query/{target_db}_{mapping_db}.log"
    threads: config["mmseqs2"]["map"]["threads"]
    shell:
        """
        mmseqs search --threads {threads} {input.query}/{wildcards.mapping_db} \
            {input.target}/{wildcards.target_db} {output.outdir}/{wildcards.mapping_db} \
            results/tmp50 --min-seq-id 0.50 {params.sensitivity} 2> {log}
        """

rule mmseqs2_convertalis:
     input:
         query=rules.make_query_databases.output.query_path,
         target=rules.make_target_databases.output.uniref50_path,
         map=rules.map_query.output.outdir
     output:
         blast="results/{database}/annotation/mmseqs2/convertalis/{target_db}_{mapping_db}_convertlis.8",
         list="results/{database}/annotation/mmseqs2/convertalis/{target_db}_{mapping_db}_convertlis.list"
     params:
         prefix=rules.mmseqs2_map_uniref50.params.prefix,
         query_prefix=rules.create_mmseqs2_query_db.params.query_prefix,
         target_prefix=rules.create_uniref50.params.uniref50_prefix
     conda: "../envs/blast.yml"
     log: "logs/{database}/annotation/mmseqs2_convertalis_blast/{target_db}_{mapping_db}.log"
     threads: config["mmseqs2"]["convertalis"]["threads"]
     shell:
         """
         mmseqs convertalis {input.query}/{params.query_prefix} \
             {input.target}/{params.target_prefix} \
             {input.map}/{params.prefix} {output.blast} --format-mode 4 \
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

checkpoint get_unaligned_sequences:
    input:
        aligned=rules.get_top_50_evals.output.unfiltered,
        all_sequences=get_queries
    output:
        unaligned="results/{database}/annotation/faSomeRecords/unmapped/{target_db}_{mapping_db}.fa"
    conda: "../envs/clustering.yml"
    shell:
        """
        cut -f1 {input.aligned} | sed '1d' | uniq > /tmp/tmp.aligned
        grep '>' {input.all_sequences} | sed "s/>//g" > /tmp/tmp.all
        grep -F -v -x -f /tmp/tmp.aligned /tmp/tmp.all > /tmp/tmp.unaligned
        faSomeRecords {input.all_sequences} /tmp/tmp.unaligned {output}
        """

#rule get_unaligned_sequences:
#    input:
#        aligned=rules.get_top_50_evals.output.unfiltered,
#        all_sequences=get_queries
#    output: "results/{database}/annotation/faSomeRecords/unmapped/{database}.fa"
#    conda: "../envs/clustering.yml"
#    shell:
#        """
#        cut -f1 {input.aligned} | sed '1d' | uniq > /tmp/tmp.aligned
#        grep '>' {input.all_sequences} | sed "s/>//g" > /tmp/tmp.all
#        grep -F -v -x -f /tmp/tmp.aligned /tmp/tmp.all > /tmp/tmp.unaligned
#        faSomeRecords {input.all_sequences} /tmp/tmp.unaligned {output}
#        """


