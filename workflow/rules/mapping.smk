include: "query.smk"
include: "target.smk"
include: "translation.smk"


def get_iterative_mapping(wildcards):
    iterative_bool=config["mapping"]["iterative_mapping"]
    if iterative_bool == "False":
        out=rules.mmseqs2_search.output[0].format(
            database=wildcards.database,
            mapping_db=wildcards.mapping_db,
            target_db=wildcards.target_db)
    else:
        out=rules.mmseqs2_merge_results.output[0].format(
            database=wildcards.database,
            mapping_db=wildcards.mapping_db,
            target_db=wildcards.target_db)
    return(out)

rule mmseqs2_search_1:
    input:
        target=rules.make_targets.output.target_path,
        query=rules.make_current_query_databases.output.query_path
    output: "results/{database}/mapping/mmseqs2_search_1/{target_db}/{mapping_db}/{mapping_db}"
    conda: "../envs/mapping.yml"
    threads: config["mmseqs2"]["map"]["threads"]
    log: "logs/{database}/mapping/search_first/{target_db}_{mapping_db}.log"
    shell:
        """
        mmseqs search {input.query} {input.target} {output} results/tmp/mmseqs2_search_first/round1 -s 1 --threads {threads} 2> {log}
        """

checkpoint mmseqs2_nohits:
    input:
        query=rules.make_current_query_databases.output.query_path,
        target=rules.make_targets.output.target_path,
        results=rules.mmseqs2_search_1.output
    output: "results/{database}/mapping/mmseqs2_nohit_{round,[2-9][0-9]*}/{target_db}/{mapping_db}/{mapping_db}"
    conda: "../envs/mapping.yml"
    shell:
        """
        mmseqs listseqs {input.query} > all_queries.txt
        mmseqs result2flat {input.query} {input.target} {input.results} hits.tsv
        cut -f1 hits.tsv | sort -u > matched.txt
        comm -23 <(sort all_queries.txt) matched.txt > nohit.txt
        mmseqs createsubdb nohit.txt {input.query} {output}
        """

rule mmseqs2_search_iter:
    input:
        query=rules.mmseqs2_nohits.output,
        target=rules.make_targets.output.target_path
    output: "results/{database}/mapping/mmseqs2_search_{round}/{target_db}/{mapping_db}/{mapping_db}"
    conda: "../envs/mapping.yml"
    threads: config["mmseqs2"]["map"]["threads"]
    params:
        s=lambda wildcards: int(wildcards.round)
    shell:
        """
        mmseqs search {input.query} {input.target} {output} tmp/round{wildcards.round} -s {params.s} --threads {threads}
        """

rule mmseqs2_merge_results:
    input:
        query=rules.make_current_query_databases.output.query_path, 
        target=rules.make_targets.output.target_path,
        iters = expand(
                rules.mmseqs2_search_iter.output,
                target_db="{target_db}",
                mapping_db="{mapping_db}",
                database="{database}",
                round=range(1, config["mapping"]["max_iter"] + 1)),
    output: "results/{database}/mapping/mmseqs2_mergedbs/{target_db}/{mapping_db}/{mapping_db}"
    conda: "../envs/mapping.yml"
    shell:
        """
        mmseqs mergedbs {input.query} {input.target} {output} results/tmp_merge/ {input.iters}
        """

# Map the amino acid sequences by similarity in the uniref50 database
# The internal prefilter module is called which is high sensitivity to
# detect high scores and ungapped alignment. This could be exchanged for
# the `mmseqs search` command for lower sensitivity.
#
# Note: We call params here from previous rules. This method is continued
#    throughout all subsequent methods for continuety. These can be abstracted
#    out to a config however it looses some of the automation that way.
# Note: We use an iterative step method desribed in the mmseqs documentation
#    for how to speed up the mapping step iteratively. See the params set in 
#    the tools.json file.
# Map queries to target database
rule mmseqs2_search:
    input:
        target=rules.make_targets.output.target_path,
        query=rules.make_current_query_databases.output.query_path,
        previous=lambda wildcards: (
            rules.database_processing_checkpoint.output[0].format(
                    database=wildcards.database,
                    mapping_db=config["files"]["fasta"].keys(),
                    target_db=get_previous_target_db(wildcards)
                )
            if get_previous_target_db(wildcards) is not None else []
        )
    output:"results/{database}/mapping/map_query/{target_db}/{mapping_db}/{mapping_db}.index"
    params:
        outdir=directory("results/{database}/mapping/map_query/{target_db}/{mapping_db}/"),
        sensitivity=config["mmseqs2"]["map"]["sensitivity"],
        split=config["mmseqs2"]["map"]["split"],
	split_memory_limit=config["mmseqs2"]["map"]["split_memory_limit"],
        min_seq_id=config["mmseqs2"]["map"]["min_seq_id"]
    conda: "../envs/mapping.yml"
    log: "logs/{database}/mapping/map_query/{target_db}_{mapping_db}.log"
    threads: config["mmseqs2"]["map"]["threads"]
    shell:
        """
        mmseqs search --threads {threads} {input.query}/{wildcards.mapping_db} \
            {input.target}/{wildcards.target_db} {params.outdir}/{wildcards.mapping_db} \
            results/tmp50_2 {params.min_seq_id} {params.split} {params.sensitivity} {params.split_memory_limit} 2> {log}
        """
