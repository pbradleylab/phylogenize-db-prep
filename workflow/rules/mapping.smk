include: "query.smk"
include: "target.smk"
include: "translation.smk"

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
rule map_query:
    input:
       target=rules.make_targets.output.target_path,
       query=rules.make_current_query_databases.output.query_path,
       # If this is not the first database, wait for previous database processing to complete
       previous_checkpoint=lambda wildcards: (
           f"results/{wildcards.database}/annotation/checkpoints/{get_previous_target_db(wildcards)}_processed.done"
           if get_previous_target_db(wildcards) is not None else []
       )
    output:
        outdir=directory("results/{database}/map_query/{target_db}/{mapping_db}/"),
        index="results/{database}/map_query/{target_db}/{mapping_db}/{mapping_db}_map.index"
    params:
        sensitivity=config["mmseqs2"]["map"]["sensitivity"]
    conda: "../envs/mapping.yml"
    log: "logs/{database}/map_query/{target_db}_{mapping_db}.log"
    threads: config["mmseqs2"]["map"]["threads"]
    shell:
        """
        mmseqs search --threads {threads} {input.query}/{wildcards.mapping_db} \
            {input.target}/{wildcards.target_db} {output.outdir}/{wildcards.mapping_db} \
            results/tmp50 --min-seq-id 0.50 {params.sensitivity} 2> {log}
        """
