include: "query_db.smk"
include: "target_db.smk"
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
rule mmseqs2_map_uniref50:
     input:
         query=rules.create_mmseqs2_query_db.output.query_path,
         target=rules.create_uniref50.output.uniref50_path
     output:
         outdir=directory("results/{database}/uniref50/mmseqs2/mapping/"),
         index="results/{database}/uniref50/mmseqs2/mapping/{database}_map.index"
     log: "logs/{database}/uniref50/mmseqs2/mapping/mmseqs2_map.log"
     params:
         prefix="{database}_map",
         query_prefix=rules.create_mmseqs2_query_db.params.query_prefix,
         target_prefix="UniRef50",
         sensitivity=config["mmseqs2"]["map"]["sensitivity"]
     threads: config["mmseqs2"]["map"]["threads"]
     conda: "../envs/mapping.yml"
     shell:
         """
         mmseqs search --threads {threads} {input.query}/{params.query_prefix} \
            {input.target}/{params.target_prefix} {output.outdir}/{params.prefix} \
            results/tmp50 --min-seq-id 0.50 {params.sensitivity} 2> {log}
         """
