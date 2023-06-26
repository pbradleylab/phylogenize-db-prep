# Map the amino acid sequences by similarity in the UniProt 90 database
# The internal prefilter module is called which is high sensitivity to
# detect high scores and ungapped alignment. This could be exchanged for
# the `mmseqs search` command for lower sensitivity.
#
# Note: We call params here from previous rules. This method is continued
#    throughout all subsequent methods for continuety. These can be abstracted
#    out to a config however it looses some of the automation that way.
rule mmseqs2_map:
     input:
         query=rules.create_mmseqs2_query_db.output.query_path,
         target=rules.create_mmseqs2_target_db.output.uniprot90_path
     output:
         out_dir=directory("results/{database}/mmseqs2/mapping/"),
         index="results/{database}/mmseqs2/mapping/{database}_map.index"
     log: "logs/{database}/mmseqs2/mapping/mmseqs2_map.log"
     params:
         prefix="{database}_map",
         query_prefix=rules.create_mmseqs2_query_db.params.query_prefix,
         target_prefix=rules.create_mmseqs2_target_db.params.uniprot90_prefix
     threads: config["mmseqs2"]["map"]["threads"]
     conda: "../envs/transformation.yml"
     shell:
         """
         mkdir -p {output.out_dir}
            mmseqs map --threads {threads} {input.query}/{params.query_prefix} \
            {input.target}/{params.target_prefix} {output.out_dir}/{params.prefix} \
            /tmp/tmp -a 2> {log}
         """