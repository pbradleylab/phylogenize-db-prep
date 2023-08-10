include: "clustering.smk"
include: "mapping.smk"


rule mmseqs2_convertalis_blast_uniref50_db:
     input:
         query=rules.create_mmseqs2_query_db.output.query_path,
         target=rules.create_uniref50.output.uniref50_path,
         map=rules.mmseqs2_map_uniref50.output.outdir
     output: "results/{database}/uniref50/mmseqs2/convertalis/{database}_convertlis.8"
     params:
         prefix=rules.mmseqs2_map_uniref50.params.prefix,
         query_prefix=rules.create_mmseqs2_query_db.params.query_prefix,
         target_prefix=rules.create_uniref50.params.uniref50_prefix
     threads: config["mmseqs2"]["convertalis"]["threads"]
     conda: "../envs/blast.yml"
     shell:
         """
         mmseqs convertalis {input.query}/{params.query_prefix} \
             {input.target}/{params.target_prefix} \
             {input.map}/{params.prefix} {output} --format-mode 4 \
             --format-output query,target,pident
         """

rule mmseqs2_convertalis_blast_uhgp50_db:
     input:
         query=rules.create_mmseqs2_query_db.output.query_path,
         target=rules.create_uhgp50.output,
         map=rules.mmseqs2_map_uhgp50.output.outdir
     output: "results/{database}/uhgp50/mmseqs2/convertalis/{database}_convertlis.8"
     params:
         prefix=rules.mmseqs2_map_uhgp50.params.prefix,
         query_prefix=rules.create_mmseqs2_query_db.params.query_prefix,
         target_prefix=rules.create_uhgp50.output
     threads: config["mmseqs2"]["convertalis"]["threads"]
     conda: "../envs/blast.yml"
     shell:
         """
         mmseqs convertalis {input.query}/{params.query_prefix} \
             {input.target}/uhgp50 \
             {input.map}/{params.prefix} {output} --format-mode 4 \
             --format-output query,target,pident
         """

rule mmseqs2_convertalis_unmapped_blast_uhgp50_db:
     input:
         query=rules.create_mmseqs2_unaligned_uhgp50_db.output.outdir,
         target=rules.create_mmseqs2_unaligned_uhgp50_db.output.outdir,
         cluster=rules.mmseqs2_linclust_uhgp50_db.output.outdir
     output: "results/{database}/uhgp50/mmseqs2/convertalis/unaligned/{database}_convertlis.8"
     params:
         prefix=rules.mmseqs2_linclust_uhgp50_db.params.prefix,
         query_prefix=rules.create_mmseqs2_unaligned_uhgp50_db.params.unaligned_prefix,
         target_prefix=rules.create_mmseqs2_unaligned_uhgp50_db.params.unaligned_prefix
     threads: config["mmseqs2"]["convertalis"]["threads"]
     conda: "../envs/blast.yml"
     shell:
         """
         mmseqs createtsv {input.query}/{params.query_prefix} \
             {input.query}/{params.query_prefix} \
             {input.cluster}/{params.prefix} {output}
         """
