include: "clustering.smk"
include: "mapping.smk"

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
