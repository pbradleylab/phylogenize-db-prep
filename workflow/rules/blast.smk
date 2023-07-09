include: "clustering.smk"
include: "mapping.smk"

rule mmseqs2_convertalis_blast:
     input:
         query=rules.create_mmseqs2_query_db.output.query_path,
         target=rules.create_mmseqs2_target_db.output.uniprot90_path,
         map=rules.mmseqs2_map.output.out_dir
     output: "results/{database}/mmseqs2/convertalis/{database}_convertlis.8"
     params:
         prefix=rules.mmseqs2_map.params.prefix,
         query_prefix=rules.create_mmseqs2_query_db.params.query_prefix,
         target_prefix=rules.create_mmseqs2_target_db.params.uniprot90_prefix
     threads: config["mmseqs2"]["convertalis"]["threads"]
     conda: "../envs/database_management.yml"
     shell:
         """
         mmseqs convertalis {input.query}/{params.query_prefix} \
             {input.target}/{params.target_prefix} \
             {input.map}/{params.prefix} {output} --format-mode 4 \
             --format-output query,target,pident
         """

rule mmseqs2_convertalis_unmapped_blast:
     input:
         query=rules.create_mmseqs2_unaligned_db.output.out_dir,
         target=rules.create_mmseqs2_target_db.output.uniprot90_path,
         cluster=rules.mmseqs2_linclust.output.out_dir
     output: "results/{database}/mmseqs2/convertalis/unaligned/{database}_convertlis.8"
     params:
         prefix=rules.mmseqs2_linclust.params.prefix,
         query_prefix=rules.create_mmseqs2_unaligned_db.params.unaligned_prefix,
         target_prefix=rules.create_mmseqs2_target_db.params.uniprot90_prefix
     threads: config["mmseqs2"]["convertalis"]["threads"]
     conda: "../envs/database_management.yml"
     shell:
         """
         mmseqs createtsv {input.query}/{params.query_prefix} \
             {input.target}/{params.target_prefix} \
             {input.cluster}/{params.prefix} {output}
         """
