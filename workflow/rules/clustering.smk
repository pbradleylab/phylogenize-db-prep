# Cluster the unaligned protein sequence database from
# mmseqs' search command.
#include: "resources.smk"
rule mmseqs2_linclust:
     input: rules.create_mmseqs2_unaligned_db.output.out_dir
     output:
         database="results/{database}/mmseqs2/linclust/unaligned_linclust.index",
         out_dir=directory("results/{database}/mmseqs2/linclust/")
     params:
         unaligned_prefix=rules.create_mmseqs2_unaligned_db.params.unaligned_prefix,
         prefix="unaligned_linclust",
         seq_id_precent=config["mmseqs2"]["linclust"]["seq_id_precent"],
         tmp_dir=config["mmseqs2"]["linclust"]["tmp_dir"]
     conda: "../envs/transformation.yml"
     log: "logs/{database}/mmseqs2/linclust/mmseqs2_linclust.log"
     threads: config["mmseqs2"]["linclust"]["threads"]
     shell:
         """
         mmseqs linclust {input}/{params.unaligned_prefix} {output.out_dir}/{params.prefix} \
            {params.tmp_dir} --min-seq-id {params.seq_id_precent} \
            --threads {threads} 2> {log}
         """
