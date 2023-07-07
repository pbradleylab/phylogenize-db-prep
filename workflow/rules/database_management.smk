include: "utilities.smk"

# Converts the database's mappings to a sam format. The unmapped (unaligned)
# sequences are then taken to generate a new database in rule:
# `samtools_get_aligned`.
rule mmseqs2_convertalis_sam:
     input:
         query=rules.create_mmseqs2_query_db.output.query_path,
         target=rules.create_mmseqs2_target_db.output.uniprot90_path,
         mapped=rules.mmseqs2_map.output.out_dir
     output: "results/{database}/mmseqs2/convertalis/{database}_convertlis.sam"
     params:
         prefix=rules.mmseqs2_map.params.prefix,
         query_prefix=rules.create_mmseqs2_query_db.params.query_prefix,
         target_prefix=rules.create_mmseqs2_target_db.params.uniprot90_prefix
     threads: config["mmseqs2"]["convertalis"]["threads"]
     conda: "../envs/database_management.yml"
     shell:
         """
         mmseqs convertalis \
             {input.query}/{params.query_prefix} \
             {input.target}/{params.target_prefix} \
             {input.mapped}/{params.prefix} \
             {output} --format-mode 1
         """

# Create a new database that is declared as temporary. This database
# holds the unaligned peptide sequences that are assumed as potential
# species specific alignments.
rule create_mmseqs2_unaligned_db:
    input: rules.get_unaligned_sequences.output
    output:
        out_dir=directory("resources/{database}/mmseqs2/unmapped/"),
        index="resources/{database}/mmseqs2/unmapped/unmapped.index",
    params:
        unaligned_prefix="unmapped"
    conda: "../envs/database_management.yml"
    log: "logs/{database}/mmseqs2/create_mmseqs2_unaligned/mmseqs2_create_mmseqs2_unaligned.log"
    shell:
        """
        mkdir -p {output.out_dir}
        mmseqs createdb {input} {output.out_dir}/{params.unaligned_prefix} --dbtype 1 2> {log}
        mmseqs createindex {output.out_dir}/{params.unaligned_prefix} \
            /tmp 2> {log}
        """

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
         map=rules.mmseqs2_map.output.out_dir
     output: "results/{database}/mmseqs2/convertalis/unaligned/{database}_convertlis.8"
     params:
         prefix=rules.mmseqs2_linclust.params.prefix,
         query_prefix=rules.create_mmseqs2_unaligned_db.params.query_prefix,
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
