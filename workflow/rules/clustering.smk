include: "mapping.smk"


rule get_unaligned_sequences:
    input: 
        aligned=rules.samtools_fasta.output,
        all_sequences=rules.combine_fasta.output
    output: "results/{database}/mapping/unmapped/{database}.fa"
    conda: "../envs/clustering.yml"
    shell:
        """
        grep '>' {input.aligned} | sed "s/>//g" > /tmp/tmp.aligned
        grep '>' {input.all_sequences} | sed "s/>//g" > /tmp/tmp.all
        grep -F -v -x -f /tmp/tmp.aligned /tmp/tmp.all > /tmp/tmp.unaligned
        faSomeRecords {input.all_sequences} /tmp/tmp.unaligned {output}
        """

# Create a new database that is declared as temporary. This database
# holds the unaligned peptide sequences that are assumed as potential
# species specific alignments.
rule create_mmseqs2_unaligned_db:
    input: rules.get_unaligned_sequences.output
    output:
        out_dir=directory("resources/{database}/mmseqs2/unmapped/"),
        index="resources/{database}/mmseqs2/unmapped/unmapped.index"
    params:
        unaligned_prefix="unmapped"
    conda: "../envs/database_management.yml"
    log: "logs/{database}/mmseqs2/create_mmseqs2_unaligned/mmseqs2_create_mmseqs2_unaligned.log"
    shell:
        """
        mkdir -p {output.out_dir}
        mmseqs createdb {input} {output.out_dir}/{params.unaligned_prefix} \
            --dbtype 1 2> {log}
        mmseqs createindex {output.out_dir}/{params.unaligned_prefix} \
            /tmp 2> {log}
        """

# Cluster the unaligned protein sequence database from
# mmseqs' search command.
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
