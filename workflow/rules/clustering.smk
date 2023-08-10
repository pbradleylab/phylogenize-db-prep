include: "mapping.smk"


rule get_unaligned_uniref50_sequences:
    input: 
        aligned=rules.get_aligned_uniref50_contigs.output,
        all_sequences=rules.combine_fasta_uniref50.output
    output: "results/{database}/uniref50/mapping/unmapped/{database}.fa"
    conda: "../envs/clustering.yml"
    shell:
        """
        grep '>' {input.all_sequences} | sed "s/>//g" > /tmp/tmp.all
        grep -F -v -x -f {input.aligned} /tmp/tmp.all > /tmp/tmp.unaligned
        faSomeRecords {input.all_sequences} /tmp/tmp.unaligned {output}
        """

# Create a new database that is declared as temporary. This database
# holds the unaligned peptide sequences that are assumed as potential
# species specific alignments.
rule create_mmseqs2_unaligned_uniref50_db:
    input: rules.get_unaligned_uniref50_sequences.output
    output:
        outdir=directory("results/{database}/uniref50/mmseqs2/unmapped/"),
        index="results/{database}/uniref50/mmseqs2/unmapped/unmapped.index"
    params:
        unaligned_prefix="unmapped"
    conda: "../envs/clustering.yml"
    log: "logs/{database}/uniref50/mmseqs2/create_mmseqs2_unaligned/mmseqs2_create_mmseqs2_unaligned.log"
    shell:
        """
        mkdir -p {output.outdir}
        mmseqs createdb {input} {output.outdir}/{params.unaligned_prefix} \
            --dbtype 1 2> {log}
        mmseqs createindex {output.outdir}/{params.unaligned_prefix} \
            /tmp 2> {log}
        """

rule mmseqs2_map_uhgp50:
     input:
         query=rules.create_mmseqs2_unaligned_uniref50_db.output.outdir,
         target=rules.create_uhgp50.params.uhgp50_path
     output:
         outdir=directory("results/{database}/mmseqs2/uhgp50/mapping/"),
         index="results/{database}/mmseqs2/uhgp50/mapping/{database}_map.index"
     log: "logs/{database}/uhgp50/mmseqs2/mapping/mmseqs2_map.log"
     params:
         prefix="{database}_map",
         query_prefix=rules.create_mmseqs2_unaligned_uniref50_db.params.unaligned_prefix,
         target_prefix="uhgp50"
     threads: config["mmseqs2"]["map"]["threads"]
     conda: "../envs/mapping.yml"
     shell:
         """
         mmseqs map -a --threads {threads} {input.query}/{params.query_prefix} \
            {input.target}/{params.target_prefix} {output.outdir}/{params.prefix} \
            --min-seq-id 0.50 /tmp 2> {log}
         """

rule get_aligned_uhgp50_contigs:
     input:
         query=rules.create_mmseqs2_query_db.output.query_path,
         target=rules.create_uhgp50.output,
         mapped=rules.mmseqs2_map_uhgp50.output.outdir
     output: "results/{database}/uhgp50/mmseqs2/convertalis/{database}_aligned.8"
     params:
         prefix=rules.mmseqs2_map_uhgp50.params.prefix,
         query_prefix=rules.create_mmseqs2_query_db.params.query_prefix,
         target_prefix="uhgp50"
     threads: config["mmseqs2"]["convertalis"]["threads"]
     conda: "../envs/mapping.yml"
     shell:
         """
         mmseqs convertalis \
             {input.query}/{params.query_prefix} \
             {input.target}/{params.target_prefix} \
             {input.mapped}/{params.prefix} \
             {output} --format-mode 4 \
             --format-output query && sed -i '1d' {output}
         """

rule get_unaligned_uhgp50_sequences:
    input: 
        aligned_uhgp50=rules.get_aligned_uhgp50_contigs.output,
        aligned_uniref50=rules.get_aligned_uniref50_contigs.output,
        all_sequences=rules.combine_fasta_uniref50.output
    output: "results/{database}/uhgp50/mapping/unmapped/{database}.fa"
    conda: "../envs/clustering.yml"
    shell:
        """
        cat {input.aligned_uniref50} {input.aligned_uhgp50} > /tmp/tmp.aligned
        grep '>' {input.all_sequences} | sed "s/>//g" > /tmp/tmp.all
        grep -F -v -x -f /tmp/tmp.aligned /tmp/tmp.all > /tmp/tmp.unaligned
        faSomeRecords {input.all_sequences} /tmp/tmp.unaligned {output}
        """

# Create a new database that is declared as temporary. This database
# holds the unaligned peptide sequences that are assumed as potential
# species specific alignments.
rule create_mmseqs2_unaligned_uhgp50_db:
    input: rules.get_unaligned_uhgp50_sequences.output
    output:
        outdir=directory("results/{database}/uhgp50/mmseqs2/unmapped/"),
        index="results/{database}/uhgp50/mmseqs2/unmapped/unmapped.index"
    params:
        unaligned_prefix="unmapped"
    conda: "../envs/clustering.yml"
    log: "logs/{database}/uhgp50/mmseqs2/create_mmseqs2_unaligned/mmseqs2_create_mmseqs2_unaligned.log"
    shell:
        """
        mkdir -p {output.outdir}
        mmseqs createdb {input} {output.outdir}/{params.unaligned_prefix} \
            --dbtype 1 2> {log}
        mmseqs createindex {output.outdir}/{params.unaligned_prefix} \
            /tmp 2> {log}
        """

# Cluster the unaligned protein sequence database from
# mmseqs' search command.
rule mmseqs2_linclust_uhgp50_db:
     input: rules.create_mmseqs2_unaligned_uhgp50_db.output.outdir
     output:
         database="results/{database}/uhgp50/mmseqs2/linclust/unaligned_linclust.index",
         outdir=directory("results/{database}/uhgp50/mmseqs2/linclust/")
     params:
         unaligned_prefix=rules.create_mmseqs2_unaligned_uhgp50_db.params.unaligned_prefix,
         prefix="unaligned_linclust",
         seq_id_precent=config["mmseqs2"]["linclust"]["seq_id_precent"],
         tmp_dir=config["mmseqs2"]["linclust"]["tmp_dir"]
     conda: "../envs/clustering.yml"
     log: "logs/{database}/uhgp50/mmseqs2/linclust/mmseqs2_linclust.log"
     threads: config["mmseqs2"]["linclust"]["threads"]
     shell:
         """
         mmseqs linclust {input}/{params.unaligned_prefix} {output.outdir}/{params.prefix} \
            {params.tmp_dir} --min-seq-id {params.seq_id_precent} \
            --threads {threads} 2> {log}
         """
