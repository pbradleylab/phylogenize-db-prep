include: "mapping.smk"

rule mmseqs2_convertalis_blast_uniref50_db:
     input:
         query=rules.create_mmseqs2_query_db.output.query_path,
         target=rules.create_uniref50.output.uniref50_path,
         map=rules.mmseqs2_map_uniref50.output.outdir
     output: 
         blast="results/{database}/uniref50/mmseqs2/convertalis/{database}_convertlis.8",
         list="results/{database}/uniref50/mmseqs2/convertalis/{database}_convertlis.list"
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
             {input.map}/{params.prefix} {output.blast} --format-mode 4 \
             --format-output query,target,pident
         cut -f1 {output.blast} | sed '1d' > {output.list}
         """

rule get_top_50_evals_uniref50:
     input: rules.mmseqs2_convertalis_blast_uniref50_db.output.blast
     output: 
         unfiltered="results/{database}/uniref50/mmseqs2/top_50/{database}_convertlis.tsv",
         tophits="results/{database}/uniref50/mmseqs2/top_50/{database}_convertlis_tophits.tsv"
     conda: "../envs/matrix.yml"
     shell:
         """
         touch {output.unfiltered}
         awk '$3>50 {{print}}' {input} > {output.unfiltered}
         python workflow/scripts/get_top_hits.py -i {output.unfiltered} -o {output.tophits}
         """

rule get_unaligned_uniref50_sequences:
    input: 
        aligned=rules.get_top_50_evals_uniref50.output.unfiltered,
        all_sequences=rules.combine_fasta_uniref50.output
    output: "results/{database}/uniref50/mapping/unmapped/{database}.fa"
    conda: "../envs/clustering.yml"
    shell:
        """
        cut -f1 {input.aligned} | sed '1d' | uniq > /tmp/tmp.aligned
        grep '>' {input.all_sequences} | sed "s/>//g" > /tmp/tmp.all
        grep -F -v -x -f /tmp/tmp.aligned /tmp/tmp.all > /tmp/tmp.unaligned
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
         outdir=directory("results/{database}/uhgp50/mmseqs2/mapping/"),
         index="results/{database}/uhgp50/mmseqs2/mapping/{database}_map.index"
     log: "logs/{database}/uhgp50/mmseqs2/mapping/mmseqs2_map.log"
     params:
         prefix="{database}_map",
         query_prefix=rules.create_mmseqs2_unaligned_uniref50_db.params.unaligned_prefix,
         target_prefix="uhgp50"
     threads: config["mmseqs2"]["map"]["threads"]
     conda: "../envs/mapping.yml"
     shell:
         """
         mmseqs map --threads {threads} {input.query}/{params.query_prefix} \
            {input.target}/{params.target_prefix} {output.outdir}/{params.prefix} \
            /tmp -a --comp-bias-corr 0 --mask 0 --min-seq-id 0.50 2> {log}
         """

rule mmseqs2_convertalis_blast_uhgp50_db:
     input:
         query=rules.create_mmseqs2_unaligned_uniref50_db.output.outdir,
         target=rules.create_uhgp50.output,
         map=rules.mmseqs2_map_uhgp50.output.outdir
     output: 
         blast="results/{database}/uhgp50/mmseqs2/convertalis/{database}_convertlis.8",
         list="results/{database}/uhgp50/mmseqs2/convertalis/{database}_convertlis.list"
     params:
         prefix=rules.mmseqs2_map_uhgp50.params.prefix,
         query_prefix=rules.create_mmseqs2_unaligned_uniref50_db.params.unaligned_prefix,
         target_prefix=rules.create_uhgp50.output
     threads: config["mmseqs2"]["convertalis"]["threads"]
     conda: "../envs/blast.yml"
     shell:
         """
         mmseqs convertalis {input.query}/{params.query_prefix} \
             {input.target}/uhgp50 \
             {input.map}/{params.prefix} {output.blast} --format-mode 4 \
             --format-output query,target,pident
         cut -f1 {output.blast} | sed '1d' > {output.list}
         """

rule get_top_50_evals_uhgp50:
     input: rules.mmseqs2_convertalis_blast_uhgp50_db.output.blast
     output: 
         unfiltered="results/{database}/uhgp50/mmseqs2/top_50/{database}_convertlis.tsv",
         tophits="results/{database}/uhgp50/mmseqs2/top_50/{database}_convertlis_tophits.tsv"
     conda: "../envs/matrix.yml"
     shell:
         """
         touch {output.unfiltered}
         awk '$3>50 {{print}}' {input} > {output.unfiltered}
         python3 workflow/scripts/get_top_hits.py -i {output.unfiltered} -o {output.tophits}
         """

rule get_unaligned_uhgp50_sequences:
    input: 
        aligned_uhgp50=rules.get_top_50_evals_uhgp50.output.unfiltered,
        aligned_uniref50=rules.get_top_50_evals_uniref50.output.unfiltered,
        all_sequences=rules.get_unaligned_uniref50_sequences.output
    output: "results/{database}/uhgp50/mapping/unmapped/{database}.fa"
    conda: "../envs/clustering.yml"
    shell:
        """
        cut -f1 {input.aligned_uhgp50} | sed '1d' | uniq > /tmp/tmp.aligned
        grep '>' {input.all_sequences} | sed "s/>//g" | uniq > /tmp/tmp.all
        grep -F -v -x -f /tmp/tmp.aligned /tmp/tmp.all | uniq > /tmp/tmp.unaligned
        faSomeRecords {input.all_sequences} /tmp/tmp.unaligned {output}
        """

# Cluster the unaligned protein sequence database from
# mmseqs' search command.
rule mmseqs2_linclust_uhgp50_db:
     input:rules.get_unaligned_uhgp50_sequences.output
     output:
         outdir=directory("results/{database}/uhgp50/mmseqs2/linclust/"),
         tsv="results/{database}/uhgp50/mmseqs2/linclust/unaligned_linclust_cluster.tsv"
     params:
         prefix="unaligned_linclust",
         tmp_dir=config["mmseqs2"]["linclust"]["tmp_dir"]
     conda: "../envs/clustering.yml"
     log: "logs/{database}/uhgp50/mmseqs2/linclust/mmseqs2_linclust.log"
     threads: config["mmseqs2"]["linclust"]["threads"]
     shell:
         """
         mmseqs easy-linclust {input[0]} {output.outdir}/{params.prefix} \
            {params.tmp_dir} --threads {threads} 2> {log}
         """

