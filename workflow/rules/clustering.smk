include: "mapping.smk"

TARGET_DBS = list(config["target_db"]["paths"].keys())

# Convert mmseqs results to blast format
rule mmseqs2_convertalis:
     input:
         query=rules.make_current_query_databases.output.query_path,
         target=rules.make_targets.output.target_path,
         linker=rules.map_query.output
     params:
         map=rules.map_query.params.outdir
     output:
         blast="results/{database}/clustering/mmseqs2/convertalis/{target_db}_{mapping_db}_convertlis.8",
         list="results/{database}/clustering/mmseqs2/convertalis/{target_db}_{mapping_db}_convertlis.list"
     conda: "../envs/clustering.yml"
     log: "logs/{database}/clustering/mmseqs2_convertalis/{target_db}_{mapping_db}.log"
     threads: config["mmseqs2"]["convertalis"]["threads"]
     shell:
         """
         mmseqs convertalis {input.query}/{wildcards.mapping_db} \
             {input.target}/{wildcards.target_db} \
             {params.map}/{wildcards.mapping_db} {output.blast} --format-mode 4 \
             --format-output query,target,pident 2> {log}
         cut -f1 {output.blast} | sed '1d' > {output.list}
         """

rule mmseqs2_convertalis_full_blast:
     input:
         query=rules.make_current_query_databases.output.query_path,
         target=rules.make_targets.output.target_path,
         linker=rules.map_query.output
     params:
         map=rules.map_query.params.outdir
     output:"results/{database}/clustering/mmseqs2/convertalis_full_blast/{target_db}_{mapping_db}_convertlis.8"
     conda: "../envs/clustering.yml"
     log: "logs/{database}/clustering/mmseqs2_convertalis/{target_db}_{mapping_db}.log"
     threads: config["mmseqs2"]["convertalis"]["threads"]
     shell:
         """
         mmseqs convertalis {input.query}/{wildcards.mapping_db} \
             {input.target}/{wildcards.target_db} \
             {params.map}/{wildcards.mapping_db} {output} --format-mode 4 \
             --format-output "query,qlen,target,tlen,evalue,bitscore,alnlen,nident" 2> {log}
             #--format-output query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits 2> {log}
         sed -i '1d' {output}
         """

# Filter results to get top hits
rule get_top_50_evals:
     input: 
         corrected=rules.mmseqs2_convertalis.output.blast,
         full=rules.mmseqs2_convertalis_full_blast.output
     output:
         unfiltered="results/{database}/clustering/mmseqs2/top_50/{target_db}_{mapping_db}_convertlis.tsv",
         tophits="results/{database}/clustering/mmseqs2/top_50/{target_db}_{mapping_db}_convertlis_tophits.tsv"
     params:
         outdir="results/{database}/clustering/mmseqs2/top_50/"
     shell:
         """
         touch {output.unfiltered}
         awk '$3>50 {{print}}' {input.corrected} > {output.unfiltered}
         python workflow/scripts/get_top_hits.py -i {output.unfiltered} -o {output.tophits}
         """

# Extract sequences that didn't map for next iteration
rule faSomeRecords:
    input:
        aligned=rules.get_top_50_evals.output.unfiltered,
        all_sequences=rules.prepare_current_iteration_input.output
    output:
        unmapped="results/{database}/clustering/faSomeRecords/unmapped/{target_db}_{mapping_db}.fa",
        cumulative_unmapped="results/{database}/clustering/faSomeRecords/cumulative_unmapped/{mapping_db}_after_{target_db}.fa",
        tmp_unaligned="results/{database}/clustering/faSomeRecords/cumulative_unmapped/{mapping_db}_after_{target_db}.unaligned",
        tmp_aligned="results/{database}/clustering/faSomeRecords/cumulative_unmapped/{mapping_db}_after_{target_db}.aligned",
        tmp_all="results/{database}/clustering/faSomeRecords/cumulative_unmapped/{mapping_db}_after_{target_db}.all"
    conda: "../envs/clustering.yml"
    shell:
        """
        cut -f1 {input.aligned} | sed '1d' | uniq > {output.tmp_aligned}
        grep '>' {input.all_sequences} | sed "s/>//g" | sed "s/ .*//g" > {output.tmp_all}
        grep -F -v -x -f {output.tmp_aligned} {output.tmp_all} > {output.tmp_unaligned}
        # Current iteration's unmapped
        faSomeRecords {input.all_sequences} {output.tmp_unaligned} {output.unmapped}
        # Cumulative unmapped for next iterations
        cp {output.unmapped} {output.cumulative_unmapped}
        """

# Combine all final unmapped sequences from the last target database
rule combine_final_unmapped_sequences:
    input:"results/{database}/clustering/faSomeRecords/cumulative_unmapped/{mapping_db}_after_" + TARGET_DBS[-1] + ".fa",
    output:"results/{database}/clustering/faSomeRecords/final_unmapped/{mapping_db}_clustered.fa"
    shell:
        """
        cat {input} > {output}
        """

# Run linclust on the final unmapped sequences
rule mmseqs2_linclust:
    input: lambda wildcards: rules.combine_final_unmapped_sequences.output
    output:
        outdir=directory("results/{database}/clustering/mmseqs2_linclust/{mapping_db}"),
        tsv="results/{database}/clustering/mmseqs2_linclust/{mapping_db}/unaligned_linclust_cluster.tsv",
        fasta="results/{database}/clustering/mmseqs2_linclust/{mapping_db}/unaligned_linclust_rep_seq.fasta"
    params:
        prefix="unaligned_linclust",
        tmp_dir=config["mmseqs2"]["linclust"]["tmp_dir"]
    conda: "../envs/clustering.yml"
    log: "logs/{database}/clustering/mmseqs2_linclust/{mapping_db}/mmseqs2_linclust.log"
    threads: config["mmseqs2"]["linclust"]["threads"]
    shell:
        """
        mmseqs easy-linclust {input} {output.outdir}/{params.prefix} \
            {params.tmp_dir} --threads {threads} 2> {log}
        """
