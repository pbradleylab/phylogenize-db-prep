include: "mapping.smk"

TARGET_DBS = list(config["target_db"]["paths"].keys())

# Convert mmseqs results to blast format
rule mmseqs2_convertalis:
     input:
         query=rules.make_current_query_databases.output.query_path,
         target=rules.make_targets.output.target_path,
         map=rules.map_query.output.outdir
     output:
         blast="results/{database}/annotation/mmseqs2/convertalis/{target_db}_{mapping_db}_convertlis.8",
         list="results/{database}/annotation/mmseqs2/convertalis/{target_db}_{mapping_db}_convertlis.list"
     conda: "envs/blast.yml"
     log: "logs/{database}/annotation/mmseqs2_convertalis_blast/{target_db}_{mapping_db}.log"
     threads: config["mmseqs2"]["convertalis"]["threads"]
     shell:
         """
         mmseqs convertalis {input.query}/{wildcards.mapping_db} \
             {input.target}/{wildcards.target_db} \
             {input.map}/{wildcards.mapping_db} {output.blast} --format-mode 4 \
             --format-output query,target,pident 2> {log}
         cut -f1 {output.blast} | sed '1d' > {output.list}
         """

# Filter results to get top hits
rule get_top_50_evals:
     input: rules.mmseqs2_convertalis.output.blast
     output:
         unfiltered="results/{database}/annotation/mmseqs2/top_50/{target_db}_{mapping_db}_convertlis.tsv",
         tophits="results/{database}/annotation/mmseqs2/top_50/{target_db}_{mapping_db}_convertlis_tophits.tsv"
     params:
         outdir="results/{database}/annotation/mmseqs2/top_50/"
     shell:
         """
         touch {output.unfiltered}
         awk '$3>50 {{print}}' {input} > {output.unfiltered}
         python workflow/scripts/get_top_hits.py -i {output.unfiltered} -o {output.tophits}
         """

# Extract sequences that didn't map for next iteration
rule get_unaligned_sequences:
    input:
        aligned=rules.get_top_50_evals.output.unfiltered,
        all_sequences=rules.prepare_current_iteration_input.output.current_input
    output:
        unmapped="results/{database}/annotation/faSomeRecords/unmapped/{target_db}_{mapping_db}.fa",
        cumulative_unmapped="results/{database}/annotation/faSomeRecords/cumulative_unmapped/{mapping_db}_after_{target_db}.fa"
    conda: "envs/clustering.yml"
    shell:
        """
        cut -f1 {input.aligned} | sed '1d' | uniq > /tmp/tmp.aligned
        grep '>' {input.all_sequences} | sed "s/>//g" > /tmp/tmp.all
        grep -F -v -x -f /tmp/tmp.aligned /tmp/tmp.all > /tmp/tmp.unaligned
        # Current iteration's unmapped
        faSomeRecords {input.all_sequences} /tmp/tmp.unaligned {output.unmapped}
        # Cumulative unmapped for next iterations
        cp {output.unmapped} {output.cumulative_unmapped}
        """

# Combine all final unmapped sequences from the last target database
rule combine_final_unmapped_sequences:
    input:
        final_unmapped=expand("results/{{database}}/annotation/faSomeRecords/cumulative_unmapped/{mapping_db}_after_{last_target_db}.fa",
                              mapping_db=config["annotation"]["mapping_databases"].keys(),
                              last_target_db=TARGET_DBS[-1]),
        # Ensure all databases have finished processing
        checkpoints=expand("results/{{database}}/annotation/checkpoints/{target_db}_processed.done",
                          target_db=TARGET_DBS)
    output:
        combined="results/{database}/annotation/faSomeRecords/final_unmapped/all_final_unmapped.fa"
    shell:
        """
        cat {input.final_unmapped} > {output.combined}
        """

# Run linclust on the final unmapped sequences
rule mmseqs2_linclust:
    input: rules.combine_final_unmapped_sequences.output.combined
    output:
        outdir=directory("results/{database}/mmseqs2/linclust/"),
        tsv="results/{database}/mmseqs2/linclust/unaligned_linclust_cluster.tsv"
    params:
        prefix="unaligned_linclust",
        tmp_dir=config["mmseqs2"]["linclust"]["tmp_dir"]
    conda: "envs/clustering.yml"
    log: "logs/{database}/mmseqs2/linclust/mmseqs2_linclust.log"
    threads: config["mmseqs2"]["linclust"]["threads"]
    shell:
        """
        mmseqs easy-linclust {input} {output.outdir}/{params.prefix} \
            {params.tmp_dir} --threads {threads} 2> {log}
        """
