include: "clustering.smk"


# Rule to collect all results at the end - depends on all databases being processed
rule collect_annotation_results:
    input:
        tophits=expand("results/{{database}}/annotation/mmseqs2/top_50/{target_db}_{mapping_db}_convertlis_tophits.tsv",
                       target_db=TARGET_DBS,
                       mapping_db=config["annotation"]["mapping_databases"].keys()),
        # Wait for all database checkpoints to complete
        checkpoints=expand("results/{{database}}/annotation/checkpoints/{target_db}_processed.done",
                          target_db=TARGET_DBS)
    output:
        combined="results/{database}/annotation/combined_results.tsv"
    shell:
        """
        # Combine all tophits files with appropriate headers
        echo -e "query\ttarget\tpident\tdatabase" > {output.combined}
        for file in {input.tophits}; do
            db_info=$(basename $file | sed 's/_convertlis_tophits.tsv//')
            awk -v db="$db_info" '{{print $0"\t"db}}' $file >> {output.combined}
        done
        """

# Generate final report
rule generate_annotation_report:
    input:
        combined=rules.collect_annotation_results.output.combined,
        unmapped=expand("results/{{database}}/annotation/faSomeRecords/unmapped/{target_db}_{mapping_db}.fa",
                        target_db=TARGET_DBS,
                        mapping_db=config["annotation"]["mapping_databases"].keys()),
        linclust=rules.mmseqs2_linclust.output.tsv
    output:
        report="results/{database}/annotation/final_report.html"
    conda: "envs/clustering.yml"
    script:
        "../scripts/generate_report.py"

# Combines species with a 90% or greater identity match to the target database, 
# and the unmapped regions to a list of species specific vectors by their centroid.
rule combine_species_hits:
    input:
        uhgp50_identity_50=rules.get_top_50_evals_uhgp50.output.tophits,
        unmapped=rules.mmseqs2_linclust_uhgp50_db.output.tsv,
        uniref50_identity_50=rules.get_top_50_evals_uniref50.output.tophits
    output: 
         txt="results/{database}/uniref50/mmseqs2/combined_species_hits/{database}.txt",
         outdir=directory("results/{database}/uniref50/mmseqs2/combined_species_hits/")
    conda: "../envs/matrix.yml"
    log: "logs/{database}/uniref50/mmseqs2/hits_50/mmseqs2_hits_50.log"
    shell:
        """
        cat {input.unmapped} | cut -f1  > /tmp/tmp.1
        cat {input.unmapped} | cut -f2  > /tmp/tmp.2 
        paste /tmp/tmp.2 /tmp/tmp.1 > {output.txt}
        cat {input.uhgp50_identity_50} | cut -f1,2 | sed '1d' >> {output.txt}
        cat {input.uniref50_identity_50} | cut -f1,2 | sed '1d' >> {output.txt}
        """

# Combine any database files run in previous iterations together for the final matrix
rule combine_hits:
    input: rules.combine_species_hits.output
    output: 
        txt="results/{database}/combined_hits/{database}.txt",
        outdir=directory("results/{database}/combined_hits/")
    params:
        search_dir="results/{database}/",
        name="{database}.txt"
    shell:
        """
        bash workflow/scripts/combine_hits.sh {params.search_dir} {params.name} {output.txt}
        """


# Create a binary matrix of the hits > 50% and those that unmapped 
# but clustered.
rule create_species_matrix:
    input: rules.combine_hits.output.outdir
    output: 
        out="results/{database}/final/species_matrix/{database}.txt",
        dup="results/{database}/final/species_matrix/{database}.dups"
    conda: "../envs/matrix.yml"
    shell:
        """
        #mkdir -p $(dirname {output.out})
        #touch {output.dup} && chmod 777 {output.dup}
        python workflow/scripts/combine_species.py \
            --output {output.out} \
            --dir {input} \
            --ext ".txt" \
            --duplicates {output.dup}
        """
