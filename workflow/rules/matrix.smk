include: "clustering.smk"


def get_16s(wildcards):
    return(config["files"]["16S"][wildcards.mapping_db]["fasta"])


# Combines species with a 90% or greater identity match to the target database, 
# and the unmapped regions to a list of species specific vectors by their centroid.
rule combine_species_hits:
    input:
        tophits=lambda wildcards: [
            f for f in rules.get_top_50_evals.output.tophits
            if f"_{wildcards.mapping_db}_" in f and f.startswith(f"results/{wildcards.database}/")],
        clustered=lambda wildcards: f"results/{wildcards.database}/mmseqs2_linclust/{wildcards.mapping_db}/unaligned_linclust_cluster.tsv",
        checkpoints=lambda wildcards: [
            f"results/{wildcards.database}/checkpoints/database_processing_checkpoint/{wildcards.mapping_db}_{target_db}_processed.done"
            for target_db in TARGET_DBS]
    output:"results/{database}/binary/combined_species_hits/{mapping_db}_{database}.tsv"
    shell:
        """
        echo -e "query\ttarget" > {output}
        cat {input.clustered} | cut -f1  > /tmp/tmp.1
        cat {input.clustered} | cut -f2  > /tmp/tmp.2
        paste /tmp/tmp.2 /tmp/tmp.1 >> {output} 
        # Combine all tophits files with appropriate headers
        for file in {input.tophits}; do
            db_info=$(basename $file | sed 's/_convertlis_tophits.tsv//')
            cut -f1,2 $file | awk -v db="`$db_info`" '{{print $0"\t"db}}' >> {output}
        done
        """

# Create the taxonomy file and generate the input for making the binary file needed
# to pass into R to compress and spearate into a list of binaries per phylum.
rule get_taxonomy:
    input: rules.combine_species_hits.output
    output: 
        out=temp("results/{database}/binary/get_taxonomy/{mapping_db}-binary-temp.csv"),
        tax="results/{database}/binary/get_taxonomy/{mapping_db}-taxonomy.csv"
    params:
        split_char=config["create_species_matrix"]["split_char"],
        tax=lambda wildcards: config["files"]["taxonomy"][wildcards.mapping_db],
        mapping=lambda wildcards: config["files"]["mapping"][wildcards.mapping_db]
    conda: "../envs/matrix.yml"
    shell:
        """
        python workflow/scripts/reformat_taxonomy.py \
            --input {input} \
            --tax {params.tax} \
            --split_char {params.split_char} \
            --output {output.out} \
            --tax_output {output.tax} \
            --mapping {params.mapping}
        """

rule get_binary:
    input: rules.get_taxonomy.output.out
    output: "results/{database}/binary/get_binary/{mapping_db}-binary.rds"
    conda: "../envs/matrix.yml"
    shell:
        """
        Rscript workflow/scripts/make_binary.R {input} {output} 
        """

rule get_tree:
    input: rules.get_taxonomy.output.out
    output: "results/{database}/binary/get_tree/{mapping_db}-tree.rds"
    params:
        tree=lambda wildcards: config["files"]["tree"][wildcards.mapping_db]
    conda: "../envs/matrix.yml"
    shell:
        """
        Rscript workflow/scripts/make_tree.R {input} {params.tree} {output}
        """

rule get_16s:
    input:get_16s
    output: "results/{database}/binary/get_16s/{mapping_db}.faa"
    conda: "../envs/matrix.yml"
    shell:
        """
        echo "to be completed"  
        """
