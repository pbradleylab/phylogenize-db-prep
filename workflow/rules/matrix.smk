include: "clustering.smk"

import shlex



def get_top_hits(wildcards):
    out=expand(rules.get_top_50_evals.output.tophits, mapping_db=wildcards.mapping_db, target_db=config["target_db"]["paths"].keys(), database=wildcards.database)
    return(out)

def get_all_hits(wildcards):
    out=expand(rules.mmseqs2_convertalis_full_blast.output, mapping_db=wildcards.mapping_db, target_db=config["target_db"]["paths"].keys(), database=wildcards.database)
    return(out)

# Combines species with a 90% or greater identity match to the target database, 
# and the unmapped regions to a list of species specific vectors by their centroid.
rule combine_species_hits:
    input:
        clustered=rules.mmseqs2_linclust.output.tsv,
        tophits=get_top_hits,
        allhits=get_all_hits
    output:
        tophits="results/{database}/binary/combined_species_hits/{mapping_db}_{database}.tsv",
        allhits="results/{database}/binary/combined_species_hits/{mapping_db}_{database}.all"
    shell:
        """
        echo -e "query\ttarget" > {output.tophits}

        for file in {input.tophits}; do
            db_info=$(basename $file | sed 's/_convertlis_tophits.tsv//')
            cut -f1,2 $file | awk -v db="$db_info" '{{print $0"\t"db}}' >> /tmp/{wildcards.mapping_db}.1
        done
        
        # MMseqs cluster tables contain representative,member; queries are members.
        awk -v db="denovo" '{{print $2"\t"$1"\t"db}}' {input.clustered} >> /tmp/{wildcards.mapping_db}.1

        # Ensure only one gene hit is used per annotation and in the order of databases selected.
        awk '!seen[$1]++' /tmp/{wildcards.mapping_db}.1 | cut -f1,2 >> {output.tophits}
        rm /tmp/{wildcards.mapping_db}.1

        cat {input.allhits} > {output.allhits}
        """

# Create the taxonomy file and generate the input for making the binary file needed
# to pass into R to compress and separate into a list of binaries per phylum.
rule get_taxonomy:
    input:
        hits=rules.combine_species_hits.output.tophits,
        mapping=lambda wildcards: config["files"]["mapping"][wildcards.mapping_db],
        taxonomy=lambda wildcards: config["files"]["taxonomy"][wildcards.mapping_db],
        script="workflow/scripts/reformat_taxonomy.py"
    output: 
        out=temporary("results/{database}/binary/get_taxonomy/{mapping_db}-binary-temp.csv"),
        tax="results/{database}/binary/get_taxonomy/{mapping_db}-taxonomy.csv"
    params:
        split_char=config["create_species_matrix"]["split_char"]
    conda: "../envs/matrix.yml"
    shell:
        """
        python {input.script:q} \
            --input {input.hits:q} \
            --tax {input.taxonomy:q} \
            --split_char {params.split_char:q} \
            --output {output.out:q} \
            --tax_output {output.tax:q} \
            --mapping {input.mapping:q}
        """

rule get_binary:
    input: 
        out=rules.get_taxonomy.output.out,
        tax=rules.get_taxonomy.output.tax,
        script="workflow/scripts/make_binary.R"
    output: "results/{database}/binary/get_binary/{mapping_db}-binary.rds"
    conda: "../envs/matrix.yml"
    shell:
        """
        Rscript {input.script:q} {input.out:q} {input.tax:q} {output:q}
        """

rule get_continuous:
    input:
        hits=rules.combine_species_hits.output.tophits,
        mapping=lambda wildcards: config["files"]["mapping"][wildcards.mapping_db],
        taxonomy=lambda wildcards: config["files"]["taxonomy"][wildcards.mapping_db],
        script="workflow/scripts/sparse_continuous_pangenomes.R"
    output:
        rds="results/{database}/binary/get_binary/get_continuous/{mapping_db}-continuous.rds",
        csv="results/{database}/binary/get_binary/get_continuous/{mapping_db}-continuous-full.csv"
    params:
        split_char=config["create_species_matrix"]["split_char"]
    threads: 8
    resources:
        mem_mb=64000
    conda: "../envs/matrix.yml"
    shell:
        """
        Rscript {input.script:q} \
            --input {input.mapping:q} \
            --genome_metadata {input.taxonomy:q} \
            --combined_hits {input.hits:q} \
            --output_rds {output.rds:q} \
            --output_csv {output.csv:q} \
            --split_char {params.split_char:q} \
            --threads {threads} --memory 16 --speedy TRUE
        """


rule get_tree:
    input: 
        tax="results/{database}/binary/get_taxonomy/{mapping_db}-taxonomy.csv",
        tree=lambda wildcards: config["files"]["tree"][wildcards.mapping_db],
        alttree=lambda wildcards: [config["files"]["alttree"][wildcards.mapping_db]["tree"]]
        if config["files"].get("alttree", {}).get(wildcards.mapping_db, {}).get("tree") else [],
        script="workflow/scripts/make_tree.R"
    output: "results/{database}/binary/get_tree/{mapping_db}-tree.rds"
    params:
        alttree_args=lambda wildcards, input: shlex.join([
            "-a", *input.alttree
        ]) if input.alttree else ""
    conda: "../envs/matrix.yml"
    shell:
        """
        Rscript {input.script:q} -i {input.tree:q} {params.alttree_args} -t {input.tax:q} -o {output:q}
        """

