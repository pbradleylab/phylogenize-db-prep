include: "clustering.smk"


def get_16s(wildcards):
    return(config["files"]["16S"][wildcards.mapping_db]["fasta"])

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
        
        cat {input.clustered} | cut -f1  > /tmp/tmp.2
        cat {input.clustered} | cut -f2  > /tmp/tmp.3
        paste /tmp/tmp.2 /tmp/tmp.3 >> /tmp/{wildcards.mapping_db}.1 

        # Ensure only one gene hit is used per annotation and in the order of databases selected.
        awk '!seen[$1]++' /tmp/{wildcards.mapping_db}.1 | cut -f1,2 >> {output.tophits}
        rm /tmp/{wildcards.mapping_db}.1 /tmp/tmp.2 /tmp/tmp.3

        cat {input.allhits} > {output.allhits}
        """

# Create the taxonomy file and generate the input for making the binary file needed
# to pass into R to compress and spearate into a list of binaries per phylum.
rule get_taxonomy:
    input: rules.combine_species_hits.output.tophits
    output: 
        out=temporary("results/{database}/binary/get_taxonomy/{mapping_db}-binary-temp.csv"),
        tax="results/{database}/binary/get_taxonomy/{mapping_db}-taxonomy.csv"
    params:
        split_char=config["create_species_matrix"]["split_char"],
        key="{mapping_db}"
    conda: "../envs/matrix.yml"
    shell:
        """
        mapping=$(python -c 'import json; print(json.load(open("config/config.json"))["files"]["mapping"]["{params.key}"])')
        tax=$(python -c 'import json; print(json.load(open("config/config.json"))["files"]["taxonomy"]["{params.key}"])')
        
        python workflow/scripts/reformat_taxonomy.py \
            --input {input} \
            --tax $tax \
            --split_char {params.split_char} \
            --output {output.out} \
            --tax_output {output.tax} \
            --mapping $mapping
        """

rule get_binary:
    input: 
        out=rules.get_taxonomy.output.out,
        tax=rules.get_taxonomy.output.tax
    output: "results/{database}/binary/get_binary/{mapping_db}-binary.rds"
    conda: "../envs/matrix.yml"
    shell:
        """
        Rscript workflow/scripts/make_binary.R {input.out} {input.tax} {output} 
        """

rule get_tree:
    input: rules.get_taxonomy.output.out
    output: "results/{database}/binary/get_tree/{mapping_db}-tree.rds"
    params:
        tree=lambda wildcards: config["files"]["tree"][wildcards.mapping_db]
    conda: "../envs/matrix.yml"
    shell:
        """
        Rscript workflow/scripts/make_tree.R {params.tree} {output}
        """

rule get_16s:
    input:get_16s
    output: "results/{database}/binary/get_16s/{mapping_db}.faa"
    conda: "../envs/matrix.yml"
    shell:
        """
        echo "to be completed"  
        """
