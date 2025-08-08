include: "clustering.smk"



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
            awk -v db="$db_info" '{{print $1"\t"$2"\t"db}}' $file >> /tmp/{wildcards.mapping_db}.1
        done
        # reverse order for clustering
        awk -v db="denovo" '{{print $2"\t"$1"\t"db}}' {input.clustered} >> /tmp/{wildcards.mapping_db}.1 

        # Ensure only one gene hit is used per annotation and in the order of databases selected.
        awk '!seen[$1]++' /tmp/{wildcards.mapping_db}.1 | cut -f1,2 >> {output.tophits}
        rm /tmp/{wildcards.mapping_db}.1 /tmp/tmp.2 /tmp/tmp.3

        cat {input.allhits} > {output.allhits}
        """

# Create the taxonomy file and generate the input for making the binary file needed
# to pass into R to compress and separate into a list of binaries per phylum.
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

rule get_continuous:
    resources:
        mem_mb=16000
    input:
        combined="results/{database}/binary/combined_species_hits/{mapping_db}_{database}.tsv",
        prot_map=config["files"]["mapping"][wildcards.mapping_db],
        genome_md=config["files"]["taxonomy"][wildcards.mapping_db]
    output:
        rds="results/{database}/binary/get_binary/{mapping_db}-continuous.rds",
        tsv="results/{database}/binary/get_binary/{mapping_db}-continuous-full.tsv"
    shell: """
        scripts/sparse_continuous_pangenomes.R -i {input.prot_map} -g {input.genome_md} \
            -c {input.combined} -r {output.rds} -t {output.tsv} -m {resources.mem_mb/1000}
    """

rule get_tree:
    input: "results/{database}/binary/get_taxonomy/{mapping_db}-taxonomy.csv"
    output: "results/{database}/binary/get_tree/{mapping_db}-tree.rds"
    params:
        tree=lambda wildcards: config["files"]["tree"][wildcards.mapping_db]
    conda: "../envs/matrix.yml"
    shell:
        """
        Rscript workflow/scripts/make_tree.R -i {params.tree} -t {input} -o {output}
        """

rule get_16s:
    input: get_16s
    output: "results/{database}/binary/get_16s/{mapping_db}.faa"
    conda: "../envs/matrix.yml"
    shell:
        """
        echo "to be completed"  
        """
