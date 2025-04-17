include: "clustering.smk"


# Combines species with a 90% or greater identity match to the target database, 
# and the unmapped regions to a list of species specific vectors by their centroid.
rule combine_species_hits:
    input:
        tophits=expand(rules.get_top_50_evals.output.tophits,
                       target_db=TARGET_DBS,
                       mapping_db=config["annotation"]["mapping_databases"].keys(),
                       database=config["database"]),
        clustered=rules.mmseqs2_linclust.output.tsv,
        # Wait for all database checkpoints to complete
        checkpoints=expand(rules.database_processing_checkpoint.output,
                          target_db=TARGET_DBS,
                          database=config["database"])
    output:
         txt="results/{database}/binary/combined_species_hits/{database}.txt",
         outdir=directory("results/{database}/binary/combine_species_hits/combined_species_hits/")
    shell:
        """
        echo -e "query\ttarget\tpident\tdatabase" > {output.txt}
        cat {input.clustered} | cut -f1  > /tmp/tmp.1
        cat {input.clustered} | cut -f2  > /tmp/tmp.2
        paste /tmp/tmp.2 /tmp/tmp.1 >> {output.txt} 
        # Combine all tophits files with appropriate headers
        for file in {input.tophits}; do
            db_info=$(basename $file | sed 's/_convertlis_tophits.tsv//')
            awk -v db="$db_info" '{{print $0"\t"db}}' $file >> {output.txt}
        done
        """

# Create a binary matrix of the hits > 50% and those that unmapped 
# but clustered.
rule create_species_matrix:
    input: rules.combine_species_hits.output.outdir
    output: 
        out="results/{database}/binary/final/species_matrix/{database}.txt",
        dup="results/{database}/binary/species_matrix/{database}.dups"
    conda: "../envs/matrix.yml"
    shell:
        """
        python workflow/scripts/combine_species.py \
            --output {output.out} \
            --dir {input} \
            --ext ".txt" \
            --duplicates {output.dup}
        """
