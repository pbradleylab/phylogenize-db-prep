include: "blast.smk"

rule get_top_90_evals:
     input: rules.mmseqs2_convertalis_blast.output
     output: "results/{database}/mmseqs2/top_90/{database}_convertlis.tsv"
     conda: "../envs/matrix.yml"
     shell:
         """
         awk '$3>90 {{print}}' {input} > {output}
         """

# Combines species with a 90% or greater identity match to the target database, 
# and the unmapped regions to a list of species specific vectors by their centroid.
rule combine_species_hits:
    input:
        unmapped=rules.mmseqs2_convertalis_unmapped_blast.output,
        identity_90=rules.get_top_90_evals.output
    output: 
         txt="results/{database}/mmseqs2/combined_species_hits/{database}.txt",
         out_dir=directory("results/{database}/mmseqs2/combined_species_hits/")
    conda: "../envs/matrix.yml"
    log: "logs/{database}/mmseqs2/hits_90/mmseqs2_hits_90.log"
    shell:
        """
        cat {input.unmapped} | cut -f1,2 > {output.txt}
        cat {input.identity_90} | cut -f1,2 | sed '1d' >> {output.txt}
        """

# Create a binary matrix of the hits > 90% and those that unmapped 
# but clustered.
rule create_species_matrix:
    input: rules.combine_species_hits.output.out_dir
    output: "results/{database}/final/species_matrix/{database}.txt"
    conda: "../envs/matrix.yml"
    log: "logs/{database}/mmseqs2/hits_90/mmseqs2_hits_90.log"
    shell:
        """
        python workflow/scripts/combine_species.py \
            --output {output} \
            --dir {input} \
            --ext ".txt"
        """
