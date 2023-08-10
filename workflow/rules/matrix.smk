include: "blast.smk"


rule get_top_50_evals_uniref50:
     input: rules.mmseqs2_convertalis_blast_uniref50_db.output
     output: "results/{database}/uniref50/mmseqs2/top_50/{database}_convertlis.tsv"
     conda: "../envs/matrix.yml"
     shell:
         """
         awk '$3>50 {{print}}' {input} > {output}
         """

rule get_top_50_evals_uhgp50:
     input: rules.mmseqs2_convertalis_blast_uhgp50_db.output
     output: "results/{database}/uhgp50/mmseqs2/top_50/{database}_convertlis.tsv"
     conda: "../envs/matrix.yml"
     shell:
         """
         awk '$3>50 {{print}}' {input} > {output}
         """

# Combines species with a 90% or greater identity match to the target database, 
# and the unmapped regions to a list of species specific vectors by their centroid.
rule combine_species_hits:
    input:
        uhgp50_identity_50=rules.get_top_50_evals_uhgp50.output,
        unmapped=rules.mmseqs2_convertalis_unmapped_blast_uhgp50_db.output,
        uniref50_identity_50=rules.get_top_50_evals_uniref50.output
    output: 
         txt="results/{database}/uniref50/mmseqs2/combined_species_hits/{database}.txt",
         outdir=directory("results/{database}/uniref50/mmseqs2/combined_species_hits/")
    conda: "../envs/matrix.yml"
    log: "logs/{database}/uniref50/mmseqs2/hits_50/mmseqs2_hits_50.log"
    shell:
        """
        cat {input.unmapped} | cut -f1,2 > {output.txt}
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
    output: "results/{database}/final/species_matrix/{database}.txt"
    conda: "../envs/matrix.yml"
    log: "logs/{database}/mmseqs2/hits_50/mmseqs2_hits_50.log"
    shell:
        """
        python workflow/scripts/combine_species.py \
            --output {output} \
            --dir {input} \
            --ext ".txt"
        """
