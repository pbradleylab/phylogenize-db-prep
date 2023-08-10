include: "query_db.smk"
include: "target_db.smk"
include: "translation.smk"

# Map the amino acid sequences by similarity in the uniref50 database
# The internal prefilter module is called which is high sensitivity to
# detect high scores and ungapped alignment. This could be exchanged for
# the `mmseqs search` command for lower sensitivity.
#
# Note: We call params here from previous rules. This method is continued
#    throughout all subsequent methods for continuety. These can be abstracted
#    out to a config however it looses some of the automation that way.
rule mmseqs2_map_uniref50:
     input:
         query=rules.create_mmseqs2_query_db.output.query_path,
         target=rules.create_uniref50.output.uniref50_path
     output:
         outdir=directory("results/{database}/mmseqs2/uniref50/mapping/"),
         index="results/{database}/mmseqs2/uniref50/mapping/{database}_map.index"
     log: "logs/{database}/uniref50/mmseqs2/mapping/mmseqs2_map.log"
     params:
         prefix="{database}_map",
         query_prefix=rules.create_mmseqs2_query_db.params.query_prefix,
         target_prefix="UniRef50"
     threads: config["mmseqs2"]["map"]["threads"]
     conda: "../envs/mapping.yml"
     shell:
         """
         mmseqs map -a --threads {threads} {input.query}/{params.query_prefix} \
            {input.target}/{params.target_prefix} {output.outdir}/{params.prefix} \
            --min-seq-id 0.50 /tmp 2> {log}
         """

# Converts the database's mappings to a sam format. The mapped (aligned)
# sequences are then taken to generate a new database in rule:
# `samtools_get_aligned`.
rule mmseqs2_convertalis_sam_uniref50:
     input:
         query=rules.create_mmseqs2_query_db.output.query_path,
         target=rules.create_uniref50.output.uniref50_path,
         mapped=rules.mmseqs2_map_uniref50.output.outdir
     output: "results/{database}/uniref50/mmseqs2/convertalis/{database}_convertlis.sam"
     params:
         prefix=rules.mmseqs2_map_uniref50.params.prefix,
         query_prefix=rules.create_mmseqs2_query_db.params.query_prefix,
         target_prefix="UniRef50"
     threads: config["mmseqs2"]["convertalis"]["threads"]
     conda: "../envs/mapping.yml"
     shell:
         """
         mmseqs convertalis \
             {input.query}/{params.query_prefix} \
             {input.target}/{params.target_prefix} \
             {input.mapped}/{params.prefix} \
             {output} --format-mode 1
         """

# Get all the sequences that are not labeled as aligned in the .sam
# file. The output is a bam as to retain headers as samtools doesn't
# keep header information which results in many downstream errors for
# programs.
rule samtools_aligned_bam_uniref50:
     input: rules.mmseqs2_convertalis_sam_uniref50.output,
     output: "results/{database}/uniref50/samtools/aligned/{database}_aligned.bam",
     log: "logs/{database}/uniref50/samtools/mapping/{database}_map.log"
     conda: "../envs/mapping.yml"
     shell:
         """
         samtools view -b {input} -o {output} 2> {log}
         """

# Get only the aligned sequences.
# Convert the aligned samples to a fasta to build a new database.
# We consider these to be temporary files as we retain the bam for QC
# afterwards since reads may have failed to align due to other reasons
# and not solely based on species sequence homology.
rule samtools_fasta_uniref50:
     input: rules.samtools_aligned_bam_uniref50.output
     output: "results/{database}/uniref50/samtools/fasta/{database}.fasta"
     conda: "../envs/mapping.yml"
     shell:
         """
         samtools fasta {input} > {output}
         """

rule get_aligned_sequences:
    input:
         mapped=rules.samtools_fasta_uniref50.output,
         all_sequences=rules.combine_fasta_uniref50.output
    output: "results/{database}/uniref50/samtools/fasta/{database}_unmapped.fasta"
    conda: "../envs/mapping.yml"
    shell:
         """
         bash workflow/scripts/get_aligned_sequences.sh {input.mapped} {input.all_sequences} {output}
         """
