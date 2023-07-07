include: "query_db.smk"
include: "target_db.smk"
include: "translation.smk"

# Map the amino acid sequences by similarity in the UniProt 90 database
# The internal prefilter module is called which is high sensitivity to
# detect high scores and ungapped alignment. This could be exchanged for
# the `mmseqs search` command for lower sensitivity.
#
# Note: We call params here from previous rules. This method is continued
#    throughout all subsequent methods for continuety. These can be abstracted
#    out to a config however it looses some of the automation that way.
rule mmseqs2_map:
     input:
         query=rules.create_mmseqs2_query_db.output.query_path,
         target=rules.create_mmseqs2_target_db.output.uniprot90_path
     output:
         out_dir=directory("results/{database}/mmseqs2/mapping/"),
         index="results/{database}/mmseqs2/mapping/{database}_map.index"
     log: "logs/{database}/mmseqs2/mapping/mmseqs2_map.log"
     params:
         prefix="{database}_map",
         query_prefix=rules.create_mmseqs2_query_db.params.query_prefix,
         target_prefix=rules.create_mmseqs2_target_db.params.uniprot90_prefix
     threads: config["mmseqs2"]["map"]["threads"]
     conda: "../envs/transformation.yml"
     shell:
         """
         mkdir -p {output.out_dir}
            mmseqs map --threads {threads} {input.query}/{params.query_prefix} \
            {input.target}/{params.target_prefix} {output.out_dir}/{params.prefix} \
            /tmp/tmp -a --min-seq-id 0.0 2> {log}
         """

# Converts the database's mappings to a sam format. The unmapped (unaligned)
# sequences are then taken to generate a new database in rule:
# `samtools_get_aligned`.
rule mmseqs2_convertalis_sam:
     input:
         query=rules.create_mmseqs2_query_db.output.query_path,
         target=rules.create_mmseqs2_target_db.output.uniprot90_path,
         mapped=rules.mmseqs2_map.output.out_dir
     output: "results/{database}/mmseqs2/convertalis/{database}_convertlis.sam"
     params:
         prefix=rules.mmseqs2_map.params.prefix,
         query_prefix=rules.create_mmseqs2_query_db.params.query_prefix,
         target_prefix=rules.create_mmseqs2_target_db.params.uniprot90_prefix
     threads: config["mmseqs2"]["convertalis"]["threads"]
     conda: "../envs/database_management.yml"
     shell:
         """
         mmseqs convertalis \
             {input.query}/{params.query_prefix} \
             {input.target}/{params.target_prefix} \
             {input.mapped}/{params.prefix} \
             {output} --format-mode 1
         """
         
# Get all the sequences that are not labeled as unaligned in the .sam
# file. The output is a bam as to retain headers as samtools doesn't
# keep header information which results in many downstream errors for
# programs.
rule samtools_aligned_bam:
     input: rules.mmseqs2_convertalis_sam.output,
     output: "results/{database}/samtools/unaligned/{database}_unaligned.bam",
     log: "logs/{database}/samtools/mapping/{database}_map.log"
     conda: "../envs/transformation.yml"
     shell:
         """
         samtools view -b {input} -o {output} 2> {log}
         """

# Get only the unaligned sequences.
# Convert the unaligned samples to a fasta to build a new database.
# We consider these to be temporary files as we retain the bam for QC
# afterwards since reads may have failed to align due to other reasons
# and not solely based on species sequence homology.
rule samtools_fasta:
     input: rules.samtools_aligned_bam.output
     output: "results/{database}/samtools/fasta/{database}.fasta"
     conda: "../envs/transformation.yml"
     shell:
         """
         samtools fasta {input} > {output}
         """

rule get_aligned_sequences:
    input:
         mapped=rules.samtools_fasta.output,
         all_sequences=rules.combine_fasta.output
    output: "results/{database}/samtools/fasta/{database}_unmapped.fasta"
    conda: "../envs/transformation.yml"
    shell:
         """
         bash workflow/scripts/get_unaligned_sequences.sh {input.mapped} {input.all_sequences} {output}
         """