include: "transformation.smk"


def get_transeq_output(wildcards):
    outputLST = []
    for subsample in pep.subsample_table.subsample.tolist():
        project = get_subsample_attributes(subsample, "project", pep)
        # Always run rules on the outside
        outputLST.append(rules.transeq.output[0].format(database=project, pangenome=subsample))
    return outputLST
        

# Combine the fasta that are translated to retrieve the unmapped alignments
# in later steps.
rule combine_fasta:
     input: get_transeq_output
     output: "results/{database}/combine_fasta/{database}.fa"
     conda: "../envs/transformation.yml"
     shell:
         """
         cat {input} > {output}
         """

rule get_unaligned_sequences:
    input:
         mapped=rules.samtools_fasta.output,
         all_sequences=rules.combine_fasta.output
    output: "results/{database}/samtools/fasta/{database}_unmapped.fasta"
    conda: "../envs/transformation.yml"
    shell:
         """
         bash workflow/scripts/get_unaligned_sequences.sh {input.mapped} {input.all_sequences} {output}
         """

rule get_top_90_evals:
     input: rules.mmseqs2_convertalis_blast.output
     output: "results/{database}/mmseqs2/top_90/{database}_convertlis.tsv"
     conda: "../envs/transformation.yml"
     shell:
         """
         awk '$3>90 {{print}}' {input} > {output}
         """
