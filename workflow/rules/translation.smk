def get_transeq_output(wildcards):
    outputLST = []
    for subsample in pep.subsample_table.subsample.tolist():
        project = get_subsample_attributes(subsample, "project", pep)
        # Always run rules on the outside
        outputLST.append(rules.transeq.output[0].format(database=project, pangenome=subsample))
    return outputLST

def get_pangenomes(wildcards):
    pangenomes = get_subsample_attributes(wildcards.pangenome, "pangenomes", pep)
    return pangenomes


# Append the name of the MIDAS identifier to the fron for later parsing.
rule append_name:
    input: get_pangenomes
    output: "results/{database}/renamed/{pangenome}.ffn"
    params:
        prefix="{pangenome}"
    conda: "../envs/translation.yml"
    log: "logs/{database}/renamed/{pangenome}.log"
    shell:
        """
        bbrename.sh in={input} out={output} prefix={params.prefix} addprefix=t 2> {log}
        """

# Translate nucleotides per genome to peptide sequences per genome.
# Please check the config to set if stop codons shouldn't convert
# from the default '*' character to an 'X' representing any animo acid.
rule transeq:
    input: rules.append_name.output
    output: "results/{database}/transeq/{pangenome}.ffn"
    params:
        clean=config["transeq"]["convert_missing_to_x"]
    log: "logs/{database}/transeq/{pangenome}.log"
    conda: "../envs/translation.yml"
    shell:
        """
        transeq {input} {output} -clean {params.clean} 2> {log}
        """
 
# Combine the fasta that are translated to retrieve the unmapped alignments
# in later steps.
rule combine_fasta_uniref50:
     input: get_transeq_output
     output: "results/{database}/combined_uniref/{database}.fa"
     conda: "../envs/translation.yml"
     shell:
         """
         cat {input} > {output}
         """






rule reduce_complexity:
     input: config["indir"]
     output: "results/{database}/combined_fasta/{database}.fa"
     params:
         database="{database}"
     conda:"../envs/translation.yml"
     shell:
         """
         bash workflow/scripts/reduce_complexity.sh {input} {output} {params.database}
         """
