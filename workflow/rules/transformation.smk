""" Add rules to this section that are related to quality control and post processing.
"""
from scripts.utils import *
configfile: "config/config.json"

def get_pangenomes(wildcards):
    pangenomes = get_subsample_attributes(wildcards.pangenome, "reads", pep)
    return pangenomes

def get_mmseq2_input(wildcards):
    outputLST = []
    for subsample in pep.subsample_table.subsample.tolist():
        project = get_subsample_attributes(subsample, "project", pep)
        # Always run rules on the outside
        outputLST.append(rules.transeq.output[0].format(database=project, pangenome=subsample))
    return outputLST

# Translate nucleotides per genome to peptide sequences per genome.
rule transeq:
    input:get_pangenomes
    output: "results/{database}/transeq/{pangenome}.ffn"
    conda: "../envs/transformation.yml"
    shell:
        """
        transeq {input} {output}
        """

# Map the protein sequences to 
rule create_mmseq2_query_db:
    input: get_mmseq2_input
    output: "resources/{database}/custom/custom.index"
    params:
        db_name="resources/{database}/custom/custom"
    conda: "../envs/transformation.yml"
    shell:
        """
        mmseqs createdb {input} {params.db_name} --dbtype 1 
        mmseqs createindex {params.db_name} /tmp
	"""

rule create_mmseq2_target_db:
    output: "resources/{database}/uniprot90/tmp/latest/uniref90.fasta.gz"
    params:
         db_name="resources/{database}/uniprot90"
    conda: "../envs/transformation.yml"
    threads: 8
    shell:
        """
        mmseqs databases UniRef90 UniRef90 {params.db_name}
        """

rule mmseq2_query:
     input: 
         rules.create_mmseq2_query_db.output,
         rules.create_mmseq2_target_db.output
     output: "results/mmseq2/db/{database}/"
     params:
         query="resources/{database}/custom/custom",
         target="resources/{database}/uniprot90/UniRef90",
         db="{database}"
     conda: "../envs/transformation.yml"
     shell:
         """
         mmseqs search {params.query} {params.target} {output}/{params.db} tmp --num-iterations 2
	 """

rule mmseq2_convert_blast:
     input: rules.mmseq2_query.output
     output: "results/mmseq2/convertalis/{database}/"
     conda: "../envs/transformation.yml"
     params:
         query="resources/{database}",
         target="resources/UniRef90",
         db="{database}"
     shell:
         """
         mmseqs convertalis {params.query} {params.target} {output}/{params.db} {input}
         """

# Cluster the unmapped protein sequences from mmseqs' search
rule mmseq2_cluster:
     input: rules.mmseq2_query.output
     output: 
         database="results/{database}/mmseq2/{database}.dbtype"
     params:
         db_name="resources/{database}",
         out_dir="results/{database}/mmseq2/{database}",
         seq_id_precent=config["mmseq2"]["seq_id_precent"],
         tmp_dir=config["mmseq2"]["tmp_dir"]
     conda: "../envs/transformation.yml"
     threads: config["mmseq2"]["threads"]
     shell:
         """
         mmseqs cluster {params.db_name} {params.out_dir} {params.tmp_dir} --min-seq-id {params.seq_id_precent} --threads {threads}
         """




# rule peptide_matrix_generation:
#     input:rules.build_pangenome_database.output
#     output: 
#     params:
#     log: 
#     resources: 
#     conda: "../envs/transformation.yml"
#     shell:
#         """
#         """
