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

rule create_mmseqs2_query_db:
    input: get_mmseqs2_input
    output: "resources/{database}/custom/custom.index"
    params:
        db_name="resources/{database}/custom/custom"
    conda: "../envs/transformation.yml"
    shell:
        """
        mmseqs createdb {input} {params.db_name} --dbtype 1 
        mmseqs createindex {params.db_name} /tmp
	"""

rule create_mmseqs2_target_db:
    output: 
         fasta="resources/{database}/uniprot90/tmp/latest/uniref90.fasta.gz",
         db_folder="resources/{database}/uniprot90"
    conda: "../envs/transformation.yml"
    threads: 32
    shell:
        """
        mmseqs databases UniRef90 UniRef90 {params.db_folder} --threads {threads}
        """

# Compute the ungapped alignment score for the target and query database
# with the highest score being returned for consecutive k-mer matches.
rule mmseqs2_pref:
     input:
         query=rules.create_mmseqs2_query_db.output,
         target=rules.create_mmseqs2_target_db.output
     output:
         out_dir="results/{database}/mmseqs2/pref/"
     params:
         prefix="{database}_pref"
     threads: config["mmseq2"]["threads"]
     conda: "../envs/clustering.yml"
     shell:
         """
         mkdir -p {output}
         mmseqs prefilter --threads {threads} {input.query} {output.target} {output}/{params.prefix}
         """

# Align the amino acid sequences by similarity in the UniProt 90 database
# The internal prefilter module is called which is low sensitivity to 
# detect high scores and ungapped alignment.
rule mmseqs2_map:
     input:
         query=rules.create_mmseqs2_query_db.output,
         target=rules.create_mmseqs2_target_db.output,
     output:
         out_dir="results/{database}/mmseqs2/mapping/",
         index="results/{database}/mmseqs2/mapping/{database}_map.index"
     log: "logs/{database}/mmseqs2/mapping/mmseqs2_map.log"
     params:
         prefix="{database}_map"
     threads: config["mmseq2"]["threads"]
     conda: "../envs/clustering.yml"
     shell:
         """
         mkdir -p {output} 
         mmseqs map --threads {threads} {input.query} {output.query} {output}/{params.prefix} /tmp/tmp \
             --comp-bias-corr 0 --mask 0 -e inf --max-seqs 300 --exact-kmer-matching 1 \
             --spaced-kmer-pattern 110111 -k 5 -a 1 --min-seq-id 1 \
             2> {log}
         """

rule mmseqs2_convertalis:
     input:
         query=rules.create_mmseqs2_query_db.output,
         target=rules.create_mmseqs2_target_db.output,
         mapped=rules.mmseqs_map.output
     output: "results/{database}/mmseqs2/convertalis/{database}_convertlis.sam"
     threads: config["mmseq2"]["threads"]
     conda: "../envs/clustering.yml"
     shell:
         """
         mmseqs convertalis {input.query} {input.target} {input.mapped} {output} --format-mode 1
         """

# Get all the sequences that are not labeled as unaligned in the .sam file
rule samtools_get_aligned:
     input:
         mapped=rules.mmseqs2_convertalis.output
     output: "results/{database}/samtools/unaligned/{database}_aligned.bam"
     log: "logs/{database}/samtools/mapping/{database}_map.log"
     conda: "../envs/transformation.yml"
     shell:
         """
         samtools view -b -F 4 {input} -o {output}
         """ 

rule samtools_fasta:
     input: rules.samtools_get_aligned.output
     output: "results/{database}/samtools/fasta/{database}.fasta"
     conda: "../envs/transformation.yml"
     shell:
         """
         samtools fasta {input} > {output}
         """

rule create_mmseqs2_unaligned_db:
    input: rules.samtools_fasta.output
    output: "resources/{database}/unmapped/unmapped.index"
    params:
        db_name="resources/{database}/unmapped/unmapped"
    conda: "../envs/transformation.yml"
    shell:
        """
        mmseqs createdb {input} {params.db_name} --dbtype 1
        mmseqs createindex {params.db_name} /tmp
        """

# Make a new database based only on unmapped sequences
# mmseqs clusterupdate DB_old DB_new CLU_old DB_updated CLU_updated tmp
# Sanity Check
# mmseqs createtsv DB_updated DB_updated CLU_updated clusters.tsv

# Cluster the unmapped protein sequences from mmseqs' search
rule mmseq2_linclust:
     input: rules.create_mmseqs2_unaligned_db.output
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
         mmseqs linclust {params.db_name} {params.out_dir} {params.tmp_dir} --min-seq-id {params.seq_id_precent} --threads {threads}
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
