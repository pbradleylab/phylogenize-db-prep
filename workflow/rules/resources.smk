""" Resources retrieval and any rule relating to the transformation of resource files    
should be placed here. 
"""
configfile: "config/config.json"


def retrieve_databases(wildcards):
    databaseLST = [] 
    if config["target_db"]["uniprot90"].lower() == "true":
        databaseLST.append(rules.create_mmseqs2_target_db.output.fasta.format(database=database))
    if config["target_db"]["uhgp90"].lower() == "true":
        databaseLST.append(rules.unpack_ughp90.output.fasta.format(database=database))

def get_mmseqs2_input(wildcards):
    outputLST = []
    for subsample in pep.subsample_table.subsample.tolist():
        project = get_subsample_attributes(subsample, "project", pep)
        # Always run rules on the outside
        outputLST.append(rules.transeq.output[0].format(database=project, pangenome=subsample))
    return outputLST
    

rule download_ughp90:
    output: "resources/{database}/uhgp90/uhgp-90.tar.gz"
    params:
       url=config["target_db"]["uhgp90_url"]
    conda: "../envs/resources.yml"
    shell:
        """
        wget -c {params.url} -O {output}
        """

rule unpack_ughp90:
    input: rules.download_ughp90.output
    output: "resources/{database}/uhgp90/uhgp-90"
    shell:
        """
        tar -zxvf {input} 
        """

rule create_mmseqs2_target_db:
    output:
         fasta="resources/{database}/uniprot90/tmp/latest/uniref90.fasta.gz",
         target_path="resources/{database}/uniprot90"
    params:
        target_prefix="UniRef90"
    conda: "../envs/transformation.yml"
    threads: config["mmseqs2"]["createdb"]["threads"]
    shell:
        """
        mmseqs databases UniRef90 {params.target_prefix} {output.target_path}/{params.target_prefix} --threads {threads}
        """

# Creates a query database, query being the database containing the
# pangenomes that is being made into a final species level protein
# binary for Phylogenize.
rule create_mmseqs2_query_db:
    input: get_mmseqs2_input
    output:
        index="resources/{database}/custom/custom.index",
        query_path=directory("resources/{database}/custom")
    params:
        query_prefix="custom"
    conda: "../envs/transformation.yml"
    log: "logs/{database}/custom/{database}.log"
    shell:
        """
        mmseqs createdb {input} {output.query_path}/{params.query_prefix} --dbtype 1 2> {log}
        mmseqs createindex {output.query_path}/{params.query_prefix} /tmp 2> {log}
        """