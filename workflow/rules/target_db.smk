def get_target(wildcards):
    out=[]
    if config["target_db"]["db"].lower() == "uhgp90":
        out = rules.create_ughp90_target_db.output
    elif config["target_db"]["db"].lower() == "unitprot90":
        out = rules.create_uniprot90_target_db.output
    return out


rule download_ughp90:
    output: "resources/{database}/uhgp90/uhgp-90.tar.gz"
    params:
       url=config["target_db"]["uhgp90_url"]
    conda:"../envs/target_db.yml"
    shell:
        """
        wget -c {params.url} -O {output}
        """

rule unpack_ughp90:
    input: rules.download_ughp90.output
    output: "resources/{database}/uhgp90/uhgp-90/uhgp.fasta"
    shell:
        """
        tar -zxvf {input} 
        """

rule create_ughp90_target_db:
    input:
        ughp90_fasta=rules.unpack_ughp90.output
    output:
        index="resources/{database}/ughp90/ughp90.index"
    params:
        ughp90_prefix="UniRef90",
        ughp90_path="resources/{database}/ughp90"
    conda: "../envs/target_db.yml"
    threads: config["mmseqs2"]["createdb"]["threads"]
    shell:
        """
        mmseqs createdb {input} {params.ughp90_path}/{params.ughp90_prefix} --dbtype 1 2> {log}
        mmseqs createindex {params.ughp90_path}/{params.ughp90_prefix} /tmp 2> {log}
        """

rule create_uniprot90_target_db:
    output:
        uniprot90_fasta="resources/{database}/uniprot90/tmp/latest/uniref90.fasta.gz",
        uniprot90_path="resources/{database}/uniprot90"
    params:
        uniprot90_prefix="UniRef90",
    conda: "../envs/target_db.yml"
    threads: config["mmseqs2"]["createdb"]["threads"]
    shell:
        """
        mmseqs databases UniRef90 {params.uniprot90_prefix} \
            {output.uniprot90_path}/{params.uniprot90_prefix} \
            --threads {threads}
        """