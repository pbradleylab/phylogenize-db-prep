def get_target(wildcards):
    out=[]
    if config["target_db"]["db"].lower() == "uhgp90":
        out = rules.create_uhgp90_target_db.output.outdir
    elif config["target_db"]["db"] == "UniRef90":
        out = rules.create_uniprot90_target_db.output.uniprot90_path
    return out


rule download_uhgp90:
    output: "resources/uhgp-90.tar.gz"
    params:
       url=config["target_db"]["uhgp90_url"]
    conda:"../envs/target_db.yml"
    shell:
        """
        wget -c {params.url} -O {output}
        """

rule unpack_uhgp90:
    input: rules.download_uhgp90.output
    output: 
        fasta="resources/references/uhgp/uhgp-90/uhgp-90.faa"
    params:
        outdir=directory("resources/references/uhgp/")
    shell:
        """
        tar -zxvf {input} -C {params.outdir}
        """

rule create_uhgp90_target_db:
    input:
        uhgp90_fasta=rules.unpack_uhgp90.output.fasta
    output:
        index="resources/uhgp90/uhgp90.index",
        outdir=directory("resources/uhgp90/")
    params:
        uhgp90_prefix="uhgp90",
        uhgp90_path="resources/uhgp90/"
    log: "logs/uhgp90/createdb/uhgp90.log"
    conda: "../envs/target_db.yml"
    threads: config["mmseqs2"]["createdb"]["threads"]
    shell:
        """
        mmseqs createdb {input} {params.uhgp90_path}/{params.uhgp90_prefix} --dbtype 1 2> {log}
        mmseqs createindex {params.uhgp90_path}/{params.uhgp90_prefix} /tmp 2> {log}
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
