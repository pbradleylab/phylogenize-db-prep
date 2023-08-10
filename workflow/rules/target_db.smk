""" Create the input databases for UHGP50 and UniRef50. UniRef50 will be used to map and cluster
the database and then the sequences that are not found at a 50% sequence identity are mapped
against UHGP50.

By Kathryn Kananen
"""

rule download_uhgp50:
    output: "resources/uhgp-50.tar.gz"
    params:
       url=config["target_db"]["uhgp50_url"]
    conda:"../envs/target_db.yml"
    shell:
        """
        wget -c {params.url} -O {output}
        """

rule unpack_uhgp50:
    input: rules.download_uhgp50.output
    output: 
        fasta="resources/references/uhgp/uhgp-50/uhgp-50.faa"
    params:
        outdir=directory("resources/references/uhgp/")
    shell:
        """
        tar -zxvf {input} -C {params.outdir}
        """

rule create_uhgp50:
    input: rules.unpack_uhgp50.output.fasta
    output: directory("resources/uhgp50")
    params:
        uhgp50_prefix="uhgp50",
        uhgp50_path="resources/uhgp50/"
    log: "logs/uhgp50/createdb/create_uhgp50/uhgp50.log"
    conda: "../envs/target_db.yml"
    threads: config["mmseqs2"]["createdb"]["threads"]
    shell:
        """
        mkdir -p {params.uhgp50_path}
        mmseqs createdb {input} {params.uhgp50_path}/{params.uhgp50_prefix} --dbtype 1 2> {log}
        """

rule index_uhgp50:
    input:
        uhgp50_fasta=rules.unpack_uhgp50.output.fasta
    output:
        index="resources/uhgp50/uhgp50.index"
    params:
        uhgp50_prefix="uhgp50",
        uhgp50_path="resources/uhgp50/"
    log: "logs/uhgp50/createdb/uhgp50/index_uhgp50.log"
    conda: "../envs/target_db.yml"
    threads: config["mmseqs2"]["createdb"]["threads"]
    shell:
        """
        mmseqs createindex {params.uhgp50_path}/{params.uhgp50_prefix} /tmp 2> {log}
        """

rule create_uniref50:
    output:
        uniref50_path=directory("resources/{database}/UniRef50"),
        uniref50_raw=directory("resources/{database}/UniRef50/raw/")
    params:
        uniref50_prefix="UniRef50",
    conda: "../envs/target_db.yml"
    threads: config["mmseqs2"]["createdb"]["threads"]
    shell:
        """
        mmseqs databases {params.uniref50_prefix} \
            {output.uniref50_path}/{params.uniref50_prefix} \
            {output.uniref50_raw} \
            --threads {threads}
        """
