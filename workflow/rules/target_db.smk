rule create_mmseqs2_target_db:
    output:
        uniprot90_fasta="resources/{database}/uniprot90/tmp/latest/uniref90.fasta.gz",
        uniprot90_path="resources/{database}/uniprot90"
    params:
        uniprot90_prefix="UniRef90",
    conda: "../envs/database_management.yml"
    threads: config["mmseqs2"]["createdb"]["threads"]
    shell:
        """
        mmseqs databases UniRef90 {params.uniprot90_prefix} \
            {output.uniprot90_path}/{params.uniprot90_prefix} \
            --threads {threads}
        """