
rule get_16s:
    input: config["files"]["16S"][wildcards.mapping_db]["fasta_dir"]
    output: "results/{database}/binary/get_16s/{mapping_db}.faa"
    conda: "../envs/16S.yml"
    run:
        log_file = f"results/{wildcards.database}/ssu.log"
        md_file = config["files"]["taxonomy"][wildcards.mapping_db]
        shell("""
            scripts/extract_ssu_from_gff.py -i {input} \
               -l {log_file}
               -m {md_file}
               -o {output}
        """)
