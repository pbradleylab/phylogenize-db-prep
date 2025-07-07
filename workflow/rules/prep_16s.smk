
rule get_16s:
    input:
      fa=config["files"]["16S"][wildcards.mapping_db]["fasta_dir"],
      md=config["files"]["taxonomy"][wildcards.mapping_db]
    output: "results/{database}/16S/initial/{mapping_db}.faa"
    conda: "../envs/16S.yml"
    run:
        log_file = f"results/{wildcards.database}/ssu.log"
        shell("""
            scripts/extract_ssu_from_gff.py -i {input.fa} \
               -l {log_file}
               -m {input.md}
               -o {output}
        """)

rule 16s_all_v_all:
    input: "results/{database}/16S/initial/{mapping_db}.faa"
    output: "results/{database}/16S/initial/{mapping_db}_all_v_all.txt"
    conda: "../envs/16S.yml"
    run:
        dbfile = f"results/{wildcards.database}/16S/initial/{wildcards.mapping_db}.udb"
        shell("""
            vsearch --makeudb_usearch {input} --output {dbfile}
            vsearch --usearch_global {input} --db {dbfile} --id 0.95 --maxaccepts=20 \ 
               --blast6out {output}
         """)

rule filter_results:
    input:
        ava="results/{database}/16S/initial/{mapping_db}_all_v_all.txt",
        fa="results/{database}/16S/initial/{mapping_db}.faa",
        md=config["files"]["taxonomy"][wildcards.mapping_db]
    output: "results/{database}/16S/filtered/{mapping_db}.faa"
    shell: """
        scripts/filter_all_v_all.R -i {input.ava} -f {input.fa} -m {input.md} -o {output}
    """