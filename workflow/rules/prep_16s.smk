
rule get_16s:
    input:
      fa=lambda wc: config["files"]["16S"][wc.mapping_db]["fasta_dir"],
      md=lambda wc: config["files"]["taxonomy"][wc.mapping_db]
    output: "results/{database}/16S/initial/{mapping_db}.faa"
    conda: "../envs/16S.yml"
    shell:
        """
            scripts/extract_ssu_from_gff.py -i {input.fa} \
               -l "results/{wildcards.database}/ssu.log" \
               -m {input.md} \
               -o {output}
        """

rule make_db:
    input: "results/{database}/16S/initial/{mapping_db}.faa"
    output: "results/{database}/16S/initial/{mapping_db}.udb"
    conda: "../envs/16S.yml"
    shell: """
            vsearch --makeudb_usearch {input} --output \
               results/{wildcards.database}/16S/initial/{wildcards.mapping_db}.udb 
           """

rule all_v_all_16s:
    input:
        fa="results/{database}/16S/initial/{mapping_db}.faa",
        db="results/{database}/16S/initial/{mapping_db}.udb"
    output: "results/{database}/16S/initial/{mapping_db}_all_v_all.txt"
    conda: "../envs/16S.yml"
    shell:
        "vsearch --usearch_global {input.fa} --db {input.db} --id 0.95 " + \
        "--maxaccepts=20 --blast6out {output}"

rule filter_results:
    input:
        ava="results/{database}/16S/initial/{mapping_db}_all_v_all.txt",
        fa="results/{database}/16S/initial/{mapping_db}.faa",
        md=lambda wc: config["files"]["taxonomy"][wc.mapping_db]
    output: "results/{database}/16S/filtered/{mapping_db}.faa"
    shell: """
        scripts/filter_all_v_all.R -i {input.ava} -f {input.fa} -m {input.md} -o {output}
    """
