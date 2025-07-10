rule get_16s:
    input:
      fa=lambda wc: config["files"]["16S"][wc.mapping_db]["fasta_dir"],
      md=lambda wc: config["files"]["taxonomy"][wc.mapping_db]
    output: "results/{database}/16S/initial/{mapping_db}.fna"
    conda: "../envs/16S.yml"
    threads: 16
    shell:
        """
            scripts/extract_ssu_from_gff.py -i {input.fa} \
               -l "results/{wildcards.database}/ssu.log" \
               -m {input.md} \
               -o {output} \
               -n {threads}
        """

rule make_db:
    input: "results/{database}/16S/initial/{mapping_db}.fna"
    output: "results/{database}/16S/initial/{mapping_db}.udb"
    conda: "../envs/16S.yml"
    threads: 16
    shell: """
            vsearch --makeudb_usearch {input} --output \
               results/{wildcards.database}/16S/initial/{wildcards.mapping_db}.udb \
               --threads {threads}
           """

rule all_v_all_16s:
    input:
        fa="results/{database}/16S/initial/{mapping_db}.fna",
        db="results/{database}/16S/initial/{mapping_db}.udb"
    output: "results/{database}/16S/initial/{mapping_db}_all_v_all.txt"
    conda: "../envs/16S.yml"
    threads: 16
    shell:
        """
            vsearch --usearch_global {input.fa} --db {input.db} --id 0.95 \
            --maxaccepts=20 --blast6out {output} --threads {threads}
        """

rule tax_filter_results:
    input:
        ava="results/{database}/16S/initial/{mapping_db}_all_v_all.txt",
        fa="results/{database}/16S/initial/{mapping_db}.fna",
        md=lambda wc: config["files"]["taxonomy"][wc.mapping_db]
    output: "results/{database}/16S/tax_filtered/{mapping_db}.fna"
    conda: "../envs/16S.yml"
    shell: """
        scripts/filter_all_v_all.R -i {input.ava} -f {input.fa} -m {input.md} -o {output}
    """

# Filter by length (sequence must be at least half the upper 95th percentile to get counted)
rule len_filter_results:
    input: "results/{database}/16S/tax_filtered/{mapping_db}.fna"
    output: "results/{database}/16S/len_filtered/{mapping_db}.fna"
    conda: "../envs/16S.yml"
    shell: """
        scripts/fasta_length_filter.R -i {input} -o {output} -f 0.5 -u 0.95
    """

rule make_16S_tree:
    input: "results/{database}/16S/len_filtered/{mapping_db}.fna"
    output: directory("results/{database}/16S/{mapping_db}-pasta/")
    conda: "../envs/16S.yml"
    shell: "run_pasta.py -i {input} -o {output}"
