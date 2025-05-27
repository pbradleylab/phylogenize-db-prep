from glob import glob
import os, peppy

include: "resources.smk"
include: "matrix.smk"


def get_all_links(wildcards):
    return glob.glob(rules.link_information.output+"*.tsv")

def get_chunks(wildcards):
    chunk_dir=checkpoints.break_fasta_apart.get(
        database=wildcards.database,
        mapping_db=wildcards.mapping_db).output[0]

    # extract chunks
    chunks=glob_wildcards(os.path.join(chunk_dir, "/*.fasta"))
    print(chunks)
    #print(rules.break_fasta_apart.output+"/{chunk}.")

def get_annotations(wildcards):
    chunk_dir=checkpoints.break_fasta_apart.get(
        database=wildcards.database,
        mapping_db=wildcards.mapping_db).output[0]

    # extract chunks
    chunks=glob_wildcards(os.path.join(chunk_dir, "*.faa"))
    out=expand(rules.get_headers.output, database=wildcards.database, mapping_db=wildcards.mapping_db, chunk=chunks)
    return(out)

checkpoint break_fasta_apart:
    input:rules.mmseqs2_linclust.output.fasta
    output:directory("results/{database}/annotation/break_fasta_apart/{mapping_db}")
    conda:"../envs/annotation.yml"
    shell:
        """
        seqkit split {input} -s 10000 -O {output}
        """

rule get_headers:
    input:get_chunks
    output:"results/{database}/annotation/get_headers/{mapping_db}/{chunk}_headers.txt"
    shell:
        """
        echo {wildcards.chunk}
        #grep '>' {input.fasta} | sed 's/>//g' > {output}
        #awk -F, '{{print FNR,$0}}' OFS='\t' {output} > tmpFile && mv tmpFile {output}
        """

rule run_anntoation:
    input:get_annotations
    output:"results/{database}/run_anntoation/{mapping_db}.done"
    shell:
        """
        echo "hit"
        """

# -------





rule make_gene_calls:
    input: 
        fasta=rules.break_fasta_apart.output[0]+"/{chunk}.faa",
        headers=rules.get_headers.output
    output:"results/{database}/annotation/make_gene_calls/{chunk}.tsv"
    params:
        indir=rules.break_fasta_apart.output
    shell:
        """
        bash make_external_gene_calls.sh {params.indir}/{wildcards.chunk}.faa {output}
        """
# ======


rule anvio_contigs_db:
    input:
        chunk=get_chunks,
        gene_calls=rules.make_gene_calls.output
    output:"results/{project}/annotation/anvio_contigs_db/{chunk}.db"
    params:
        indir=rules.break_fasta_apart.output
    conda:"../envs/annotation.yml"
    threads: config["anvio"]["contigs_db"]["threads"]
    shell:
        """
        anvi-gen-contigs-database -T {threads} -o {output} -f {params.indir}/{wildcards.chunk}.faa --external-gene-calls {input.gene_calls}
        """

rule anvio_kegg_annotation:
    input:
        db=rules.anvio_contigs_db.output,
        kofam=rules.anvio_setup_kegg_kofams.output
    output:"results/temp/{database}/annotation/anvio_kegg_annotation/{mapping_db}/{chunk}_anvio_run_kegg_kofams.0"
    conda:"../envs/anvio.yml"
    log: "logs/{database}/annotation/anvio_run_kegg_kofams/{mapping_db}_{chunk}.log"
    threads: config["anvio"]["run_kegg_annotations"]["threads"]
    shell:
        """
        anvi-run-kegg-kofams -c {input.db} --kegg-data-dir {input.kofam} -T {threads} --just-do-it 2> {log}
        touch {output}
        """

rule anvio_export_functions:
    input:
        kegg=rules.anvio_kegg_annotation.output,
        db=rules.anvio_contigs_db.output,
    output: "results/{database}/annotation/anvio_export_functions/{mapping_db}_{chunk}.tsv"
    conda:"../envs/annotation.yml"
    shell:
        """
        anvi-export-functions -c {input.db} -o {output}
        """

#----

rule link_information:
    input:
        linker=rules.get_headers.output,
        gene_calls=rules.make_gene_calls.output,
        functions=rules.anvio_export_functions.output
    output:"results/{database}/annotation/link_information/{chunk}_merged.tsv"
    shell:
        """
        python merge_annotations.py {input.functions} {input.gene_calls} {input.linker} {output}
        """

rule combine_link_info:
    input:get_all_links
    output:"results/{database}/annotation/combine_link_info/ko_merged.tsv"
    shell:
        """
        {head -n 1 {input[0]};
        for f in {input[1:]}; do
            tail -n +2 "$f";
        done} > {output}
        """



rule link_sequences:
    input:rules.combine_link_info.output
    output:"results/{database}/annotation/link_sequences/annotations.tsv"
    conda:"../envs/annotation.yml"
    params:
        node_heads="dummy"#config["node_heads"]
    shell:
        """
        python sequence_to_annotation_linker.py {input} {params.node_heads} {output}
        """

rule link_nodes:
    input:rules.link_sequences.output
    output:"results/{database}/annotation//link_nodes/node_annotations.tsv"
    conda:"../envs/annotation.yml"
    params:
        node_heads="dummy"#config["node_heads"]
    shell:
        """
        python sequence_to_annotation_linker.py {params.node_heads} {input} {output}
        """
