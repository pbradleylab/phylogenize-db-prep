from glob import glob
import os, peppy
from pathlib import Path

include: "resources.smk"
if annotate_only:
    include: "target.smk"
else:
    include: "matrix.smk"


def get_annotation_input(wildcards, pattern):
    out = pattern.format(database=wildcards.database, mapping_db=wildcards.mapping_db,
        target_db=getattr(wildcards, "target_db", ""))
    if annotate_only and not os.path.isfile(out):
        raise ValueError("annotate_only requires an existing mapping intermediate: " + out +
            ". Run the full workflow first with annotate_only=false.")
    return out


def get_all_links(wildcards):
    chunks=get_chunks(wildcards)
    out=expand(rules.link_information.output,
        database=wildcards.database,
        mapping_db=wildcards.mapping_db,
        chunk=chunks)
    return(out)

def get_chunks(wildcards):
    chunk_dir=checkpoints.break_fasta_apart.get(
        database=wildcards.database,
        mapping_db=wildcards.mapping_db).output[0]
    
    chunks=sorted(p.stem for p in Path(chunk_dir).glob("*.fasta"))
    return(chunks)
    
def get_fasta(wildcards):
    chunks=get_chunks(wildcards)
    out=expand(rules.break_fasta_apart.output[0]+"/{chunk}.fasta", 
        database=wildcards.database, 
        mapping_db=wildcards.mapping_db,
        chunk=chunks)
    return(out)


def get_fasta_chunk(wildcards):
    chunk_dir=checkpoints.break_fasta_apart.get(
        database=wildcards.database,
        mapping_db=wildcards.mapping_db).output[0]
    return str(Path(chunk_dir) / (wildcards.chunk + ".fasta"))


def get_emapper_chunks(wildcards):
    return expand(rules.emapper.output,
        database=wildcards.database,
        mapping_db=wildcards.mapping_db,
        chunk=get_chunks(wildcards))


def get_chunk_headers(wildcards):
    return expand(rules.get_headers.output.og_header,
        database=wildcards.database,
        mapping_db=wildcards.mapping_db,
        chunk=get_chunks(wildcards))


rule emapper:
    input:
        fasta=get_fasta_chunk,
        diamond=rules.download_diamond.output,
        annotations=rules.download_eggnog_annotations.output,
        taxa=rules.download_eggnog_taxa.output
    output:"results/{database}/annotation/emapper/{mapping_db}/{chunk}/go_annotations.emapper.annotations"
    params:
        prefix="go_annotations",
        outdir=lambda wildcards, output: str(Path(output[0]).parent),
        go_evidence=config["emapper"]["go_evidence"],
        data_dir="resources/eggnog"
    conda: "../envs/eggnog.yml"
    log: "logs/{database}/annotation/emapper/{mapping_db}/{chunk}.log"
    threads:config["emapper"]["threads"]
    shell:
        """
        mkdir -p {params.outdir:q}
        emapper.py -i {input.fasta:q} \
           -o {params.prefix:q} \
           --output_dir {params.outdir:q} \
           --cpu {threads} \
           --go_evidence {params.go_evidence} \
           --data_dir {params.data_dir:q} \
           --override > {log:q} 2>&1
        """

rule combine_emapper:
    input:
        annotations=get_emapper_chunks,
        fastas=get_fasta,
        linkers=get_chunk_headers
    output:"results/{database}/annotation/combine_emapper/{mapping_db}/go_annotations.emapper.annotations"
    script:"../scripts/combine_emapper_annotations.py"

rule make_go_db:
    input:
        annotations=rules.combine_emapper.output,
        script="workflow/scripts/make_go_db.py",
        helpers="workflow/scripts/annotation_io.py"
    output:"results/{database}/annotation/make_go_db/{mapping_db}.3"
    shell:
        """
        python {input.script:q} {input.annotations:q} {output:q}
        """

rule parse_blast:
    input:
        db=rules.make_go_db.output,
        blast=lambda wildcards: get_annotation_input(wildcards,
            "results/{database}/clustering/mmseqs2/convertalis_full_blast/{target_db}_{mapping_db}_convertlis.8")
    output:"results/{database}/annotation/parse_blast/{target_db}_{mapping_db}.tsv"
    conda: "../envs/eggnog.yml"
    shell:
        """
        python workflow/scripts/parse_blast.py \
	    {input.blast:q} \
	    {input.db:q} \
	    {output:q} \
	    2
        """

rule reannotate_blast_results:
    input:
        processed_go=rules.parse_blast.output,
        phenotype=rules.download_go_ontology.output,
        script="workflow/scripts/annotate_go_terms.py"
    output:"results/{database}/annotation/reannotate_blast_results/{target_db}_{mapping_db}.tsv"
    conda: "../envs/eggnog.yml"
    shell:
        """
        python {input.script:q} --obo {input.phenotype:q} --go_ids {input.processed_go:q} --out {output:q}
        """

rule export_go_annotations:
    input:
        annotations=rules.combine_emapper.output,
        hits=lambda wildcards: get_annotation_input(wildcards,
            "results/{database}/binary/combined_species_hits/{mapping_db}_{database}.tsv"),
        ontology=rules.download_go_ontology.output,
        script="workflow/scripts/export_go_annotations.py",
        metadata="workflow/scripts/annotate_go_terms.py",
        helpers="workflow/scripts/annotation_io.py"
    output:
        families="results/{database}/annotation/export_go_annotations/{mapping_db}/family_annotations.tsv",
        evidence="results/{database}/annotation/export_go_annotations/{mapping_db}/evidence.tsv"
    params:
        go_evidence=config["emapper"]["go_evidence"]
    conda: "../envs/eggnog.yml"
    shell:
        """
        python {input.script:q} --annotations {input.annotations:q} --hits {input.hits:q} \
            --obo {input.ontology:q} --families {output.families:q} --evidence {output.evidence:q} \
            --go-evidence {params.go_evidence:q}
        """

rule link_go_nodes:
    input:
        annotations=rules.export_go_annotations.output.families,
        hits=lambda wildcards: get_annotation_input(wildcards,
            "results/{database}/binary/combined_species_hits/{mapping_db}_{database}.tsv"),
        pairs=lambda wildcards: get_annotation_input(wildcards,
            "results/{database}/binary/get_taxonomy/{mapping_db}-binary-temp.csv"),
        script="workflow/scripts/sequence_to_annotation_linker.py"
    output:"results/{database}/annotation/link_go_nodes/{mapping_db}/node_annotations.tsv"
    conda:"../envs/annotation.yml"
    shell:
        """
        python {input.script:q} {input.hits:q} {input.annotations:q} {input.pairs:q} {output:q}
        """

rule get_fasta_seqs_to_annotate:
    input: 
        mapped=lambda wildcards: get_annotation_input(wildcards,
            "results/{database}/binary/combined_species_hits/{mapping_db}_{database}.tsv"),
        clustered=lambda wildcards: get_annotation_input(wildcards,
            "results/{database}/clustering/mmseqs2_linclust/{mapping_db}/unaligned_linclust_rep_seq.fasta"),
        references=list(config["target_db"]["paths"].values()),
        script="workflow/scripts/select_annotation_sequences.py",
        helpers="workflow/scripts/annotation_io.py"
    output: 
        fasta=temporary("results/{database}/annotation/get_fasta_seqs_to_annotate/{mapping_db}.fasta"),
        mapping=temporary("results/{database}/annotation/get_fasta_seqs_to_annotate/{mapping_db}.txt")
    conda: "../envs/clustering.yml"
    shell:
        """
        python {input.script:q} --hits {input.mapped:q} --references {input.references:q} \
            --denovo {input.clustered:q} --output {output.fasta:q} --identifiers {output.mapping:q}
        """

rule anvio_reformat_fasta:
    input:rules.get_fasta_seqs_to_annotate.output.fasta
    output:
        fasta="results/{database}/annotation/anvio_reformat_fasta/{mapping_db}.fasta",
        report="results/{database}/annotation/anvio_reformat_fasta/{mapping_db}.names.tsv"
    conda:"../envs/annotation.yml"
    shell:
        """
        anvi-script-reformat-fasta --simplify-names {input:q} -o {output.fasta:q} --report-file {output.report:q}
        """

checkpoint break_fasta_apart:
    input:rules.anvio_reformat_fasta.output.fasta
    output:directory("results/{database}/annotation/break_fasta_apart/{mapping_db}")
    conda:"../envs/annotation.yml"
    shell:
        """
        seqkit split {input:q} -s 10000 -O {output:q}
        """

rule get_headers:
    input:
        report=rules.anvio_reformat_fasta.output.report,
        fasta=get_fasta_chunk,
        script="workflow/scripts/make_annotation_headers.py",
        helpers="workflow/scripts/annotation_io.py"
    output:
        og_header="results/{database}/annotation/get_headers/{mapping_db}/{chunk}_header.txt",
        temp=temporary("results/{database}/annotation/get_headers/{mapping_db}/{chunk}.tmp")
    conda:"../envs/annotation.yml"
    shell:
        """
        python {input.script:q} --fasta {input.fasta:q} --report {input.report:q} \
            --output {output.og_header:q} --identifiers {output.temp:q}
        """

rule make_gene_calls:
    input:
        fasta=get_fasta_chunk,
        script="workflow/scripts/make_external_gene_calls.sh"
    output:"results/{database}/annotation/make_gene_calls/{mapping_db}/{chunk}.tsv"
    shell:
        """
        bash {input.script:q} {input.fasta:q} {output:q}
        """

rule anvio_contigs_db:
    input:
        gene_calls=rules.make_gene_calls.output,
        faa=get_fasta_chunk
    output:"results/{database}/annotation/anvio_contigs_db/{mapping_db}{chunk}.db"
    conda:"../envs/annotation.yml"
    threads: config["anvio"]["contigs_db"]["threads"]
    shell:
        """
        anvi-gen-contigs-database --allow-amino-acid-contig-db -T {threads} -o {output:q} -f {input.faa:q} --external-gene-calls {input.gene_calls:q}
        """

rule anvio_kegg_annotation:
    input:
        db=rules.anvio_contigs_db.output,
        kofam=rules.anvio_setup_kegg_kofams.output
    output:"results/temp/{database}/annotation/anvio_kegg_annotation/{mapping_db}/{chunk}_anvio_run_kegg_kofams.0"
    conda:"../envs/annotation.yml"
    log: "logs/{database}/annotation/anvio_run_kegg_kofams/{mapping_db}_{chunk}.log"
    threads: config["anvio"]["run_kegg_annotations"]["threads"]
    shell:
        """
        anvi-run-kegg-kofams -c {input.db:q} --kegg-data-dir {input.kofam:q} -T {threads} --just-do-it 2> {log:q}
        touch {output:q}
        """

rule anvio_export_functions:
    input:
        kegg=rules.anvio_kegg_annotation.output,
        db=rules.anvio_contigs_db.output,
    output: "results/{database}/annotation/anvio_export_functions/{mapping_db}_{chunk}.tsv"
    conda:"../envs/annotation.yml"
    shell:
        """
        anvi-export-functions -c {input.db:q} -o {output:q}
        """

rule link_information:
    input:
        linker=rules.get_headers.output.og_header,
        gene_calls=rules.make_gene_calls.output,
        functions=rules.anvio_export_functions.output,
        script="workflow/scripts/merge_annotations.py",
        helpers="workflow/scripts/annotation_io.py"
    output:"results/{database}/annotation/link_information/{mapping_db}/{chunk}_merged.tsv"
    shell:
        """
        python {input.script:q} {input.functions:q} {input.gene_calls:q} {input.linker:q} {output:q}
        """

rule combine_link_info:
    input:
        chunks=get_all_links,
        script="workflow/scripts/combine_ko_annotations.py",
        helpers="workflow/scripts/annotation_io.py",
        columns="workflow/scripts/merge_annotations.py"
    output:"results/{database}/annotation/combine_link_info/{mapping_db}/ko_merged.tsv"
    shell:
        """
        python {input.script:q} --output {output:q} {input.chunks:q}
        """

rule link_nodes:
    input:
        seqs=rules.combine_link_info.output,
        merged=lambda wildcards: get_annotation_input(wildcards,
            "results/{database}/binary/combined_species_hits/{mapping_db}_{database}.tsv"),
        internal_map=lambda wildcards: get_annotation_input(wildcards,
            "results/{database}/binary/get_taxonomy/{mapping_db}-binary-temp.csv"),
        script="workflow/scripts/sequence_to_annotation_linker.py"
    output:"results/{database}/annotation/link_nodes/{mapping_db}/node_annotations.tsv"
    conda:"../envs/annotation.yml"
    shell:
        """
        python {input.script:q} {input.merged:q} {input.seqs:q} {input.internal_map:q} {output:q}
        """

rule run_annotation:
    input:
        ko=rules.link_nodes.output,
        go=rules.link_go_nodes.output
    output:"results/{database}/annotation/run_anntoation/{mapping_db}.done"
    shell:
        """
        touch {output:q}
        """
