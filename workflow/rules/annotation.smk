from glob import glob
import os, peppy

include: "resources.smk"
include: "matrix.smk"


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
    
    chunks=[p.stem for p in Path(chunk_dir).glob("*.fasta")]
    return(chunks)
    
def get_fasta(wildcards):
    chunks=get_chunks(wildcards)
    out=expand(rules.break_fasta_apart.output[0]+"/{chunk}.fasta", 
        database=wildcards.database, 
        mapping_db=wildcards.mapping_db,
        chunk=chunks)
    return(out)

def get_ko(wildcards):
    chunk_dir=checkpoints.break_fasta_apart.get(
        database=wildcards.database,
        mapping_db=wildcards.mapping_db).output[0]
    
    chunks=[p.stem for p in Path(chunk_dir).glob("*.fasta")]
    out=expand(rules.link_nodes.output, 
        database=wildcards.database,
        mapping_db=wildcards.mapping_db,
        chunk=chunks)
    return(out)

rule emapper:
    input:
        targets=get_targets,
        diamond=rules.download_diamond.output
    output:"results/{database}/annotation/emapper/{target_db}/go_annotations.emapper.annotations"
    params:
        prefix="go_annotations",
        go_evidence=config["emapper"]["go_evidence"]
    conda: "../envs/annotation.yml"
    threads:config["emapper"]["threads"]
    shell:
        """
        outdir=$(dirname {output})
        emapper.py -i {input.targets} \
           -o $outdir/{params.prefix} \
           --cpu {threads} \
           --go_evidence {params.go_evidence} \
           --data_dir {input.diamond}
        """

rule make_go_db:
    input:rules.emapper.output
    output:"results/{database}/annotation/make_go_db/{target_db}.3"
    shell:
        """
        cut -f1,10 {input} | grep -v '##' | grep -v "-" > {output}
        """

rule parse_blast:
    input:
        db=rules.make_go_db.output,
        blast=rules.mmseqs2_convertalis_full_blast.output
    output:"results/{database}/annotation/parse_blast/{target_db}_{mapping_db}.tsv"
    conda: "../envs/annotation.yml"
    shell:
        """
        python workflow/scripts/parse_blast.py \
	    {input.blast} \
	    {input.db} \
	    {output} \
	    2
        """

rule reannotate_blast_results:
    input:
        processed_go=rules.parse_blast.output
        phenotype=rule.download_go_ontology.output
    output:"results/{database}/annotation/reannotate_blast_results/{target_db}_{mapping_db}.tsv"
    conda: "../envs/annotation.yml"
    shell:
        """
        python workflow/scripts/annotate_go_terms.py --obo {input.phenotype} --go_ids {input.processed_go} --out {output} 
        """

rule get_fasta_seqs_to_annotate:
    input: 
        mapped=rules.combine_species_hits.output,
        clustered=rules.mmseqs2_linclust.output.fasta
    output: 
        fasta=temporary("results/{database}/annotation/get_fasta_seqs_to_annotate/{mapping_db}.fasta"),
        mapping=temporary("results/{database}/annotation/get_fasta_seqs_to_annotate/{mapping_db}.txt")
    params:
        database=list(config["target_db"]["paths"].values())
    conda: "../envs/clustering.yml"
    shell:
        """
        cut -f2 {input.mapped} > {output.mapping}
        for f in {params.database}
        do
            faSomeRecords $f {output.mapping} stdout >> {output.fasta}
        done
        cat {input.clustered} >> {output.fasta}
        """

rule anvio_reformat_fasta:
    input:rules.get_fasta_seqs_to_annotate.output.fasta
    output:"results/{database}/annotation/anvio_reformat_fasta/{mapping_db}.fasta"
    conda:"../envs/annotation.yml"
    shell:
        """
        anvi-script-reformat-fasta --simplify-names {input} -o {output} 
        """

# remember to change get_header if -s is changed in this checkpoint
checkpoint break_fasta_apart:
    input:rules.anvio_reformat_fasta.output
    output:directory("results/{database}/annotation/break_fasta_apart/{mapping_db}")
    conda:"../envs/annotation.yml"
    shell:
        """
        seqkit split {input} -s 10000 -O {output}
        """

rule get_headers:
    input:
        clustered=rules.get_fasta_seqs_to_annotate.output.fasta,
        fasta=get_fasta
    output:
        og_header="results/{database}/annotation/get_headers/{mapping_db}/{chunk}_header.txt",
        temp=temporary("results/{database}/annotation/get_headers/{mapping_db}/{chunk}.tmp")
    conda:"../envs/annotation.yml"
    shell:
        """
        echo "gene_callers_id\tlinker_info" > {output.og_header} 
        in=$(echo {input.fasta} | tr ' ' '\\n' | grep '{wildcards.chunk}\\.fasta')
        grep ">" "$in" | sed 's/>c_//g' | sed 's/^0*//' > {output.temp}
        
        start=$(head -n1 {output.temp})
        end=$(tail -n1 {output.temp})
        
        seqkit range -r $start:$end {input.clustered} | grep '>' | sed 's/>//g' | sed 's/ .*//g' > {output.temp}
        awk -F, '{{print FNR-1,$0}}' OFS='\t' {output.temp} >> {output.og_header}
        """

rule make_gene_calls:
    input: get_fasta
    output:"results/{database}/annotation/make_gene_calls/{mapping_db}/{chunk}.tsv"
    shell:
        """
        in=$(echo {input} | tr ' ' '\\n' | grep '{wildcards.chunk}\\.fasta')
        bash workflow/scripts/make_external_gene_calls.sh $in {output}
        """

rule anvio_contigs_db:
    input:
        gene_calls=rules.make_gene_calls.output,
        faa=get_fasta
    output:"results/{database}/annotation/anvio_contigs_db/{mapping_db}{chunk}.db"
    conda:"../envs/annotation.yml"
    threads: config["anvio"]["contigs_db"]["threads"]
    shell:
        """
        in=$(echo {input.faa} | tr ' ' '\n' | grep '{wildcards.chunk}.fasta')
        anvi-gen-contigs-database --allow-amino-acid-contig-db -T {threads} -o {output} -f $in --external-gene-calls {input.gene_calls}
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

rule link_information:
    input:
        linker=rules.get_headers.output.og_header,
        gene_calls=rules.make_gene_calls.output,
        functions=rules.anvio_export_functions.output
    output:"results/{database}/annotation/link_information/{mapping_db}/{chunk}_merged.tsv"
    shell:
        """
        python workflow/scripts/merge_annotations.py {input.functions} {input.gene_calls} {input.linker} {output}
        """

rule combine_link_info:
    input:get_all_links
    output:"results/{database}/annotation/combine_link_info/{mapping_db}/ko_merged.tsv"
    shell:
        """
        head -n 1 {input}[0] > {output}

        for f in {input}; do
            tail -n +2 "$f" >> {output}
        done
        """

rule link_nodes:
    input:
        seqs=rules.combine_link_info.output,
        merged=rules.combine_species_hits.output,
        internal_map=rules.get_taxonomy.output.out
    output:"results/{database}/annotation/link_nodes/{mapping_db}/node_annotations.tsv"
    conda:"../envs/annotation.yml"
    params:
    shell:
        """
        python workflow/scripts/sequence_to_annotation_linker.py {input.merged} {input.seqs} {input.internal_map} {output}
        """

rule run_annotation:
    input:get_ko
    output:"results/{database}/annotation/run_anntoation/{mapping_db}.done"
    shell:
        """
        touch {output}
        """
