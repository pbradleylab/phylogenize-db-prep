def get_fasta_for_pangenome(wildcards):
    indir=config["modes"]["pangenome"]["indir"]
    chunks=[p.stem for p in Path(indir).glob("*.fna")]
    return(chunks)

def get_fasta(wildcards):
    chunks=get_chunks(wildcards)
    out=expand(rules.make_pangenome_contigs_db.output[0]+"/{fasta}.fna",
        database=wildcards.database,
        chunk=chunks)
    return(out)


checkpoint make_pangenome_contigs_db:
    input:get_fasta
    output:"results/{database}/pangenome/anvio_contigs_db/{chunk}.db"
    conda:"../envs/pangenome.yml"
    threads: config["anvio"]["contigs_db"]["threads"]
    shell:
        """
        in=$(echo {input.faa} | tr ' ' '\n' | grep '{wildcards.chunk}.fasta')
        anvi-gen-contigs-database --allow-amino-acid-contig-db -T {threads} -o {output} -f $in --external-gene-calls {input.gene_calls}
        """
