""" Create the input databases for UHGP50 and UniRef50. UniRef50 will be used to map and cluster
the database and then the sequences that are not found at a 50% sequence identity are mapped
against UHGP50.

By Kathryn Kananen
"""


rule phylogenize_setup:
    output:"results/temp/{database}/phylogenize_setup.done"
    log: "logs/{database}/resources/phylogenize_setup.log"
    conda:"../envs/phylogenetics.yml"
    shell:
        """
        Rscript -e "phylogenize::install.data.figshare()" 2> {log}
        """

rule anvio_setup_kegg_kofams:
    output:directory("resources/anvio/kofams/")
    conda:"../envs/annotation.yml"
    shell:
        """
        anvi-setup-kegg-data --kegg-data-dir {output}
        """

rule download_diamond:
    output:"resources/eggnog/eggnog_proteins.dmnd"
    params:
        data_dir="resources/eggnog"
    conda:"../envs/eggnog.yml"
    shell:
        """
        mkdir -p {params.data_dir:q}
        python workflow/scripts/download_eggnog_data.py --db diamond -y --data_dir {params.data_dir:q}
        test -s {output:q}
        """

rule download_eggnog_annotations:
    output:"resources/eggnog/eggnog.db"
    params:
        data_dir="resources/eggnog"
    conda:"../envs/eggnog.yml"
    shell:
        """
        mkdir -p {params.data_dir:q}
        python workflow/scripts/download_eggnog_data.py --db annotations -y --data_dir {params.data_dir:q}
        test -s {output:q}
        """

rule download_eggnog_taxa:
    output:
        db="resources/eggnog/eggnog.taxa.db",
        traversal="resources/eggnog/eggnog.taxa.db.traverse.pkl"
    params:
        data_dir="resources/eggnog"
    conda:"../envs/eggnog.yml"
    shell:
        """
        mkdir -p {params.data_dir:q}
        python workflow/scripts/download_eggnog_data.py --db taxa -y --data_dir {params.data_dir:q}
        for db in {output:q}; do
            test -s "$db"
        done
        """

rule download_go_ontology:
    output:"resources/go_ontology/obo/go.obo"
    conda:"../envs/eggnog.yml"
    shell:
        """
        mkdir $(dirname {output})
        wget -O http://purl.obolibrary.org/obo/go.obo {output}
        """
