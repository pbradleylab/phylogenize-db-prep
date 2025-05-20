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
