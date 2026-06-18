

rule phylogenize:
    input:
        abundance_file=config["phylogenize_abundance_file"],
        metadata_file=config["metadata_file"]
    output:"results/{project}/phylogenetics/phylogenize/"+config["phylogenize"]["which_phenotype"]+".html"
    params:
        outdir=config["phylogenize"]["out_dir"],
        db=config["phylogenize"]["db"],
        taxon_level=config["phylogenize"]["taxon_level"],
        type_16S=config["phylogenize"]["type_16S"],
        which_phenotype=config["phylogenize"]["which_phenotype"],
        input_format=config["phylogenize"]["input_format"],
        sample_column=config["phylogenize"]["sample_column"],
        vsearch_bin=config["phylogenize"]["vsearch_bin"],
        abundance_method=config["phylogenize"]["diff_abund_method"]
    log:"results/{project}/phylogenetics/phylogenize.log"
    conda:"../envs/phylogenetics.yml"
    threads:config["phylogenize"]["threads"]
    shell:
        """
        Rscript workflows/scripts/phylogenize_setup.R 2> {log
        """
