library(phylogenize)
library(optparse)

option_list <- list(
  make_option(c("-o", "--output_file"),       type="character",   help="Output file name"),
  make_option(c("-d", "--out_dir"),           type="character",   help="Output directory"),
  make_option(c("-b", "--db"),                type="character",   help="Database to use (e.g. uhgp, gtdb, mouse-gut)"),
  make_option(c("-t", "--taxon_level"),       type="character",   help="Taxonomic level (e.g. phylum)"),
  make_option(c("-s", "--type_16S"),          type="logical",     help="Is this 16S data? TRUE or FALSE"),
  make_option(c("-p", "--which_phenotype"),   type="character",   help="Phenotype type (e.g. abundance, prevalence, specificity)"),
  make_option(c("-a", "--abundance_file"),    type="character",   help="Path to abundance file"),
  make_option(c("-m", "--metadata_file"),     type="character",   help="Path to metadata file"),
  make_option(c("-f", "--input_format"),      type="character",   help="Input format (e.g. tabular)"),
  make_option(c("-e", "--which_envir"),       type="character",   help="Environment label (e.g. RAGKO)"),
  make_option(c("-c", "--sample_column"),     type="character",   help="Sample column name"),
  make_option(c("-v", "--vsearch_bin"),       type="character",   help="Path to vsearch binary"),
  make_option(c("-n", "--ncl"),               type="integer",     help="Number of cores to use"),
  make_option(c("-x", "--diff_abund_method"), type="character", help="Which differential abundance method to use")
)

arg_parser <- OptionParser(option_list=option_list)
args <- parse_args(arg_parser)

phylogenize::render.report(
    output_file=args$output_file,
    out_dir=args$outdir,
    db=args$db,
    taxon_level=args$taxon_level,
    type_16S=args$type_16S,
    which_phenotype=args$which_phenotype,
    diff_abund_method=args$diff_abund_method,
    abundance_file=args$abundance_file,
    metadata_file=args$metadata_file, 
    input_format=args$input_format,
    which_envir=args$which_envir,
    sample_column=args$sample_column,
    vsearch_bin=args$vsearch_bin,
    ncl=args$ncl
)
