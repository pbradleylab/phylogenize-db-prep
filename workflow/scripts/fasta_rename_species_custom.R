#!/usr/bin/env Rscript
# Rename species in a FASTA file where the names are "<gene_id>;;<genus>;;<species_id>"
library(optparse)
library(seqinr)
library(purrr)
library(tibble)
library(dplyr)
library(readr)
library(stringr)

opt_list <- list(
  make_option(c("-i", "--input"), type="character", help="path to FASTA file to rename"),
  make_option(c("-m", "--metadata"), type="character", help="path to metadata file"),
  make_option(c("-g", "--genome_column"), type="character", default="MGnify_accession", help="which column in metadata to match on"),
  make_option(c("-d", "--id_column"), type="character", default="species_id", help="which column to grab new IDs from"),
  make_option(c("-r", "--regex_match"), type="character", default="MGYG-HGUT-", help="regex to apply to genome-column to make it match fasta file"),
  make_option(c("-R", "--regex_sub"), type="character", default="MGYG0000", help="replacement text for regex"),
  make_option(c("-o", "--output"), type="character", help="path to output file")
)
prs <- OptionParser(option_list = opt_list)
p <- parse_args(prs)

# read data and split off species IDs
fa <- read.fasta(p$input, seqtype="DNA", as.string=TRUE, forceDNAtolower = FALSE)
md <- read_tsv(p$metadata)
species_ids <- map_chr(names(fa), ~ str_split_1(string=.x, pattern=";;")[3])

# apply regex 
md_t <- md
md_t[, p$genome_column] <- gsub(p$regex_match, p$regex_sub, md_t[, p$genome_column])

# make mapping
sid_vec <- enframe(species_ids, name="fullname", value=p$genome_column) %>%
  left_join(md_t) %>%
  select(fullname, all_of(p$id_column)) %>%
  deframe()
print(head(sid_vec))

# map IDs
new_species_ids <- map_chr(names(fa), ~ {
  sp <- str_split_1(string=.x, pattern=";;")
  if (!(.x %in% names(sid_vec))) { warning(paste0("Not found: ", .x))}
  paste(c(sp[1:2], sid_vec[.x]), collapse=';;')
})

# output
write.fasta(fa, new_species_ids, file.out = p$output)
