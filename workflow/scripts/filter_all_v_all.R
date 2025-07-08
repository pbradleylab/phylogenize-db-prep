#!/usr/bin/env Rscript
library(optparse)
library(readr)
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)

opt_list <- list(
  make_option(c("-i", "--input"), type="character", help="path to all-versus-all results (in blast6out format)"),
  make_option(c("-f", "--fasta"), type="character", help="path to FASTA file of SSU rRNAs"),
  make_option(c("-m", "--metadata"), type="character", help="path to TSV containing genome metadata"),
  make_option(c("-o", "--output"), type="character", help="path to output file")
)
p <- OptionParser(option_list = opt_list)
parse_args(p)

md <- read_tsv(p$metadata) %>%
  separate_wider_delim(Lineage,
    delim=";",
    names=c(
      "domain", "phylum", "class",
      "order", "family", "genus",
      "species"))

tax <- select(md, Species_rep, domain:species) %>% distinct()

all_v_all <- read_tsv(p$input, col_names=FALSE) %>%
  mutate(q_species_id = gsub("(.*);;(.*);;(.*)", "\\3", X1)) %>%
  mutate(d_species_id = gsub("(.*);;(.*);;(.*)", "\\3", X2)) %>%
  relocate(q_species_id, d_species_id)

ava_tax <- left_join(all_v_all, tax, by=c("q_species_id"="Species_rep")) %>%
  left_join(., tax, by=c("d_species_id"="Species_rep"), suffix=c("_q", "_d"))

# in progress -
# still need to adapt from notes
