#!/usr/bin/env Rscript
# Based on minimap2 output, create a .csv in gene_presence_absence.csv format (with at least the important columns)
library(optparse)
library(readr)
library(dplyr)
library(tibble)
library(tidyr)

opt_list <- list(
  make_option(c("-i", "--input"), type="character", help="path to tabular minimap2 output file"),
  make_option(c("-o", "--output"), type="character", help="path to output file")
)
prs <- OptionParser(option_list = opt_list)
p <- parse_args(prs)

mm <- read_tsv(p$input, col_names=FALSE, col_types="cdddccddddddc")
wide <- mm %>%
  group_by(X6) %>%
  slice_max(X12) %>%
  separate_wider_delim(X1, delim='_', cols_remove=FALSE, names=c("genome","gene")) %>%
  select(X6,  genome, X1) %>%
  pivot_wider(names_from=genome, values_from=X1, values_fn=\(x) paste(x, collapse=';'), values_fill="") %>%
  mutate(`Non-unique Gene name` = "", Annotation = "") %>% # dummy columns to match old .csv
  relocate(Gene, `Non-unique Gene name`, Annotation)
write_csv(wide, p$output)
