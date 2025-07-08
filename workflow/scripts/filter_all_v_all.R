#!/usr/bin/env Rscript
library(optparse)
library(seqinr)
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
prs <- OptionParser(option_list = opt_list)
prs_list <- parse_args(prs)
p <- prs_list$options

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

very_unlikely <- ava_tax %>%
  group_by(X1) %>%
  count(agree=ifelse(family_d==family_q, 'Y', 'N')) %>%
  pivot_wider(values_from=n, names_from=agree, values_fill=0) %>%
  mutate(T = Y + N) %>% filter((Y / T) < 0.5, T > 2) %>%
  pull(X1)

pare_by <- function(res, exclude=NULL, lvl="family", frac=0.9, min_n=2) {
  level_d <- paste0(lvl, "_d")
  level_q <- paste0(lvl, "_q")
  res2 <- res %>%
    filter(!(X1 %in% exclude)) %>%
    filter(!(X2 %in% exclude)) %>%
    group_by(X1) %>%
    count(agree=ifelse(
      !!(rlang::sym(level_d)) == !!(rlang::sym(level_q)),
      'Y',
      'N')
    ) %>%
    pivot_wider(values_from=n, names_from=agree, values_fill=0)
  new_exc <- res2 %>%
    mutate(T = Y + N) %>%
    filter((Y / T) <= frac, T > min_n) %>%
    pull(X1)
  c(exclude, new_exc)
}

to_exclude <- pare_by(ava_tax, lvl="family", frac=0.8, min_n=2)
prev_excluded <- c()
round <- 0
message(paste0("Round ", round, ": ", length(to_exclude), " sequences excluded"))
while(length(setdiff(to_exclude, prev_excluded)) > 1) {
  round <- round + 1
  prev_excluded <- to_exclude
  to_exclude <- pare_by(ava_tax, exclude=prev_excluded, lvl="family", frac=0.8, min_n=5)
  to_exclude <- pare_by(ava_tax, exclude=to_exclude, lvl="genus", frac=0.5, min_n=5)
  message(paste0("Round ", round, ": ", length(setdiff(to_exclude, prev_excluded)), " sequences excluded"))
}
message("Stopping...")
filtered <- ava_tax %>% filter(!(X1 %in% already_excluded), !(X2 %in% already_excluded))

message(paste0("Removed a total of ", length(already_excluded), " out of ", nrow(ava_tax), " sequences"))

fa <- read.fasta(p$fasta, seqtype="DNA", as.string=TRUE, forceDNAtolower = FALSE)
fa <- fa[intersect(filtered$X1, names(fa))]
write.fasta(fa, names(fa), file.out = p$output)
