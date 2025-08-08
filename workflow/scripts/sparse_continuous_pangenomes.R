#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(duckplyr)
  library(Matrix)
  library(optparse)
})

opt_list <- list(
  make_option(c("-i", "--input"), type="character", help="path to protein family clustering .tsv file"),
  make_option(c("-g", "--genome_metadata"), type="character", help="path to genome metadata .tsv file"),
  make_option(c("-c", "--combined_hits"), type="character", default=NA, help="path to combined_species_hits .tsv file"),
  make_option(c("-r", "--output_rds"), default="output.rds", type="character", help="path to output .rds file containing sparse per-phylum matrices"),
  make_option(c("-C", "--output_csv"), default="output.csv", type="character", help="path to output .tsv file containing all data"),
  make_option(c("-m", "--memory"), default=16, type="numeric",
    help="max memory (in Gb) for DuckDB")
)

prs <- OptionParser(option_list = opt_list)
p <- parse_args(prs)

# p <- list(
#   input = ~/Data/honeybee-gut/v1.0.1/protein_catalogue/protein_catalogue-50.tsv",
#   genome_metadata = "~/Data/honeybee-gut/v1.0.1/genomes-all_metadata.tsv",
#   output_rds = "test_output.rds",
#   output_csv = "test_output.csv"
# )

# cap memory usage to avoid getting oom-killed on a shared machine
db_exec(paste0("PRAGMA memory_limit = '", p$memory,"GB'"))

# read data
prot_cat <- duckplyr::read_csv_duckdb(p$input, options = list(delim="\t", header=FALSE))
genome_md <- duckplyr::read_csv_duckdb(p$genome_metadata,  options = list(delim="\t"))
prot_cat_split <- prot_cat |> mutate(Genome = gsub("([^_]+)_.*", "\\1", column1))
if (!is.na(p$combined_hits)) {
  combined_hits <- duckplyr::read_csv_duckdb(p$combined_hits, options = list(delim="\t"))
  prot_cat_split <- prot_cat_split |>
    left_join(combined_hits, by=c("column0"="query"))
  missing <- prot_cat_split |> select(column0, target) |> filter(is.na(target)) |> distinct()
  missing_gfs <- nrow(collect(missing))
  if (nrow(collect(missing) > 1)) {
    warning(paste0("Warning: ", (missing_gfs), " gene families were not mapped"))
    print(missing)
  }
  prot_cat_split <- select(prot_cat_split, column0=target, Genome, column1)
}
genome_tax <- genome_md %>% select(Genome, Lineage, cluster=Species_rep) %>% separate_wider_delim(Lineage, names=c("domain","phylum","class","order","family","genus","species"), delim=";")

# use duckplyr to compute frac pangenomes
prot_cat_join <- left_join(prot_cat_split, genome_tax)
collapse_genes_by_family_and_genome <- prot_cat_join %>%
  select(-column1) %>%
  distinct()
count_genomes_per_species_per_family <- collapse_genes_by_family_and_genome %>%
  summarize(.by=c(domain:species, cluster, column0),
            genomes_with_family = n())
count_total_genomes_per_species <- collapse_genes_by_family_and_genome %>%
  select(Genome, cluster) %>%
  distinct() %>%
  summarize(.by=c(cluster), total_genomes = n())
fractional_pangenomes <- full_join(count_genomes_per_species_per_family,
                                   count_total_genomes_per_species) %>%
  mutate(frac_observed = genomes_with_family / total_genomes)

message(paste0("Writing out fractional pangenomes to ", p$output_csv))
duckplyr::compute_csv(fractional_pangenomes, p$output_csv)

# compute per-phylum sparse matrices
phyla <- fractional_pangenomes %>% select(phylum) %>% distinct() %>% pull(phylum)
names(phyla) <- gsub("p__", "", phyla)

message(paste0("Computing per-phylum sparse matrices..."))
per_phylum_matrices <- lapply(phyla, function(p) {
  message(p)
  sub_frac <- fractional_pangenomes |>
    filter(phylum==p)
  gf_nums <- sub_frac |>
    select(column0) |>
    distinct() |>
    mutate(n_genefam = row_number()) |>
    collect()
  cl_nums <- sub_frac |>
    select(cluster) |>
    distinct() |>
    mutate(n_cluster = row_number()) |>
    collect()
  numbered_subfrac <- left_join(sub_frac, gf_nums, by="column0") |>
    left_join(cl_nums, by="cluster")
  sparse_mtx_tbl <- numbered_subfrac |>
    select(n_cluster, n_genefam, frac_observed) |>
    collect()
  mtx <- Matrix::sparseMatrix(
    i=sparse_mtx_tbl[["n_genefam"]],
    j=sparse_mtx_tbl[["n_cluster"]],
    x=sparse_mtx_tbl[["frac_observed"]],
    dimnames=list(gf_nums[["column0"]],
                  cl_nums[["cluster"]])
  )
  mtx
})

message(paste0("Writing out per-phylum sparse matrices to ", p$output_rds))
saveRDS(per_phylum_matrices, p$output_rds)
