#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(Matrix)
  library(optparse)
})

opt_list <- list(
  make_option(c("-i", "--input"), type="character", help="path to fractional pangenome .csv file"),
  make_option(c("-r", "--output_rds"), default="output.rds", type="character", help="path to output .rds file containing sparse per-phylum matrices")
)

prs <- OptionParser(option_list = opt_list)
p <- parse_args(prs)


fractional_pangenomes <- read_csv(p$input)
# compute per-phylum sparse matrices
phyla <- fractional_pangenomes %>% select(phylum) %>% distinct() %>% pull(phylum)
names(phyla) <- gsub("p__", "", phyla)

message(paste0("Computing per-phylum sparse matrices..."))
per_phylum_matrices <- lapply(phyla, function(p) {
  message(p)
  gc()
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
