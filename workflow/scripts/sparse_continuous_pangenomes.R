#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(duckplyr)
  library(Matrix)
  library(optparse)
})

opt_list <- list(
  make_option(c("-i", "--input"), type="character", help="path to protein family clustering .tsv file"),
  make_option(c("-g", "--genome_metadata"), type="character", help="path to genome metadata .tsv file"),
  make_option(c("-c", "--combined_hits"), type="character", default=NA, help="path to combined_species_hits .tsv file"),
  make_option(c("-r", "--output_rds"), default="output.rds", type="character", help="path to output .rds file containing sparse per-phylum matrices"),
  make_option(c("-C", "--output_csv"), default="output.csv", type="character", help="path to output .csv file containing all data"),
  make_option(c("-m", "--memory"), default=16, type="numeric",
    help="max memory (in Gb) for DuckDB"),
  make_option(c("-s", "--speedy"), default=FALSE, type="logical", help="Make sparse matrices in memory (uses more memory but MUCH faster)"),
  make_option("--split_char", default="_", type="character", help="delimiter between a member genome ID and its protein ID"),
  make_option("--threads", default=1, type="integer", help="number of DuckDB and CSV reader threads")
)

prs <- OptionParser(option_list = opt_list)
p <- parse_args(prs)

if (is.null(p$input) || is.null(p$genome_metadata)) {
  stop("Both --input and --genome_metadata are required")
}
if (!is.finite(p$memory) || p$memory <= 0 || is.na(p$threads) || p$threads < 1 || p$split_char == "") {
  stop("Memory and threads must be positive and the protein delimiter must not be empty")
}

# Escape the configured delimiter before using it in a regular expression.
delimiter_chars <- strsplit(p$split_char, "", fixed=TRUE)[[1]]
regex_chars <- c("\\", ".", "^", "$", "|", "?", "*", "+", "(", ")", "[", "]", "{", "}")
member_pattern <- paste0(paste0(ifelse(delimiter_chars %in% regex_chars, paste0("\\", delimiter_chars), delimiter_chars), collapse=""), ".*$")

# cap memory usage to avoid getting oom-killed on a shared machine
duckplyr::db_exec(paste0("PRAGMA memory_limit = '", p$memory,"GB'"))
duckplyr::db_exec(paste0("PRAGMA threads = ", p$threads))

# read data
prot_cat <- duckplyr::read_csv_duckdb(p$input, options = list(delim="\t", header=FALSE, all_varchar=TRUE))
genome_md <- duckplyr::read_csv_duckdb(p$genome_metadata, options = list(delim="\t", all_varchar=TRUE))
source_genomes <- genome_md |> select(Genome, Species_rep, Lineage) |> collect()
source_genomes$Genome <- sub("\\.[0-9]+$", "", source_genomes$Genome)
if (nrow(source_genomes) == 0 || anyNA(source_genomes) || any(source_genomes$Genome == "") || any(source_genomes$Species_rep == "") || anyDuplicated(source_genomes$Genome)) {
  stop("Genome metadata contains missing values or ambiguous IDs after version normalization")
}
expected_species <- unique(source_genomes$Species_rep)
genome_md <- genome_md |> mutate(Genome = gsub("\\.[0-9]+$", "", Genome))
prot_cat_split <- prot_cat |>
  mutate(Genome = gsub("\\.[0-9]+$", "", gsub(member_pattern, "", column1)))
if (!is.na(p$combined_hits)) {
  combined_hits <- duckplyr::read_csv_duckdb(p$combined_hits, options = list(delim="\t", all_varchar=TRUE)) |>
    select(query, target) |>
    distinct()
  invalid_hits <- combined_hits |>
    filter(is.na(query) | is.na(target) | query == "" | target == "") |>
    collect()
  conflicting_hits <- combined_hits |>
    summarize(.by=query, families=n()) |>
    filter(families > 1) |>
    collect()
  if (nrow(invalid_hits) || nrow(conflicting_hits)) {
    stop("Combined hits contain missing IDs or a protein centroid maps to more than one target family")
  }
  prot_cat_split <- prot_cat_split |>
    left_join(combined_hits, by=c("column0"="query"))
  missing <- prot_cat_split |> select(column0, target) |> filter(is.na(target)) |> distinct()
  missing_gfs <- nrow(collect(missing))
  if (missing_gfs > 0) {
    stop(paste0("Protein centroids missing from combined hits: ", missing_gfs))
  }
  prot_cat_split <- select(prot_cat_split, column0=target, Genome, column1)
}
genome_tax <- genome_md %>% select(Genome, Lineage, cluster=Species_rep) %>% separate_wider_delim(Lineage, names=c("domain","phylum","class","order","family","genus","species"), delim=";")
genome_tax <- genome_tax |> mutate(phylum=if_else(phylum %in% c("p__", ""), "p__unclassified", phylum))

# use duckplyr to compute frac pangenomes
prot_cat_join <- left_join(prot_cat_split, genome_tax, by="Genome")
unmatched_genomes <- prot_cat_join |>
  filter(is.na(cluster)) |>
  select(Genome) |>
  distinct() |>
  collect()
if (nrow(unmatched_genomes) > 0) {
  stop(paste0("Member genomes missing from metadata: ", paste(head(unmatched_genomes$Genome, 20), collapse=", ")))
}
collapse_genes_by_family_and_genome <- prot_cat_join %>%
  select(-column1) %>%
  distinct()
count_genomes_per_species_per_family <- collapse_genes_by_family_and_genome %>%
  summarize(.by=c(domain:species, cluster, column0),
            genomes_with_family = n())
# The denominator includes every source genome, including those without a family.
count_total_genomes_per_species <- genome_tax %>%
  select(Genome, cluster) %>%
  distinct() %>%
  summarize(.by=c(cluster), total_genomes = n())
fractional_pangenomes <- full_join(count_genomes_per_species_per_family,
                                   count_total_genomes_per_species, by="cluster") %>%
  mutate(frac_observed = genomes_with_family / total_genomes)

observed_species <- fractional_pangenomes |> select(cluster) |> distinct() |> collect() |> pull(cluster)
if (!setequal(expected_species, observed_species) ||
    nrow(collect(fractional_pangenomes |> filter(is.na(column0) | is.na(frac_observed) | frac_observed <= 0 | frac_observed > 1))) > 0) {
  stop("Continuous matrix does not cover every source species")
}

message(paste0("Writing out fractional pangenomes to ", p$output_csv))
dir.create(dirname(p$output_csv), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(p$output_rds), recursive=TRUE, showWarnings=FALSE)
duckplyr::compute_csv(fractional_pangenomes, p$output_csv)

if (p$speedy) {
  message("Re-reading fractional pangenomes...")
  rm(fractional_pangenomes)
  gc()
  fractional_pangenomes <- read_csv(
    p$output_csv,
    col_types=cols(.default=col_character(), genomes_with_family=col_double(), total_genomes=col_double(), frac_observed=col_double()),
    na=character(), num_threads=p$threads, show_col_types=FALSE
  )
}

# compute per-phylum sparse matrices
phyla <- if (p$speedy) {
  sort(unique(fractional_pangenomes$phylum))
} else {
  fractional_pangenomes |> select(phylum) |> distinct() |> collect() |> pull(phylum) |> sort()
}
names(phyla) <- sub("^p__", "", phyla)

message(paste0("Computing per-phylum sparse matrices..."))
per_phylum_matrices <- lapply(phyla, function(phylum_name) {
  message(phylum_name)
  sub_frac <- if (p$speedy) {
    fractional_pangenomes[fractional_pangenomes$phylum == phylum_name, c("column0", "cluster", "frac_observed")]
  } else {
    fractional_pangenomes |>
      filter(phylum==phylum_name) |>
      select(column0, cluster, frac_observed) |>
      collect()
  }
  # Derive indices and labels from the same vectors to avoid deferred row ordering.
  families <- sort(unique(sub_frac$column0))
  species <- sort(unique(sub_frac$cluster))
  mtx <- Matrix::sparseMatrix(
    i=match(sub_frac$column0, families),
    j=match(sub_frac$cluster, species),
    x=sub_frac$frac_observed,
    dims=c(length(families), length(species)),
    dimnames=list(families, species)
  )
  mtx
})

message(paste0("Writing out per-phylum sparse matrices to ", p$output_rds))
saveRDS(per_phylum_matrices, p$output_rds)
