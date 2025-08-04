#!/usr/bin/env Rscript
library(ape)
library(castor)
library(optparse)
library(readr)
library(dplyr)

opt_list <- list(
  make_option(c("-i", "--input"), type="character", help="path to initial tree file"),
  make_option(c("-a", "--altinput"), type="character", default=NULL, help="(optional) path to secondary file"),
  make_option(c("-t", "--taxonomy"), type="character", help="path to taxonomy .csv table"),
  make_option(c("-o", "--output"), type="character", help="path to output rds file")
)
prs <- OptionParser(option_list = opt_list)
p <- parse_args(prs)

tree <- ape::read.tree(p$input)
tax <- readr::read_csv(p$taxonomy)

red_tree <- castor::date_tree_red(castor::root_at_midpoint(tree))
if (!red_tree$success) { stop("Error: couldn't date tree using RED") }

phyla <- dplyr::select(tax, "phylum") %>% dplyr::distinct() %>% dplyr::pull("phylum")
red_subtrees <- purrr::map(phyla, ~ {
  species <- dplyr::filter(tax, phylum %in% .x) %>%
    dplyr::pull("cluster")
  ape::keep.tip(red_tree$tree,
    intersect(red_tree$tree$tip.label, species))
}) %>%
  setNames(phyla) %>%
  (\(.) .[which(!(purrr::map_lgl(., is.null)))])


if (!is.null(p$altinput)) {
  alt_tree <- castor::date_tree_red(castor::root_at_midpoint(tree))
  red_subtrees_alt <- purrr::map(phyla, ~ {
    species <- dplyr::filter(tax, phylum %in% .x) %>%
      dplyr::pull("cluster")
    ape::keep.tip(alt_tree$tree,
      intersect(alt_tree$tree$tip.label, species))
  }) %>%
    setNames(phyla) %>%
    (\(.) .[which(!(purrr::map_lgl(., is.null)))])
  
  for (name in names(red_subtrees_alt)) {
    red_subtrees[[name]] <- red_subtrees_alt[[name]]
  }
  message("Secondary input loaded.")
} else {
  message("No secondary input provided.")
}


saveRDS(red_subtrees, p$output)
