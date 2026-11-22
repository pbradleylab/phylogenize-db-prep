#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(ape)
  library(castor)
  library(optparse)
  library(readr)
})

opt_list <- list(
  make_option(c("-i", "--input"), type="character", help="path to initial tree file"),
  make_option(c("-a", "--altinput"), type="character", default=NA, help="optional additional tree, for example archaeal genomes"),
  make_option("--alternate_input", type="character", default=NA, help="alias for --altinput"),
  make_option(c("-n", "--altname"), type="character", default=NA, help="legacy option; subtree names are taken from taxonomy phyla"),
  make_option(c("-t", "--taxonomy"), type="character", help="path to taxonomy .csv table"),
  make_option(c("-o", "--output"), type="character", help="path to output rds file"),
  make_option("--allow_missing_species", type="character", default=NA, help="CSV listing exactly the audited source species without tree tips, in a cluster column")
)
prs <- OptionParser(option_list = opt_list)
p <- parse_args(prs)

for (option in c("input", "taxonomy", "output")) {
  if (is.null(p[[option]]) || is.na(p[[option]]) || !nzchar(p[[option]])) {
    stop(paste0("Missing required --", option))
  }
}

tax <- readr::read_csv(p$taxonomy, col_types=cols(.default=col_character()), na=character(), show_col_types=FALSE)
if (!all(c("cluster", "phylum") %in% names(tax))) {
  stop("Taxonomy must contain cluster and phylum columns")
}
species_keys <- sub("\\.[0-9]+$", "", tax$cluster)
if (!nrow(tax) || anyNA(tax[c("cluster", "phylum")]) ||
    any(trimws(tax$cluster) == "") || any(trimws(tax$phylum) == "") || anyDuplicated(species_keys)) {
  stop("Taxonomy has missing or ambiguous species IDs or phyla")
}

alternate <- unique(c(p$altinput, p$alternate_input))
alternate <- alternate[!is.na(alternate) & alternate != ""]
if (length(alternate) > 1) {
  stop("Supply only one alternate tree using --altinput or --alternate_input")
}
tree_paths <- c(p$input, alternate)
trees <- lapply(tree_paths, function(tree_path) {
  tree <- ape::read.tree(tree_path)
  if (!inherits(tree, "phylo")) {
    stop(paste0("Input must contain exactly one phylogenetic tree: ", tree_path))
  }
  if (is.null(tree$edge.length) || any(!is.finite(tree$edge.length)) || any(tree$edge.length < 0)) {
    stop(paste0("Tree must have finite, nonnegative branch lengths: ", tree_path))
  }
  # Normalize matching keys only; every output tip retains its canonical taxonomy ID.
  matched <- match(sub("\\.[0-9]+$", "", tree$tip.label), species_keys)
  tree$tip.label[!is.na(matched)] <- tax$cluster[matched[!is.na(matched)]]
  if (anyNA(tree$tip.label) || any(tree$tip.label == "") || anyDuplicated(tree$tip.label)) {
    stop(paste0("Tree tips are missing or ambiguous after species-ID normalization: ", tree_path))
  }
  tree
})

covered <- unlist(lapply(trees, function(tree) intersect(tree$tip.label, tax$cluster)), use.names=FALSE)
if (anyDuplicated(covered)) {
  stop("Species occur in more than one input tree")
}
missing_species <- setdiff(tax$cluster, covered)
allowed_missing <- character()
if (!is.na(p$allow_missing_species)) {
  exclusions <- readr::read_csv(p$allow_missing_species, col_types=cols(.default=col_character()), na=character(), show_col_types=FALSE)
  if (!"cluster" %in% names(exclusions)) {
    stop("Audited source-tree exclusions must have a cluster column")
  }
  allowed_missing <- exclusions$cluster
}
if (anyNA(allowed_missing) || any(trimws(allowed_missing) == "") || anyDuplicated(allowed_missing) ||
    !setequal(missing_species, allowed_missing)) {
  stop(paste0("Unexpected species missing from input trees; supply an alternate tree or an exact audited --allow_missing_species list. Missing species: ", paste(head(missing_species, 20), collapse=", ")))
}
if (length(allowed_missing)) {
  warning(paste0(length(allowed_missing), " audited source species have no tree tips and remain unavailable for phylogenetic analyses."))
}

# Normalize each complete source tree before selecting per-phylum subtrees.
red_trees <- lapply(seq_along(trees), function(index) {
  red_tree <- castor::date_tree_red(castor::root_at_midpoint(trees[[index]]))
  if (!isTRUE(red_tree$success)) {
    stop(paste0("Couldn't date tree using RED: ", tree_paths[index]))
  }
  red_tree$tree
})
phyla <- unique(tax$phylum)
red_subtrees <- setNames(lapply(phyla, function(phylum) {
  species <- setdiff(tax$cluster[tax$phylum == phylum], allowed_missing)
  if (!length(species)) return(NULL)
  candidates <- Filter(function(tree) all(species %in% tree$tip.label), red_trees)
  if (length(candidates) != 1) {
    stop(paste0("Phylum must be covered by exactly one input tree: ", phylum))
  }
  ape::keep.tip(candidates[[1]], species)
}), phyla)
red_subtrees <- Filter(Negate(is.null), red_subtrees)

dir.create(dirname(p$output), recursive=TRUE, showWarnings=FALSE)
saveRDS(red_subtrees, p$output)
