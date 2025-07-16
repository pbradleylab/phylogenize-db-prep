#!/usr/bin/env Rscript
library(tidyverse)
library(tidytree)
library(treeio)
library(castor)

opt_list <- list(
  make_option(c("-i", "--input"), type="character", help="path to initial tree file"),
  make_option(c("-d", "--delimiter"), type="character", default=";;",
    help="delimiter for separating gene, genus, and species_id in tip names")
  make_option(c("-o", "--output"), type="character", help="path to filtered tree file")
)
prs <- OptionParser(option_list = opt_list)
p <- parse_args(prs)

tree <- treeio::read.newick(p$input)
tips_by_species <- enframe(tree$tip.label, name="t_num", value="tip") %>%
  separate_wider_delim(tip,
    delim=p$delimiter, names=c("gene", "genus", "species_id"),
    cols_remove=FALSE)
tips_by_genus <- tips_by_species %>% group_by(genus) %>% nest()
message("Calculating mean distances to the same species...")
all_means <- tips_by_species %>%
  group_by(species_id) %>%
  nest() %>%
  mutate(means = map(.progress=TRUE, data, ~ { 
    rowMeans(get_all_pairwise_distances(p_tree, .x$tip))
  }))
message("Calculating mean distances to different genera (via sampling)...")
bg_dist <- Reduce(c, map(.progress=TRUE, 1:100, \(.) {
  smp <- sample(tips_by_genus$genus, 100)
  smp_tips <- map_chr(filter(tips_by_genus, genus %in% smp)$data, ~ { 
    .x[sample(nrow(.x), 1), "tip"] %>% pull(tip)
  })
  rowMeans(get_all_pairwise_distances(p_tree, only_clades=smp_tips))
}))
message("Filtering...")
lower_bg_lim <- quantile(bg_dist, 0.05)
filtered_means <- all_means %>%
  mutate(below_bg = map(means, ~ names(.x[.x <= lower_bg_lim])))
filtered_tips <- filtered_means %>%
  select(below_bg) %>%
  unnest(below_bg) %>%
  filter(!is.na(below_bg)) %>%
  pull(below_bg)
filtered_tree <- tidytree::keep.tip(tree, filtered_tips)
treeio::write.tree(filtered_tree, p$output)
