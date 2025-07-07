library(Matrix)
library(ape)


# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

tree <- args[1]; output <- args[2]
tree <- ape::read.tree(tree)
saveRDS(tree, output)
