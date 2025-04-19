library(Matrix)
library(ape)


# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

df <- args[1]; tree <- args[2]; output <- args[3]
df <- read.csv(df, sep=",")
tree <- ape::read.tree(tree)

match <- match(tree$tip.label, df$Species_rep)
tree$tip.label <- df$cluster[match]

saveRDS(tree, output)
