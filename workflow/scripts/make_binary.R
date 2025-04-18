library(Matrix)

# Convert targets to row indices and cluster to column indices
create_sparse_matrix <- function(df) {
  sparse_matrix <- sparseMatrix(
    i = as.numeric(factor(df$target)),
    j = as.numeric(factor(df$cluster)),
    x = TRUE,
    dimnames = list(levels(factor(df$target)), levels(factor(df$cluster)))
  )
  # Convert to a logical sparse matrix (TRUE/FALSE)
  as(sparse_matrix, "lgCMatrix")
}


# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

df <- args[1]; output <- args[2]
df <- read.csv(df, sep=",")

phylum_split <- split(df, df$phylum)
phylum_split <- lapply(phylum_split, function(df) df[, c("cluster", "target")])
phylum_split <- phylum_split[sapply(phylum_split, function(df) nrow(df) > 1)]
sparse_matrices <- lapply(phylum_split, create_sparse_matrix)

saveRDS(sparse_matrices, output)
