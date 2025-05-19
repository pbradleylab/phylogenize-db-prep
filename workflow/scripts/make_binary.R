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

df <- args[1] 
tax <- args[2]
output <- args[3]
df <- read.csv(df, sep=",")
tax <- read.csv(tax, sep=",")

result <- merge(df, tax[c("cluster", "phylum")], by = "cluster", all.x = TRUE)

tmp <- unique(result[c("query","phylum","cluster")])
names(tmp)[names(tmp) == "query"] <- "other"
tmp <- merge(result[c("query","target","other")], tmp, by="other", all.x=TRUE)

phylum_split <- split(tmp, tmp$phylum)
phylum_split <- lapply(phylum_split, function(df) df[, c("cluster", "target")])
phylum_split <- phylum_split[sapply(phylum_split, function(df) nrow(df) > 1)]
sparse_matrices <- lapply(phylum_split, create_sparse_matrix)

saveRDS(sparse_matrices, output)
