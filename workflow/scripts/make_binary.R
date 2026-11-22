suppressPackageStartupMessages({
  library(Matrix)
  library(readr)
})

# Convert targets to row indices and cluster to column indices
create_sparse_matrix <- function(df) {
  sparse_matrix <- Matrix::sparseMatrix(
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
if (length(args) != 3) {
  stop("Usage: make_binary.R SPECIES_FAMILIES.csv TAXONOMY.csv OUTPUT.rds")
}

df <- args[1] 
tax <- args[2]
output <- args[3]
df <- as.data.frame(readr::read_csv(df, col_types=cols(.default=col_character()), na=character(), show_col_types=FALSE))
tax <- as.data.frame(readr::read_csv(tax, col_types=cols(.default=col_character()), na=character(), show_col_types=FALSE))

if (!setequal(names(df), c("cluster", "target")) || length(names(df)) != 2) {
  stop("Species/family input must contain only cluster and target columns; regenerate legacy intermediates with reformat_taxonomy.py")
}
df <- unique(df[c("cluster", "target")])
if (nrow(df) == 0 || anyNA(df) || any(trimws(df$cluster) == "") || any(trimws(df$target) == "")) {
  stop("Species/family input contains empty or missing IDs")
}
if (!all(c("cluster", "phylum") %in% names(tax))) {
  stop("Taxonomy must contain cluster and phylum columns")
}
if (anyNA(tax[c("cluster", "phylum")]) || any(trimws(tax$cluster) == "") || any(trimws(tax$phylum) == "")) {
  stop("Taxonomy contains missing species IDs or phyla")
}
if (anyDuplicated(tax$cluster) || !setequal(df$cluster, tax$cluster)) {
  stop("Species/family input and taxonomy must contain exactly the same unique species IDs")
}

result <- merge(df, tax[c("cluster", "phylum")], by = "cluster", all.x = TRUE)

phylum_split <- split(result, result$phylum)
phylum_split <- lapply(phylum_split, function(df) df[, c("cluster", "target")])
sparse_matrices <- lapply(phylum_split, create_sparse_matrix)

dir.create(dirname(output), recursive=TRUE, showWarnings=FALSE)
saveRDS(sparse_matrices, output)
