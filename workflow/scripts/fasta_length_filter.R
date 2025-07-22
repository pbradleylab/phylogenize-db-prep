#!/usr/bin/env Rscript
# Filter a FASTA file by length using a variety of possible cutoffs. If multiple cutoffs apply, the strictest will be chosen
library(optparse)
library(seqinr)
library(purrr)

opt_list <- list(
  make_option(c("-i", "--input"), type="character", help="path to FASTA file to filter"),
  make_option(c("-l", "--len"), type="integer", default = 0, help="minimum length in bases (inclusive)"),
  make_option(c("-q", "--quantile"), type="numeric", default = 0, help="minimum quantile of lengths (inclusive)"),
  make_option(c("-s", "--stdev"), type="numeric", default = 10, help="minimum as number of standard deviations below the mean (inclusive)"),
  make_option(c("-f", "--frac_max"), type="numeric", default = 0, help="minimum as fraction of max-length sequence"),
  make_option(c("-u", "--quant_max"), type="numeric", default = 1, help="quantile that counts as the maximum (for frac_max)"),
  make_option(c("-o", "--output"), type="character", help="path to output file")
)
prs <- OptionParser(option_list = opt_list)
p <- parse_args(prs)

fa <- read.fasta(p$input, seqtype="DNA", as.string=TRUE, forceDNAtolower = FALSE)
lengths <- map_int(fa, nchar)
min_quant <- quantile(na.omit(lengths), p$quantile)
mean_len <- mean(na.omit(lengths))
sd_len <- sd(na.omit(lengths))
by_sd_cutoff <- floor(max(mean_len - (p$stdev * sd_len), 0))
max_for_frac <- quantile(na.omit(lengths), p$quant_max)
frac_max_cutoff <- floor(max_for_frac * p$frac_max)
cutoffs <- c(length=p$len,
  quantile=min_quant,
  stdev=by_sd_cutoff,
  frac_max=frac_max_cutoff)
strictest_cutoff <- names(which.max(cutoffs))
message(
  paste0("Strictest cutoff was ", strictest_cutoff,
  " ( ", cutoffs[strictest_cutoff], " bases)"))
fa <- fa[(lengths >= max(cutoffs))]
message(
  paste0(length(lengths) - length(fa),
  " out of ", length(lengths),
  " sequences filtered")
)
if (length(fa) == 0) {
  warning("Warning: no sequences passed the cutoff; file will be empty")
}
write.fasta(fa, names(fa), file.out = p$output)
