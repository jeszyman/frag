#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# ratio_row_normalize.R — per-library z-score of the DELFI ratio profile across
# genomic bins (zero mean, unit SD per library).
# Input : ratios_tsv — long DELFI ratios (library, chr, start, end, fract)
# Output: output_csv — library x bin row-normalized matrix
# -----------------------------------------------------------------------------

# =============================================================================
# SECTION: PACKAGES
# =============================================================================

packages <- c("tidyverse")
suppressPackageStartupMessages(
  invisible(lapply(packages, require, character.only = TRUE))
)
source("scripts/frag_checks.R")

# =============================================================================
# SECTION: ARGUMENT PARSING
# =============================================================================

# Positional args (repo convention): Snakemake passes input then output.
args <- commandArgs(trailingOnly = TRUE)
ratios_tsv <- args[1]
output_csv <- args[2]

# =============================================================================
# SECTION: WORK FUNCTIONS
# =============================================================================

#' Row-normalize (per-library z-score) the DELFI ratio profile across bins.
#'
#' @param ratios Tibble (library, chr, start, end, fract).
#' @return Matrix of row-normalized ratios (library rows, bin columns).
compute_ratio_norm <- function(ratios) {
  # Pivot to library x bin matrix
  ratio_mat <- ratios %>%
    unite(bin, chr, start, end, sep = "_") %>%
    select(library, bin, fract) %>%
    pivot_wider(names_from = bin, values_from = fract) %>%
    column_to_rownames("library")

  # Drop columns with any NA, drop non-autosomal
  ratio_mat <- ratio_mat[, !apply(ratio_mat, 2, anyNA)]
  ratio_mat <- ratio_mat[, !grepl("^chrX_|^chrY_|^chrM_", colnames(ratio_mat))]

  # Row-wise z-score normalization
  t(scale(t(ratio_mat)))
}

# =============================================================================
# SECTION: BODY
# =============================================================================

ratios <- read_tsv(ratios_tsv)
check_nonempty(ratios, "ratios")
assert_cols(ratios, c("library", "chr", "start", "end", "fract"), "ratios")
log_n(ratios, "ratios read")

ratio_norm <- compute_ratio_norm(ratios)

write.csv(as.data.frame(ratio_norm), file = output_csv, row.names = TRUE)
