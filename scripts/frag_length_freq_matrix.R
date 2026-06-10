#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# frag_length_freq_matrix.R — library x length count and frequency matrices.
# Input : hist_files — space-separated per-library length_hist TSVs
# Output: out_counts — library x length count matrix (CSV)
#         out_freqs  — row-normalized frequency matrix (CSV)
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

# Positional args (repo convention): Snakemake passes inputs then outputs.
args <- commandArgs(trailingOnly = TRUE)
hist_files <- args[1]
out_counts <- args[2]
out_freqs  <- args[3]

# =============================================================================
# SECTION: WORK FUNCTIONS
# =============================================================================

#' Read and stack per-library length histograms, tagging each with its library.
#'
#' @param files Character vector of length_hist TSV paths.
#' @return Tibble (length, count, library).
read_all_hists <- function(files) {
  map_dfr(files, function(f) {
    lib_id <- str_extract(basename(f), "^[^.]+")
    read_tsv(f, col_types = "ii") %>%
      mutate(library = lib_id)
  })
}

#' Pivot stacked histograms to a library x length count matrix.
#'
#' @param all_hists Tibble (length, count, library).
#' @return Data frame with library rownames, one column per length.
build_count_mat <- function(all_hists) {
  all_hists %>%
    pivot_wider(names_from = length, values_from = count, values_fill = 0L) %>%
    column_to_rownames("library")
}

# =============================================================================
# SECTION: BODY
# =============================================================================

files <- str_split(hist_files, " ")[[1]]

all_hists <- read_all_hists(files)
check_nonempty(all_hists, "all_hists")
assert_cols(all_hists, c("length", "count", "library"), "all_hists")
log_n(all_hists, "histograms stacked")

count_mat <- build_count_mat(all_hists)
write.csv(count_mat, out_counts)

freq_mat <- count_mat / rowSums(count_mat)
write.csv(freq_mat, out_freqs)
