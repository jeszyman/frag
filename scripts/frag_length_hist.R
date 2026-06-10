#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# frag_length_hist.R
# Fragment-length histogram (count per length) within [start_bp, end_bp], zero-filled.
# Input : frag_bed — fragment BED (chr, start, end, gc, length)
# Output: out_tsv  — length, count
# -----------------------------------------------------------------------------

# =============================================================================
# SECTION: PACKAGES
# =============================================================================

packages <- c("data.table")
suppressPackageStartupMessages(
  invisible(lapply(packages, require, character.only = TRUE))
)
source("scripts/frag_checks.R")

# =============================================================================
# SECTION: ARGUMENT PARSING
# =============================================================================

# Positional args (repo convention): Snakemake passes inputs then outputs.
args <- commandArgs(trailingOnly = TRUE)
frag_bed <- args[1]
out_tsv  <- args[2]
start_bp <- as.integer(args[3])
end_bp   <- as.integer(args[4])

# =============================================================================
# SECTION: WORK FUNCTIONS
# =============================================================================

#' Histogram of fragment lengths within [start_bp, end_bp], zero-filled.
#'
#' @param frags data.table with column `length`.
#' @param start_bp,end_bp Integer inclusive length bounds.
#' @return data.table (length, count) for every length in range.
compute_length_hist <- function(frags, start_bp, end_bp) {
  hist_dt <- frags[length >= start_bp & length <= end_bp, .N, by = length]
  setnames(hist_dt, "N", "count")
  all_lengths <- data.table(length = start_bp:end_bp)
  hist_dt <- merge(all_lengths, hist_dt, by = "length", all.x = TRUE)
  hist_dt[is.na(count), count := 0L]
  setorder(hist_dt, length)
  hist_dt
}

# =============================================================================
# SECTION: BODY
# =============================================================================

frags <- fread(frag_bed, col.names = c("chr", "start", "end", "gc", "length"))
check_nonempty(frags, "frags")
log_n(frags, "frags read")

hist_dt <- compute_length_hist(frags, start_bp, end_bp)
log_n(hist_dt, "length histogram")

fwrite(hist_dt, out_tsv, sep = "\t")
