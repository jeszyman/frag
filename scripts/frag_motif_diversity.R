#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# frag_motif_diversity.R — per-library end-motif Shannon-entropy diversity score.
# Input : motif_tsv  — motif x library fraction matrix
# Output: output_csv — library, mds (entropy normalized by log(256))
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

# Positional args (repo convention): Snakemake passes input then output.
args <- commandArgs(trailingOnly = TRUE)
motif_tsv  <- args[1]
output_csv <- args[2]

# =============================================================================
# SECTION: WORK FUNCTIONS
# =============================================================================

#' Shannon-entropy motif diversity per library, normalized by log(256).
#'
#' @param motifs data.table; column 1 is the motif, remaining columns libraries.
#' @return data.table (library, mds).
compute_motif_diversity <- function(motifs) {
  lib_cols <- names(motifs)[-1]
  mds_list <- lapply(lib_cols, function(lib) {
    fracs <- motifs[[lib]]
    fracs[is.na(fracs)] <- 0
    fracs[fracs == 0] <- 1e-10
    data.table(library = lib, mds = -sum(fracs * log(fracs)) / log(256))
  })
  rbindlist(mds_list)
}

# =============================================================================
# SECTION: BODY
# =============================================================================

motifs <- fread(motif_tsv)
check_nonempty(motifs, "motifs")
log_n(motifs, "motifs read")

mds_dt <- compute_motif_diversity(motifs)
log_n(mds_dt, "motif diversity")

fwrite(mds_dt, output_csv)
