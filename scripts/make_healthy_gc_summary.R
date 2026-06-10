#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# make_healthy_gc_summary.R — median GC-strata fragment fraction across healthy libs.
# Input : healthy_libs_str — space-separated paths to per-library gc_distro CSVs
# Output: healthy_med_file — RDS of (gc_strata, med_frag_fract)
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
healthy_libs_str <- args[1]
healthy_med_file <- args[2]

# =============================================================================
# SECTION: WORK FUNCTIONS
# =============================================================================

#' Read one per-library GC distro CSV.
#'
#' @param gc_csv Path to a gc_distro CSV.
#' @return Data frame (library_id, gc_strata, fract_frags).
read_in_gc <- function(gc_csv) {
  read.csv(gc_csv, header = TRUE)
}

#' Median fragment fraction per GC stratum across libraries.
#'
#' @param healthy_all Combined per-library GC distro rows.
#' @return Tibble (gc_strata, med_frag_fract).
summarize_healthy_gc <- function(healthy_all) {
  healthy_all %>%
    group_by(gc_strata) %>%
    summarise(med_frag_fract = median(fract_frags))
}

# =============================================================================
# SECTION: BODY
# =============================================================================

healthy_libs_distros <- unlist(strsplit(healthy_libs_str, " "))

healthy_list <- lapply(healthy_libs_distros, read_in_gc)
healthy_all  <- do.call(rbind, healthy_list)
check_nonempty(healthy_all, "healthy_all")
assert_cols(healthy_all, c("gc_strata", "fract_frags"), "healthy_all")
log_n(healthy_all, "healthy gc rows")

healthy_med <- summarize_healthy_gc(healthy_all)
log_n(healthy_med, "healthy med")

saveRDS(healthy_med, file = healthy_med_file)
