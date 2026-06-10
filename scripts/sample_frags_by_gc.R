#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# sample_frags_by_gc.R — GC-weighted resample of fragments toward the healthy
# GC profile (with replacement).
# Input : healthy_med  — RDS of healthy median GC-strata fractions
#         frag_file    — fragment BED (chr, start, end, gc_raw, len)
# Output: sampled_file — resampled fragments (chr, start, end, len, gc_strata)
# NOTE: the resampling is NOT seeded (pre-existing); output is stochastic.
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

# Positional args (repo convention): Snakemake passes inputs then output.
args <- commandArgs(trailingOnly = TRUE)
healthy_med <- args[1]
frag_file <- args[2]
sampled_file <- args[3]

# =============================================================================
# SECTION: WORK FUNCTIONS
# =============================================================================

#' Annotate fragments with GC stratum and join healthy median fractions.
#'
#' @param frag_bed Data frame (chr, start, end, gc_raw, len).
#' @param healthy_fract Tibble (gc_strata, med_frag_fract).
#' @return Tibble of fragments with gc_strata + med_frag_fract.
annotate_gc <- function(frag_bed, healthy_fract) {
  frag_bed %>%
    mutate(gc_strata = round(gc_raw, 2)) %>%
    left_join(healthy_fract, by = "gc_strata")
}

#' GC-weighted resample of fragments (with replacement). Not seeded.
#'
#' @param frag Tibble from annotate_gc().
#' @return Tibble (chr, start, end, len, gc_strata).
resample_by_gc <- function(frag) {
  frag %>%
    filter(!is.na(med_frag_fract)) %>%
    slice_sample(., n = nrow(.), weight_by = med_frag_fract, replace = TRUE) %>%
    select(chr, start, end, len, gc_strata)
}

# =============================================================================
# SECTION: BODY
# =============================================================================

healthy_fract <- readRDS(healthy_med)
frag_file <- read.table(frag_file, sep = "\t", header = FALSE)

frag_bed <- frag_file
names(frag_bed) <- c("chr", "start", "end", "gc_raw", "len")
check_nonempty(frag_bed, "frag_bed")
log_n(frag_bed, "frags read")

frag <- annotate_gc(frag_bed, healthy_fract)

# Strata stats computed in the original but unused by the resample below;
# preserved verbatim to keep behavior identical.
stratatotake <- frag$gc_strata[which.max(frag$med_frag_fract)]
fragsinmaxstrata <- length(which(frag$gc_strata == stratatotake))
fragstotake <- round(fragsinmaxstrata / stratatotake)

sampled <- resample_by_gc(frag)

write.table(sampled, sep = "\t", col.names = FALSE, row.names = FALSE,
            quote = FALSE, file = sampled_file)
