#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# make_ratios.R
# Per-bin DELFI short/long fragment ratio, centered within each library.
# Input : frags_tsv  — long counts (library, len_class, chr, start, end, gc, count)
# Output: ratios_tsv — library, chr, start, end, fract, ratio.centered
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
frags_tsv  <- args[1]
ratios_tsv <- args[2]
# Interactive testing: comment out the block above and hard-code paths here.
# frags_tsv  <- "tests/full/frags/chr22.frag_counts.tsv"
# ratios_tsv <- "tests/full/frags/chr22.ratios.tsv"

# =============================================================================
# SECTION: WORK FUNCTIONS
# =============================================================================

#' Compute per-bin short/long ratio, centered within each library.
#'
#' @param frags Tibble with columns library, len_class, chr, start, end, count.
#' @return Tibble with columns library, chr, start, end, fract, ratio.centered.
compute_ratios <- function(frags) {
  frags %>%
    mutate_at(vars(start, end, count), as.numeric) %>%
    pivot_wider(names_from = len_class, values_from = count,
                values_fn = function(x) mean(x)) %>%
    mutate(fract = short / long) %>%
    select(library, chr, start, end, fract) %>%
    group_by(library) %>%
    mutate(ratio.centered = scale(fract, scale = FALSE)[, 1]) %>%
    ungroup()
}

# =============================================================================
# SECTION: BODY
# =============================================================================

frags <- read_tsv(frags_tsv)
check_nonempty(frags, "frags")
assert_cols(frags, c("library", "len_class", "chr", "start", "end", "count"), "frags")
log_n(frags, "frags read")

ratios <- compute_ratios(frags)
log_n(ratios, "ratios computed")

write_tsv(ratios, file = ratios_tsv)
