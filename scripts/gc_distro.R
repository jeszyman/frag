#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# gc_distro.R — per-library fragment GC-strata distribution.
# Input : bed_file    — fragment BED (chr, start, end, gc_raw, len)
# Output: distro_file — library_id, gc_strata, fract_frags
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
bed_file <- args[1]
distro_file <- args[2]

# =============================================================================
# SECTION: WORK FUNCTIONS
# =============================================================================

#' GC-strata fragment distribution for one library.
#'
#' @param bed Data frame (chr, start, end, gc_raw, len).
#' @param bed_file Path of the BED (library_id parsed from its basename).
#' @return Tibble (library_id, gc_strata, fract_frags).
compute_gc_distro <- function(bed, bed_file) {
  bed %>%
    mutate(gc_strata = round(gc_raw, 2)) %>%
    count(gc_strata) %>%
    mutate(fract_frags = n / sum(n)) %>%
    mutate(library_id = sub("\\.[^.]+\\.frag\\.bed$", "", basename(bed_file))) %>%
    select(library_id, gc_strata, fract_frags)
}

# =============================================================================
# SECTION: BODY
# =============================================================================

bed <- read.table(bed_file, sep = "\t")
names(bed) <- c("chr", "start", "end", "gc_raw", "len")
check_nonempty(bed, "bed")
assert_cols(bed, c("chr", "start", "end", "gc_raw", "len"), "bed")
log_n(bed, "bed read")

distro <- compute_gc_distro(bed, bed_file)
log_n(distro, "gc distro")

write.csv(distro, file = distro_file, row.names = FALSE)
