#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# delfi_arm_zscores.R — DELFI arm-level z-scores (Mathios et al. 2021).
# LOESS GC-corrected fragment counts per chromosome arm, z-scored against a
# healthy reference cohort.
# Input : counts_tsv   — per-bin counts (library, chr, start, end, gc, count)
#         cytoband_tsv — cytoBand (chr, start, end, band, stain)
#         healthy_list — file of healthy library_ids, one per line
# Output: output_csv   — libraries x chr_arm matrix of z-scores
# -----------------------------------------------------------------------------

# =============================================================================
# SECTION: PACKAGES
# =============================================================================

packages <- c("tidyverse", "data.table")
suppressPackageStartupMessages(
  invisible(lapply(packages, require, character.only = TRUE))
)
source("scripts/frag_checks.R")

# =============================================================================
# SECTION: ARGUMENT PARSING
# =============================================================================

# Positional args (repo convention): Snakemake passes inputs then output.
args <- commandArgs(trailingOnly = TRUE)
counts_tsv    <- args[1]
cytoband_tsv  <- args[2]
healthy_list  <- args[3]
output_csv    <- args[4]

# =============================================================================
# SECTION: CONSTANTS
# =============================================================================

# Non-acrocentric autosomal arms (acrocentric p-arms lack usable sequence).
autosomes     <- paste0("chr", 1:22)
noacro_chroms <- c(paste0("chr", 1:12), paste0("chr", 16:20))
acro_chroms   <- setdiff(autosomes, noacro_chroms)

# =============================================================================
# SECTION: WORK FUNCTIONS
# =============================================================================

#' Sum short + long counts to a total per bin.
#' @param counts Tibble (library, chr, start, end, gc, count).
#' @return Tibble of per-bin total counts.
summarize_bin_counts <- function(counts) {
  counts %>%
    group_by(library, chr, start, end, gc) %>%
    summarize(count = sum(count), .groups = "drop")
}

#' Build non-acrocentric arm boundaries (wide) from cytobands.
#' @param cytobands Tibble (chr, start, end, band, stain).
#' @return Tibble of armstart_/armend_ per chromosome.
build_noacro <- function(cytobands) {
  cytobands %>%
    mutate(arm = substr(band, 1, 1)) %>%
    group_by(chr, arm) %>%
    summarize(armstart = min(start), armend = max(end), .groups = "drop") %>%
    filter(chr %in% noacro_chroms | (chr %in% acro_chroms & arm == "q")) %>%
    pivot_wider(names_from = arm, values_from = c(armstart, armend))
}

#' LOESS GC correction (from Lucas workflow).
#' @param counts_mult Numeric vector of log-counts.
#' @param gc Numeric vector of bin GC fractions.
#' @return Numeric vector of LOESS-predicted GC trend.
gcCorrectLoess <- function(counts_mult, gc) {
  if (length(counts_mult) < 4) {
    warning("Fewer than 4 bins -- skipping LOESS GC correction")
    return(rep(0, length(counts_mult)))
  }
  upper <- quantile(counts_mult, 0.99, na.rm = TRUE)
  trend_points <- counts_mult > 0 & counts_mult < upper
  trend_counts <- counts_mult[trend_points]
  trend_gc     <- gc[trend_points]
  n_sample     <- min(length(trend_counts), 10000L)
  samp         <- sample(seq_along(trend_counts), n_sample)
  include      <- c(which(gc == min(gc)), which(gc == max(gc)))
  trend_counts <- c(counts_mult[include], trend_counts[samp])
  trend_gc     <- c(gc[include], trend_gc[samp])
  initial_trend <- suppressWarnings(loess(trend_counts ~ trend_gc))
  i <- seq(min(gc, na.rm = TRUE), max(gc, na.rm = TRUE), by = 0.001)
  final_trend <- suppressWarnings(loess(predict(initial_trend, i) ~ i))
  predict(final_trend, gc)
}

#' Annotate bins with arm and filter to non-acrocentric arms.
#' @param counts Tibble of per-bin counts.
#' @param noacro Wide arm-boundary tibble from build_noacro().
#' @return Tibble (library, chr, start, end, gc, arm, count).
annotate_arms <- function(counts, noacro) {
  counts %>%
    filter(chr %in% autosomes) %>%
    left_join(noacro, by = "chr") %>%
    mutate(arm = ifelse(end < armstart_q, "p", "q")) %>%
    filter(chr %in% noacro_chroms | (chr %in% acro_chroms & arm == "q")) %>%
    select(library, chr, start, end, gc, arm, count)
}

#' Per-library LOESS GC correction + centering (data.table, by reference).
#' @param arm_counts Tibble of arm-annotated counts.
#' @return data.table with adjusted_cent column.
loess_correct <- function(arm_counts) {
  arm_counts_DT <- setDT(arm_counts)
  arm_counts_DT[, counts_mult := log2(count + 1)]
  arm_counts_DT[, loess_pred := gcCorrectLoess(counts_mult, gc), by = "library"]
  arm_counts_DT[, adjusted := counts_mult - loess_pred]
  arm_counts_DT[, adjusted_cent := adjusted - median(adjusted, na.rm = TRUE), by = "library"]
  arm_counts_DT
}

#' Aggregate corrected counts to arm-level means, globally centered.
#' @param arm_counts_DT data.table from loess_correct().
#' @return data.table (library, chr, arm, armmean).
aggregate_arms <- function(arm_counts_DT) {
  arms <- arm_counts_DT[, .(armmean = mean(adjusted_cent, na.rm = TRUE)),
                        by = .(library, chr, arm)]
  arms[, armmean := armmean - mean(armmean, na.rm = TRUE)]
  arms
}

#' Z-score arm means against the healthy cohort; return libraries x chr_arm matrix.
#' @param arms_tib Tibble of arm means.
#' @param healthy_libs Character vector of healthy library_ids.
#' @return Data frame with library rownames and one column per chr_arm.
zscore_arms <- function(arms_tib, healthy_libs) {
  healthy_arms <- arms_tib %>%
    filter(library %in% healthy_libs) %>%
    group_by(chr, arm) %>%
    summarize(mean = mean(armmean, na.rm = TRUE),
              sd = sd(armmean, na.rm = TRUE), .groups = "drop")

  arms_tib %>%
    left_join(healthy_arms, by = c("chr", "arm")) %>%
    mutate(z = (armmean - mean) / sd) %>%
    select(-mean, -sd, -armmean) %>%
    unite(chr_arm, chr, arm, sep = "_") %>%
    pivot_wider(names_from = chr_arm, values_from = z) %>%
    column_to_rownames("library")
}

# =============================================================================
# SECTION: BODY
# =============================================================================

counts <- read_tsv(counts_tsv)
check_nonempty(counts, "counts")
assert_cols(counts, c("library", "chr", "start", "end", "gc", "count"), "counts")
log_n(counts, "counts read")

cytobands <- read_tsv(cytoband_tsv,
                      col_names = c("chr", "start", "end", "band", "stain"))
check_nonempty(cytobands, "cytobands")

healthy_libs <- read_lines(healthy_list)

counts <- summarize_bin_counts(counts)
noacro <- build_noacro(cytobands)

# Seed before the LOESS sampling so GC correction is reproducible.
set.seed(42)

arm_counts    <- annotate_arms(counts, noacro)
arm_counts_DT <- loess_correct(arm_counts)
arms          <- aggregate_arms(arm_counts_DT)
arms_tib      <- as_tibble(arms)
arm_z         <- zscore_arms(arms_tib, healthy_libs)

write.csv(arm_z, file = output_csv, row.names = TRUE)
