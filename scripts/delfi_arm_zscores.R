#!/usr/bin/env Rscript
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/frag/frag.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-24 13:15:21
# ============================================================

# DELFI arm-level z-scores (Mathios et al. 2021)
# Computes LOESS GC-corrected fragment counts per chromosome arm,
# then z-scores against healthy reference cohort.

args <- commandArgs(trailingOnly = TRUE)
counts_tsv    <- args[1]
cytoband_tsv  <- args[2]
healthy_list  <- args[3]
output_csv    <- args[4]

library(tidyverse)
library(data.table)

counts    <- read_tsv(counts_tsv)
cytobands <- read_tsv(cytoband_tsv,
                       col_names = c("chr", "start", "end", "band", "stain"))
healthy_libs <- read_lines(healthy_list)

# Sum short + long to get total counts per bin
counts <- counts %>%
  group_by(library, chr, start, end, gc) %>%
  summarize(count = sum(count), .groups = "drop")

# Define non-acrocentric arms
autosomes     <- paste0("chr", 1:22)
noacro_chroms <- c(paste0("chr", 1:12), paste0("chr", 16:20))
acro_chroms   <- setdiff(autosomes, noacro_chroms)

noacro <- cytobands %>%
  mutate(arm = substr(band, 1, 1)) %>%
  group_by(chr, arm) %>%
  summarize(armstart = min(start), armend = max(end), .groups = "drop") %>%
  filter(chr %in% noacro_chroms | (chr %in% acro_chroms & arm == "q")) %>%
  pivot_wider(names_from = arm, values_from = c(armstart, armend))

# LOESS GC correction (from Lucas workflow)
set.seed(42)

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

# Annotate arms, filter to non-acrocentric
arm_counts <- counts %>%
  filter(chr %in% autosomes) %>%
  left_join(noacro, by = "chr") %>%
  mutate(arm = ifelse(end < armstart_q, "p", "q")) %>%
  filter(chr %in% noacro_chroms | (chr %in% acro_chroms & arm == "q")) %>%
  select(library, chr, start, end, gc, arm, count)

# LOESS correction per library
arm_counts_DT <- setDT(arm_counts)
arm_counts_DT[, counts_mult := log2(count + 1)]
arm_counts_DT[, loess_pred := gcCorrectLoess(counts_mult, gc), by = "library"]
arm_counts_DT[, adjusted := counts_mult - loess_pred]
arm_counts_DT[, adjusted_cent := adjusted - median(adjusted, na.rm = TRUE), by = "library"]

# Aggregate to arm level
arms <- arm_counts_DT[, .(armmean = mean(adjusted_cent, na.rm = TRUE)),
                       by = .(library, chr, arm)]
arms[, armmean := armmean - mean(armmean, na.rm = TRUE)]

arms_tib <- as_tibble(arms)

# Z-score against healthy cohort
healthy_arms <- arms_tib %>%
  filter(library %in% healthy_libs) %>%
  group_by(chr, arm) %>%
  summarize(mean = mean(armmean, na.rm = TRUE),
            sd = sd(armmean, na.rm = TRUE), .groups = "drop")

arm_z <- arms_tib %>%
  left_join(healthy_arms, by = c("chr", "arm")) %>%
  mutate(z = (armmean - mean) / sd) %>%
  select(-mean, -sd, -armmean) %>%
  unite(chr_arm, chr, arm, sep = "_") %>%
  pivot_wider(names_from = chr_arm, values_from = z) %>%
  column_to_rownames("library")

write.csv(arm_z, file = output_csv, row.names = TRUE)
