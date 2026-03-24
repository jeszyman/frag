#!/usr/bin/env Rscript
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/frag/frag.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-24 13:15:21
# ============================================================

# Row-normalize DELFI ratios: z-score each library's ratio profile
# across genomic bins (zero mean, unit SD per library).

args <- commandArgs(trailingOnly = TRUE)
ratios_tsv <- args[1]
output_csv <- args[2]

library(tidyverse)

ratios <- read_tsv(ratios_tsv)

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
ratio_norm <- t(scale(t(ratio_mat)))

write.csv(as.data.frame(ratio_norm), file = output_csv, row.names = TRUE)
