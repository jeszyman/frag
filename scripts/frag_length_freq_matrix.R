#!/usr/bin/env Rscript
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/frag/frag.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-13 16:35:03
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
hist_files <- args[1]
out_counts <- args[2]
out_freqs  <- args[3]

library(tidyverse)

files <- str_split(hist_files, " ")[[1]]

all_hists <- map_dfr(files, function(f) {
  lib_id <- str_extract(basename(f), "^[^.]+")
  read_tsv(f, col_types = "ii") %>%
    mutate(library = lib_id)
})

count_mat <- all_hists %>%
  pivot_wider(names_from = length, values_from = count, values_fill = 0L) %>%
  column_to_rownames("library")

write.csv(count_mat, out_counts)

freq_mat <- count_mat / rowSums(count_mat)

write.csv(freq_mat, out_freqs)
