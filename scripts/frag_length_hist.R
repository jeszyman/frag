#!/usr/bin/env Rscript
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/frag/frag.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-13 15:30:47
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
frag_bed <- args[1]
out_tsv  <- args[2]
start_bp <- as.integer(args[3])
end_bp   <- as.integer(args[4])

library(tidyverse)

frags <- read_tsv(frag_bed, col_names = c("chr", "start", "end", "gc", "length"),
                  col_types = "ciidi")

hist_df <- frags %>%
  filter(length >= start_bp, length <= end_bp) %>%
  count(length, name = "count") %>%
  right_join(tibble(length = start_bp:end_bp), by = "length") %>%
  replace_na(list(count = 0L)) %>%
  arrange(length)

write_tsv(hist_df, out_tsv)
