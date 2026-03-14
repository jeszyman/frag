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
frags_tsv <- args[1]
ratios_tsv <- args[2]

library(tidyverse)

frags <- read_tsv(frags_tsv)

ratios <- frags %>%
  mutate_at(vars(start, end, count), as.numeric) %>%
  pivot_wider(names_from = len_class, values_from = count,
              values_fn = function(x) mean(x)) %>%
  mutate(fract = short / long) %>%
  select(library, chr, start, end, fract) %>%
  group_by(library) %>%
  mutate(ratio.centered = scale(fract, scale = FALSE)[, 1])

write_tsv(ratios, file = ratios_tsv)
