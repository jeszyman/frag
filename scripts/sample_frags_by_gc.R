#!/usr/bin/env Rscript
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/frag/frag.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-24 13:15:21
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
healthy_med <- args[1]
frag_file <- args[2]
sampled_file <- args[3]

library(tidyverse)

healthy_fract <- readRDS(healthy_med)
frag_file <- read.table(frag_file, sep = "\t", header = FALSE)

frag_bed <- frag_file
names(frag_bed) <- c("chr", "start", "end", "gc_raw", "len")

frag <- frag_bed %>%
  mutate(gc_strata = round(gc_raw, 2)) %>%
  left_join(healthy_fract, by = "gc_strata")

stratatotake <- frag$gc_strata[which.max(frag$med_frag_fract)]
fragsinmaxstrata <- length(which(frag$gc_strata == stratatotake))
fragstotake <- round(fragsinmaxstrata / stratatotake)

sampled <- frag %>%
  filter(!is.na(med_frag_fract)) %>%
  slice_sample(., n = nrow(.), weight_by = med_frag_fract, replace = TRUE) %>%
  select(chr, start, end, len, gc_strata)

write.table(sampled, sep = "\t", col.names = FALSE, row.names = FALSE,
            quote = FALSE, file = sampled_file)
