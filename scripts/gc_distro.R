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
bed_file <- args[1]
distro_file <- args[2]

library(tidyverse)

bed <- read.table(bed_file, sep = "\t")
names(bed) <- c("chr", "start", "end", "gc_raw", "len")

distro <- bed %>%
  mutate(gc_strata = round(gc_raw, 2)) %>%
  count(gc_strata) %>%
  mutate(fract_frags = n / sum(n)) %>%
  mutate(library_id = gsub("_frag.bed", "", gsub("^.*lib", "lib", bed_file))) %>%
  select(library_id, gc_strata, fract_frags)

write.csv(distro, file = distro_file, row.names = FALSE)
