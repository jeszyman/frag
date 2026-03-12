#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
healthy_libs_str <- args[1]
healthy_med_file <- args[2]

library(tidyverse)

healthy_libs_distros <- unlist(strsplit(healthy_libs_str, " "))

read_in_gc <- function(gc_csv) {
  read.csv(gc_csv, header = TRUE)
}

healthy_list <- lapply(healthy_libs_distros, read_in_gc)

healthy_all <- do.call(rbind, healthy_list)

healthy_med <- healthy_all %>%
  group_by(gc_strata) %>%
  summarise(med_frag_fract = median(fract_frags))

saveRDS(healthy_med, file = healthy_med_file)
