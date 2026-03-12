#!/usr/bin/env Rscript
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
