#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
motif_str <- args[1]
motif_tsv <- args[2]

library(tidyverse)

possible_motifs <- expand.grid(rep(list(c("A", "G", "T", "C")), 4)) %>%
  as_tibble() %>%
  mutate(motif = paste0(Var1, Var2, Var3, Var4)) %>%
  select(motif) %>%
  arrange(motif)

motif_files <- strsplit(motif_str, " ")[[1]]
names(motif_files) <- substr(gsub("^.*lib", "lib", motif_files), 1, 6)

ingest_motif <- function(motif_file) {
  read_tsv(motif_file, col_names = c("motif")) %>%
    group_by(motif) %>%
    summarise(count = n()) %>%
    mutate(fract = count / sum(count)) %>%
    select(motif, fract)
}

motif_tibs <- lapply(motif_files, ingest_motif)

motifs <- bind_rows(motif_tibs, .id = "library") %>%
  pivot_wider(names_from = library, values_from = fract) %>%
  filter(motif %in% possible_motifs$motif)

write_tsv(motifs, motif_tsv)
