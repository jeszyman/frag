#!/usr/bin/env Rscript
# Motif diversity score: Shannon entropy normalized by log(256).
# Uses data.table for efficient I/O.

args <- commandArgs(trailingOnly = TRUE)
motif_tsv  <- args[1]
output_csv <- args[2]

library(data.table)

motifs <- fread(motif_tsv)
lib_cols <- names(motifs)[-1]

mds_list <- lapply(lib_cols, function(lib) {
  fracs <- motifs[[lib]]
  fracs[is.na(fracs)] <- 0
  fracs[fracs == 0] <- 1e-10
  data.table(library = lib, mds = -sum(fracs * log(fracs)) / log(256))
})

mds_dt <- rbindlist(mds_list)
fwrite(mds_dt, output_csv)
