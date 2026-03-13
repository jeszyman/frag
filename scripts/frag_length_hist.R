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
frag_bed <- args[1]
out_tsv  <- args[2]
start_bp <- as.integer(args[3])
end_bp   <- as.integer(args[4])

library(data.table)

frags <- fread(frag_bed, col.names = c("chr", "start", "end", "gc", "length"))

hist_dt <- frags[length >= start_bp & length <= end_bp, .N, by = length]
setnames(hist_dt, "N", "count")
all_lengths <- data.table(length = start_bp:end_bp)
hist_dt <- merge(all_lengths, hist_dt, by = "length", all.x = TRUE)
hist_dt[is.na(count), count := 0L]
setorder(hist_dt, length)

fwrite(hist_dt, out_tsv, sep = "\t")
