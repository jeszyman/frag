#!/usr/bin/env Rscript
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/frag/frag.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-25 09:43:51
# ============================================================

# Fragment length distribution overlay plot.
# Plots all libraries' frequency profiles, optionally colored by group.

args <- commandArgs(trailingOnly = TRUE)
freq_csv    <- args[1]
output_pdf  <- args[2]
samples_tsv <- if (length(args) >= 3) args[3] else NULL

source("~/repos/science/R/figure_schema.R")
library(tidyverse)

freq_mat <- read.csv(freq_csv, row.names = 1, check.names = FALSE)

freq_long <- freq_mat %>%
  rownames_to_column("library") %>%
  pivot_longer(-library, names_to = "length", values_to = "freq") %>%
  mutate(length = as.integer(length)) %>%
  filter(length <= 400)

if (!is.null(samples_tsv) && file.exists(samples_tsv)) {
  samples <- read_tsv(samples_tsv, col_types = cols(
    library_id = col_character(), group = col_character()
  ))
  freq_long <- freq_long %>%
    left_join(samples, by = c("library" = "library_id"))
  color_aes <- "group"
} else {
  color_aes <- "library"
}

p <- ggplot(freq_long, aes(x = length, y = freq,
                            color = .data[[color_aes]], group = library)) +
  geom_vline(xintercept = 167, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 334, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_line(linewidth = 0.7, alpha = 0.85) +
  scale_color_brewer(palette = "Dark2") +
  scale_x_continuous(breaks = c(100, 200, 300, 400)) +
  labs(x = "Fragment length (bp)", y = "Frequency", color = NULL) +
  theme_scifig() +
  theme(legend.position.inside = c(0.85, 0.7)) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))

ggsave(output_pdf, p, width = PLOT_W, height = PLOT_H)
