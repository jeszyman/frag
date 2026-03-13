#!/usr/bin/env Rscript
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/frag/frag.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-13 15:30:47
# ============================================================

# Fragment length distribution overlay plot.
# Plots all libraries' frequency profiles, optionally colored by group.

args <- commandArgs(trailingOnly = TRUE)
freq_csv    <- args[1]
output_pdf  <- args[2]
samples_tsv <- if (length(args) >= 3) args[3] else NULL

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
} else {
  freq_long$group <- "all"
}

p <- ggplot(freq_long, aes(x = length, y = freq,
                            color = group, group = library)) +
  geom_line(linewidth = 0.3, alpha = 0.4) +
  theme_bw() +
  labs(x = "Fragment length (bp)", y = "Frequency",
       title = "Fragment length distributions") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.85, 0.85),
        legend.background = element_rect(linewidth = 0.3, color = "black")) +
  guides(color = guide_legend(override.aes = list(alpha = 1, linewidth = 1)))

if (length(unique(freq_long$group)) == 1) {
  p <- p + scale_color_manual(values = "grey30") + theme(legend.position = "none")
}

ggsave(output_pdf, p, width = 8, height = 5)
