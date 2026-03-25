#!/usr/bin/env Rscript
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/frag/frag.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-25 10:04:32
# ============================================================

# Genome-wide DELFI ratio profile plot.
# One panel per library, ratio.centered across genomic bins.

args <- commandArgs(trailingOnly = TRUE)
ratios_tsv <- args[1]
output_pdf <- args[2]

source("~/repos/science/R/figure_schema.R")
if (!"Arial" %in% names(grDevices::pdfFonts()))
  theme_scifig <- function(base_size = BASE_SIZE) {
    theme_bw(base_size = base_size) +
      theme(text = element_text(family = "Helvetica"),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            strip.background = element_rect(fill = "grey95", colour = NA),
            legend.background = element_rect(fill = NA), legend.key = element_rect(fill = NA),
            plot.title = element_blank(), plot.margin = margin(4, 4, 4, 4, "pt"))
  }
library(tidyverse)

ratios <- read_tsv(ratios_tsv)

# Create genomic coordinate for x-axis (numeric chromosome sort)
ratios <- ratios %>%
  mutate(chr_num = as.integer(gsub("chr", "", chr))) %>%
  filter(!is.na(chr_num)) %>%
  arrange(library, chr_num, start) %>%
  group_by(library) %>%
  mutate(genome_pos = row_number()) %>%
  ungroup()

# Alternating chromosome shading
chr_bounds <- ratios %>%
  group_by(library, chr_num) %>%
  summarize(xmin = min(genome_pos), xmax = max(genome_pos), .groups = "drop") %>%
  distinct(chr_num, xmin, xmax) %>%
  mutate(shade = chr_num %% 2 == 0)

n_libs <- length(unique(ratios$library))
plot_height <- max(4, 2 * n_libs)

p <- ggplot(ratios, aes(x = genome_pos, y = ratio.centered)) +
  geom_rect(data = filter(chr_bounds, shade),
            aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "grey90", alpha = 0.5) +
  geom_line(linewidth = 0.3, alpha = 0.8, color = "steelblue") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  facet_wrap(~library, ncol = 1, scales = "fixed") +
  labs(x = "Genomic position", y = "Centered ratio (short/long)") +
  theme_scifig() +
  theme(strip.background = element_rect(fill = "white"),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(output_pdf, p, width = 12, height = plot_height)
