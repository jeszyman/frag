#!/usr/bin/env Rscript
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/frag/frag.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-13 15:30:47
# ============================================================

# Arm z-score heatmap.
# Rows = libraries, columns = chromosome arms, fill = z-score.

args <- commandArgs(trailingOnly = TRUE)
armz_csv   <- args[1]
output_pdf <- args[2]

library(tidyverse)

armz <- read.csv(armz_csv, row.names = 1, check.names = FALSE)

armz_long <- armz %>%
  rownames_to_column("library") %>%
  pivot_longer(-library, names_to = "arm", values_to = "z") %>%
  mutate(
    chr_num = as.integer(gsub("chr(\\d+)_.*", "\\1", arm)),
    arm_letter = gsub("chr\\d+_", "", arm)
  ) %>%
  arrange(chr_num, arm_letter)

arm_order <- armz_long %>%
  distinct(arm, chr_num, arm_letter) %>%
  arrange(chr_num, arm_letter) %>%
  pull(arm)

armz_long$arm <- factor(armz_long$arm, levels = arm_order)

z_limit <- max(3, ceiling(quantile(abs(armz_long$z), 0.95, na.rm = TRUE)))

n_libs <- length(unique(armz_long$library))
plot_height <- max(4, 0.4 * n_libs + 2)
plot_width <- max(8, 0.25 * length(arm_order) + 2)

p <- ggplot(armz_long, aes(x = arm, y = library, fill = z)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                        midpoint = 0, limits = c(-z_limit, z_limit),
                        name = "Z-score") +
  theme_minimal() +
  labs(x = "Chromosome arm", y = "Library",
       title = "Arm-level z-scores") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        axis.text.y = element_text(size = 8),
        panel.grid = element_blank())

ggsave(output_pdf, p, width = plot_width, height = plot_height)
