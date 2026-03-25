#!/usr/bin/env Rscript
# Arm z-score heatmap.
# Rows = libraries (clustered), columns = chromosome arms, fill = z-score.

args <- commandArgs(trailingOnly = TRUE)
armz_csv   <- args[1]
output_pdf <- args[2]

schema_path <- path.expand("~/repos/science/R/figure_schema.R")
if (file.exists(schema_path)) {
  source(schema_path)
} else {
  PLOT_DPI <- 300; PLOT_W <- 6; PLOT_H <- 4.5; BASE_SIZE <- 18
}
theme_scifig <- function(base_size = BASE_SIZE) {
  fam <- if ("Arial" %in% names(grDevices::pdfFonts())) "Arial" else "Helvetica"
  theme_bw(base_size = base_size) +
    theme(text = element_text(family = fam),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "grey95", colour = NA),
          legend.background = element_rect(fill = NA), legend.key = element_rect(fill = NA),
          plot.title = element_blank(), plot.margin = margin(4, 4, 4, 4, "pt"))
}
library(tidyverse)

armz <- read.csv(armz_csv, row.names = 1, check.names = FALSE)

# Cluster libraries by z-score profile
if (nrow(armz) > 2) {
  lib_order <- rownames(armz)[hclust(dist(armz))$order]
} else {
  lib_order <- rownames(armz)
}

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
armz_long$library <- factor(armz_long$library, levels = lib_order)

z_limit <- max(abs(armz_long$z), na.rm = TRUE)

p <- ggplot(armz_long, aes(x = arm, y = library, fill = z)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "firebrick",
                        midpoint = 0, limits = c(-z_limit, z_limit),
                        oob = scales::squish,
                        na.value = "grey80",
                        name = "Z-score") +
  coord_fixed() +
  labs(x = "Chromosome arm", y = NULL) +
  theme_scifig(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 9),
        legend.key.height = unit(0.8, "cm"),
        legend.key.width = unit(0.4, "cm"))

ggsave(output_pdf, p, width = 12, height = 4)
