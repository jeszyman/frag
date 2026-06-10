#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# plot_ratio_profile.R — genome-wide DELFI ratio profile, one panel per library.
# Input : ratios_tsv — long DELFI ratios (library, chr, start, ..., ratio.centered)
# Output: output_pdf — faceted ratio profile (PDF)
# -----------------------------------------------------------------------------

# =============================================================================
# SECTION: PACKAGES
# =============================================================================

packages <- c("tidyverse")
suppressPackageStartupMessages(
  invisible(lapply(packages, require, character.only = TRUE))
)
source("scripts/frag_checks.R")

# =============================================================================
# SECTION: ARGUMENT PARSING
# =============================================================================

# Positional args (repo convention): Snakemake passes input then output.
args <- commandArgs(trailingOnly = TRUE)
ratios_tsv <- args[1]
output_pdf <- args[2]

# =============================================================================
# SECTION: FIGURE THEME
# =============================================================================

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

# =============================================================================
# SECTION: WORK FUNCTIONS
# =============================================================================

#' Add a per-library genomic x-coordinate (numeric chromosome sort).
#'
#' @param ratios Tibble of long DELFI ratios.
#' @return Tibble with chr_num and genome_pos columns.
add_genome_pos <- function(ratios) {
  ratios %>%
    mutate(chr_num = as.integer(gsub("chr", "", chr))) %>%
    filter(!is.na(chr_num)) %>%
    arrange(library, chr_num, start) %>%
    group_by(library) %>%
    mutate(genome_pos = row_number()) %>%
    ungroup()
}

#' Alternating-chromosome shading bounds.
#'
#' @param ratios Tibble from add_genome_pos().
#' @return Tibble (chr_num, xmin, xmax, shade).
compute_chr_bounds <- function(ratios) {
  ratios %>%
    group_by(library, chr_num) %>%
    summarize(xmin = min(genome_pos), xmax = max(genome_pos), .groups = "drop") %>%
    distinct(chr_num, xmin, xmax) %>%
    mutate(shade = chr_num %% 2 == 0)
}

#' Build the faceted ratio-profile plot.
#'
#' @param ratios Tibble from add_genome_pos().
#' @param chr_bounds Tibble from compute_chr_bounds().
#' @return A ggplot object.
build_ratio_plot <- function(ratios, chr_bounds) {
  ggplot(ratios, aes(x = genome_pos, y = ratio.centered)) +
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
}

# =============================================================================
# SECTION: BODY
# =============================================================================

ratios <- read_tsv(ratios_tsv)
check_nonempty(ratios, "ratios")
log_n(ratios, "ratios read")

ratios <- add_genome_pos(ratios)
chr_bounds <- compute_chr_bounds(ratios)

n_libs <- length(unique(ratios$library))
plot_height <- max(4, 2 * n_libs)

p <- build_ratio_plot(ratios, chr_bounds)
ggsave(output_pdf, p, width = 12, height = plot_height)
