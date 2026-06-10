#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# plot_frag_length_overlay.R — fragment-length frequency overlay across libraries.
# Input : freq_csv    — library x length frequency matrix
#         output_pdf  — overlay plot (PDF)
#         samples_tsv — optional (library_id, group) for coloring by group
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

# Positional args (repo convention): Snakemake passes inputs then output.
args <- commandArgs(trailingOnly = TRUE)
freq_csv    <- args[1]
output_pdf  <- args[2]
samples_tsv <- if (length(args) >= 3) args[3] else NULL

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

#' Reshape the frequency matrix to long form (length <= 400), optionally joining groups.
#'
#' @param freq_mat Data frame (library rownames, one column per length).
#' @param samples_tsv Optional path to (library_id, group) TSV.
#' @return Tibble (library, length, freq[, group]).
build_freq_long <- function(freq_mat, samples_tsv) {
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
  }
  freq_long
}

#' Build the fragment-length overlay plot.
#'
#' @param freq_long Long-form frequency tibble.
#' @param color_aes Column name to color by ("group" or "library").
#' @return A ggplot object.
build_overlay_plot <- function(freq_long, color_aes) {
  ggplot(freq_long, aes(x = length, y = freq,
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
}

# =============================================================================
# SECTION: BODY
# =============================================================================

freq_mat <- read.csv(freq_csv, row.names = 1, check.names = FALSE)
check_nonempty(freq_mat, "freq_mat")
log_n(freq_mat, "freq matrix")

color_aes <- if (!is.null(samples_tsv) && file.exists(samples_tsv)) "group" else "library"
freq_long <- build_freq_long(freq_mat, samples_tsv)

p <- build_overlay_plot(freq_long, color_aes)
ggsave(output_pdf, p, width = PLOT_W, height = PLOT_H)
