#!/usr/bin/env Rscript

# -----------------------------------------------------------------------------
# end_motif_mat.R — per-library 4-mer end-motif fraction matrix.
# Input : motif_str — space-separated paths to per-library motif files
# Output: motif_tsv — motif x library fraction matrix (all 256 4-mers)
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
motif_str <- args[1]
motif_tsv <- args[2]

# =============================================================================
# SECTION: WORK FUNCTIONS
# =============================================================================

#' All 256 possible 4-mers, sorted.
#'
#' @return Tibble (motif).
build_possible_motifs <- function() {
  expand.grid(rep(list(c("A", "G", "T", "C")), 4)) %>%
    as_tibble() %>%
    mutate(motif = paste0(Var1, Var2, Var3, Var4)) %>%
    select(motif) %>%
    arrange(motif)
}

#' Read one library's motif calls and return per-motif fraction.
#'
#' @param motif_file Path to a motif file (one motif per line).
#' @return Tibble (motif, fract).
ingest_motif <- function(motif_file) {
  read_tsv(motif_file, col_names = c("motif")) %>%
    check_nonempty("motifs") %>%
    group_by(motif) %>%
    summarise(count = n()) %>%
    mutate(fract = count / sum(count)) %>%
    select(motif, fract)
}

# =============================================================================
# SECTION: BODY
# =============================================================================

possible_motifs <- build_possible_motifs()

motif_files <- strsplit(motif_str, " ")[[1]]
# Extract library_id from {library_id}.{ref_name}.motifs.txt
names(motif_files) <- basename(motif_files) %>%
  sub("\\.[^.]+\\.motifs\\.txt$", "", .)

motif_tibs <- lapply(motif_files, ingest_motif)

motifs <- bind_rows(motif_tibs, .id = "library") %>%
  pivot_wider(names_from = library, values_from = fract) %>%
  filter(motif %in% possible_motifs$motif)
log_n(motifs, "motif matrix")

write_tsv(motifs, motif_tsv)
