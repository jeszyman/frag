library(tidyverse)
library(furrr)

DATA_DIR <- "/mnt/gcs/lilli/ewing/analysis_bundle/permutation_bootstrap_jeff"

freq_wide        <- read_tsv(file.path(DATA_DIR, "per_bam_fragment_frequencies.tsv"),    show_col_types = FALSE)
lilli_enrichment <- read_tsv(file.path(DATA_DIR, "fragment_enrichment_curve.tsv"),       show_col_types = FALSE)
lilli_perm_tau   <- read_tsv(file.path(DATA_DIR, "per_base_permutation_thresholds.tsv"), show_col_types = FALSE)
meta             <- read_tsv(file.path(DATA_DIR, "sample_selection_metadata.tsv"),       show_col_types = FALSE)
sig_bins         <- read_tsv(file.path(DATA_DIR, "significant_fragment_bins.tsv"),       show_col_types = FALSE)

freq_long <- freq_wide |>
  pivot_longer(-sample, names_to = "length_str", values_to = "frequency") |>
  mutate(fragment_length = as.integer(str_remove(length_str, "bp_"))) |>
  select(-length_str) |>
  left_join(meta |> select(sample, selected, tf), by = "sample")

cancer_ids  <- meta |> filter(selected == "cancer")  |> pull(sample)
healthy_ids <- meta |> filter(selected == "healthy") |> pull(sample)

cat("Cancer samples: ", length(cancer_ids), "\n")
cat("Healthy samples:", length(healthy_ids), "\n")
cat("Length bins:    ", n_distinct(freq_long$fragment_length), "\n")
