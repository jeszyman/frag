if (!exists("freq_long"))    source("00_load_data.R")
if (!exists("our_perm_tau")) our_perm_tau <- read_tsv("data/our_perm_tau.tsv", show_col_types = FALSE)

N_BOOT  <- 1000
EPSILON <- 1e-8
plan(multisession, workers = parallel::detectCores() - 1)
set.seed(42)

run_boot <- function(i) {
  c_boot <- sample(cancer_ids,  length(cancer_ids),  replace = TRUE)
  h_boot <- sample(healthy_ids, length(healthy_ids), replace = TRUE)
  p_c <- freq_long |> filter(sample %in% c_boot) |>
    group_by(fragment_length) |> summarise(p_c = mean(frequency), .groups = "drop")
  p_h <- freq_long |> filter(sample %in% h_boot) |>
    group_by(fragment_length) |> summarise(p_h = mean(frequency), .groups = "drop")
  p_c |> left_join(p_h, by = "fragment_length") |>
    mutate(enrichment = p_c / (p_h + EPSILON)) |>
    left_join(our_perm_tau, by = "fragment_length") |>
    mutate(exceeds = enrichment > threshold_95) |>
    select(fragment_length, exceeds)
}

boot_results <- future_map_dfr(1:N_BOOT, run_boot,
                               .options = furrr_options(seed = TRUE),
                               .progress = TRUE)

sel_freq <- boot_results |>
  group_by(fragment_length) |>
  summarise(selection_freq = mean(exceeds), .groups = "drop")

cat("Max selection frequency:", round(max(sel_freq$selection_freq), 3), "\n")
cat("Bins >= 0.50:", sum(sel_freq$selection_freq >= 0.50), "\n")
cat("Bins >= 0.70:", sum(sel_freq$selection_freq >= 0.70), "\n")
cat("Bins >= 0.90:", sum(sel_freq$selection_freq >= 0.90), "\n")

dir.create("data", showWarnings = FALSE)
write_tsv(sel_freq, "data/02_standard_selection_freq.tsv")
