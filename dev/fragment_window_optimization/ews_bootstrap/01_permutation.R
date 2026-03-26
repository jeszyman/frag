if (!exists("freq_long")) source("00_load_data.R")

N_PERM  <- 5000
EPSILON <- 1e-8

all_samples <- c(cancer_ids, healthy_ids)
n_cancer    <- length(cancer_ids)

plan(multisession, workers = parallel::detectCores() - 1)
set.seed(42)

perm_results <- future_map_dfr(1:N_PERM, function(i) {
  perm_c <- sample(all_samples, n_cancer, replace = FALSE)
  perm_h <- setdiff(all_samples, perm_c)
  p_c <- freq_long |> filter(sample %in% perm_c) |>
    group_by(fragment_length) |> summarise(p_c = mean(frequency), .groups = "drop")
  p_h <- freq_long |> filter(sample %in% perm_h) |>
    group_by(fragment_length) |> summarise(p_h = mean(frequency), .groups = "drop")
  p_c |> left_join(p_h, by = "fragment_length") |>
    mutate(enrichment_null = p_c / (p_h + EPSILON), perm_id = i) |>
    select(perm_id, fragment_length, enrichment_null)
}, .options = furrr_options(seed = TRUE), .progress = TRUE)

our_perm_tau <- perm_results |>
  group_by(fragment_length) |>
  summarise(threshold_95 = quantile(enrichment_null, 0.95),
            threshold_99 = quantile(enrichment_null, 0.99),
            .groups = "drop")

obs_enrichment <- freq_long |>
  group_by(fragment_length, selected) |>
  summarise(mean_freq = mean(frequency), .groups = "drop") |>
  pivot_wider(names_from = selected, values_from = mean_freq) |>
  mutate(enrichment = cancer / (healthy + EPSILON)) |>
  left_join(our_perm_tau, by = "fragment_length") |>
  mutate(sig_95 = enrichment > threshold_95,
         sig_99 = enrichment > threshold_99)

cat("Bins exceeding 95th percentile:", sum(obs_enrichment$sig_95),
    "| range:", range(filter(obs_enrichment, sig_95)$fragment_length), "\n")
cat("Bins exceeding 99th percentile:", sum(obs_enrichment$sig_99),
    "| range:", range(filter(obs_enrichment, sig_99)$fragment_length), "\n")

dir.create("data", showWarnings = FALSE)
write_tsv(our_perm_tau,   "data/our_perm_tau.tsv")
write_tsv(obs_enrichment, "data/our_obs_enrichment.tsv")
# Save raw null at representative bins for stability diagnostic
perm_results |>
  filter(fragment_length %in% c(70, 80, 90, 100, 110, 120)) |>
  write_tsv("data/perm_null_representative.tsv")
plan(sequential)
