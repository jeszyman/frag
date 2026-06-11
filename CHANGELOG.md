# Changelog

All notable changes to frag. Format loosely follows Keep a Changelog. The repo
is versioned as a whole with SemVer git tags (`vX.Y.Z`).

## [1.0.0] - 2026-06-11

First tagged release: a biopipe-standard-conformant cfDNA fragmentomics pipeline.

- Restructured to the biopipe module pattern; migrated the fragment-length (GC
  normalization chain) and end-motif pipelines from the archived `cfdna` repo
  into `workflows/frag.smk`.
- Migrated feature extraction from `nf1_fragmentomics`: fragment-length
  histograms, DELFI arm z-scores, ratio row-normalization, and QC plots (all R).
- Migrated from `candetect/delfi`: NMF fragment-length features, F-profiles
  (NMF + NNLS on motifs), motif-diversity score (Shannon entropy).
- Added GitHub Actions CI (test-data, smk-dry, smk-run) and the
  `frag_check_ids` identifier-compatibility gate.
- Brought the repo structure in line with the biopipe standard: CLAUDE.md /
  AGENTS.md / CHANGELOG.md, test sample sheet at `config/samples.tsv`, DAG figure
  at `docs/frag.dag.png`; removed legacy `test/`, the org LaTeX-preview cache, the
  `frag.md` stub, the orphaned `count_scale.R`, and the tangle preamble headers.
- Added the `frag_checks.{R,py}` data-integrity library and conformed all 13 R
  analysis scripts plus the 2 Python scripts to the code style guide, guarding
  read inputs.
- Seeded `sample_frags_by_gc` (`set.seed(42)`) for reproducibility; made
  `figure_schema.R` optional for CI with relative test-config paths.
- Full-genome validation: 6-sample run (3 healthy, 3 cancer) on ncbi_hg38 — all
  steps pass.
- Added the test data set (WGS cfDNA, PRJNA326698, chr22).
