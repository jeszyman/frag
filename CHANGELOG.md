# Changelog

All notable changes to frag. Format loosely follows Keep a Changelog. The repo
is versioned as a whole with SemVer git tags (`vX.Y.Z`).

## [Unreleased]
- Add `frag_check_ids` identifier-compatibility gate.
- Remove tangle preamble headers from all scripts.
- Make `figure_schema.R` optional for CI; use relative paths in the test config
  for CI compatibility.
- Bring repo structure in line with the biopipe standard: add CLAUDE.md /
  AGENTS.md / CHANGELOG.md, migrate the test sample sheet to `config/samples.tsv`,
  relocate the DAG figure to `docs/frag.dag.png`, remove legacy `test/`, the
  org LaTeX-preview cache, and the `frag.md` stub.
- Add the `frag_checks.{R,py}` data-integrity library and conform all 13 R
  analysis scripts plus the 2 Python scripts to the code style guide, guarding
  read inputs; remove the orphaned `count_scale.R`.
- Seed `sample_frags_by_gc` (`set.seed(42)`) so the GC-weighted resample is
  reproducible.

## [1.0.0] — proposed first release tag
- Restructured to the biopipe module pattern; migrated the fragment-length (GC
  normalization chain) and end-motif pipelines from the archived `cfdna` repo
  into `workflows/frag.smk`.
- Migrated feature extraction from `nf1_fragmentomics`: fragment-length
  histograms, DELFI arm z-scores, ratio row-normalization, and QC plots (all R).
- Migrated from `candetect/delfi`: NMF fragment-length features, F-profiles
  (NMF + NNLS on motifs), motif-diversity score (Shannon entropy).
- Added GitHub Actions CI (test-data, smk-dry, smk-run).
- Full-genome validation: 6-sample run (3 healthy, 3 cancer) on ncbi_hg38 — all
  steps pass.
- Added the test data set (WGS cfDNA, PRJNA326698, chr22).

<!-- To cut the release: finalize the [Unreleased] section, then:
     git tag v1.0.0 && git push && git push --tags -->
