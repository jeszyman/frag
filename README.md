Fragmentomics pipeline for cell-free DNA whole-genome sequencing analysis. Extracts fragment length features, DELFI ratios, end-motif profiles, and NMF signatures from paired-end WGS data.  

`frag.org` is the source of truth: the workflows, scripts, configs, tests and tools in this repository are tangled from it. The tangled files are committed complete, so the pipeline is runnable as committed without Emacs; edit `frag.org` and re-tangle rather than editing a tangled file.  

![img](resources/figures/frag.dag.png)  


# Continuous Integration

[![test-data](https://img.shields.io/github/actions/workflow/status/jeszyman/frag/test-data.yaml?branch=master&label=test-data)](https://github.com/jeszyman/frag/actions/workflows/test-data.yaml)

[![smk-dry](https://img.shields.io/github/actions/workflow/status/jeszyman/frag/smk-dry.yaml?branch=master&label=smk-dry)](https://github.com/jeszyman/frag/actions/workflows/smk-dry.yaml)

[![smk-run](https://img.shields.io/github/actions/workflow/status/jeszyman/frag/smk-run.yaml?branch=master&label=smk-run)](https://github.com/jeszyman/frag/actions/workflows/smk-run.yaml)


# Including from your wrapper

`workflows/frag.smk` contains rules only and reads no config. Every path and parameter is a Python variable the including wrapper defines before its `include:` line; the list, one line per variable, is the comment at the top of `workflows/frag.smk`. `workflows/test.smk` is the complete example: copy it, keep its variable block, and replace `rule all` with the outputs your project needs.  


# Configuration

A run takes one YAML config passed with `--configfile`. `config/common.yaml` is the template for a new project (every key, with placeholder paths); `config/test.yaml` configures the committed chr22 fixture. `config/config.schema.yaml` declares every key the wrappers read, and the wrapper validates the loaded config against it before building the DAG, printing every violation and every undeclared key together. The sample sheet (`sample-tsv-path`) is a TSV with columns `r1_basename`, `r2_basename` (relative to `directories.inputs`), `library_id` and `cohort`; `config/samples.tsv` is the template and `config/samples_test.tsv` the fixture sheet.  

End motifs are the 4 reference bases at the 5′ end of every paired, primary, mapped, non-duplicate read with MAPQ ≥ 30, strand-aware (see End motifs in `frag.org`). `end_motif.max_ends` caps the ends counted per library by read-name-hash subsampling, which keeps mates together and samples the whole genome; `null` (the default) counts every end. `end_motif.seed` sets the subsample seed.  

The conda environment is `envs/frag-conda-env.yaml` (every package pinned); run with `--use-conda`.  


# Testing

-   Fixture run: `tests/full/inputs/` holds a chr22 reference subset and 60,000 read pairs from each of four PRJNA326698 cfDNA WGS libraries (the wrapper uses lib001 and lib002), rebuilt by `tools/get_test_data.sh`. From the repository root, `snakemake --configfile config/test.yaml --use-conda --cores 4` runs it (the root `Snakefile` includes `workflows/test.smk`); add `--dry-run` to plan only.
-   Pinned outputs: `tools/check_expected_outputs.sh tests/full` compares 32 deterministic fixture outputs to `tests/expected_outputs.md5`.
-   Unit tests: `python -m pytest -q tests/unit` (known-answer tests of the end-motif extractor and the motif matrix; the extractor tests skip without samtools and bedtools on PATH).
-   Continuous integration (GitHub Actions, on push to master and on pull requests): smk-dry (DAG dry run and unit tests), smk-run (full fixture run and pinned-output check), test-data (fixture rebuild).


# Change log

-   Development since last tag  
    -   <span class="timestamp-wrapper"><span class="timestamp">[2026-09-24 Thu] </span></span> End motifs are counted strand-aware from the reference at the 5′ end of every paired, primary, mapped, non-duplicate read with MAPQ ≥ 30 (`samtools view -q 30 -f 1 -F 3340`, both mates). The previous extraction keyed the end on read number (read 1's leftmost 4 bases, read 2's rightmost 4 bases reverse-complemented), so for every pair whose read 1 maps to the minus strand it returned the read's 3′ end, interior to the fragment whenever the fragment is longer than the read (about 40 % of motifs in Axel's and Peter's data); its subsample was truncated in coordinate order, and it used MAPQ 60 and kept duplicates. `scripts/frag_end_motifs.sh` and `scripts/frag_motif_matrix.py` are adapted from Axel's extractor (commit 7a5767ab1); rule `frag_sample_motifs` is now `frag_end_motifs`, per-library outputs are `motifs/{library_id}.{ref_name}.motif_counts.tsv` (256 motifs plus `OTHER`), and a motif x library count matrix is written beside `all_motifs.tsv`. Config `end_motif.max_ends` (null counts every end) and `end_motif.seed` replace `n_motif` and `n_reads`. This changes every motif, motif-diversity and F-profile output, and changes the module interface, so v2.0.0 is the intended next tag. On the chr22 fixture lib001 counts 3,068 ends (CC-end share 0.121) and lib002 2,979 (0.105). On the full-genome healthy library kh\_01 (Validation) every one of 192,998,342 qualifying reads is counted, CCCA, CCTG and CCAG lead the table (the old extraction had AAAA second), the CC-end share is 0.157 against 0.112 for the old extraction, and a 10 M `max_ends` subsample reproduces the all-ends fractions (r = 0.99998) with every chromosome sampled in proportion, where the old extraction sampled only chr1 to chr11. Verified: 10 known-answer unit tests pass; fixture run completes; kh\_01 runs exit 0.
    -   <span class="timestamp-wrapper"><span class="timestamp">[2026-09-24 Thu] </span></span> Repository brought to the current biopipe standard. The conda env moves to `envs/frag-conda-env.yaml` with every package pinned (snakemake 9.14.5, samtools 1.24, bedtools 2.31.1, plus jsonschema and pytest). `workflows/frag.smk` reads no config: the wrappers bind every value as a variable, listed in the module preamble. `config/config.schema.yaml` and `check_config()` validate the config at startup. `config/common.yaml` and `config/samples.tsv` are templates; the fixture sheet is `config/samples_test.tsv`; a one-line root `Snakefile` includes the test wrapper; fixture logs and benchmarks go under `tests/full/`. Reference data files move to `resources/data/`, the DAG to `resources/figures/frag.dag.png`; tracked `docs/` files are untracked; `LICENSE` (MIT) and `CITATION.cff` are added; `config/int_test.yaml` (read by nothing) and the copied reference implementation folder are removed. The full-genome kh\_01 wrapper and config get org blocks under Validation, which also records the kh\_01 end-motif comparison. 32 deterministic fixture outputs are pinned in `tests/expected_outputs.md5` and checked by `tools/check_expected_outputs.sh`. CI triggers on push to master and pull requests, and runs the unit tests and the pin check. This change log replaces `CHANGELOG.md`. Verified: unit tests pass; fixture run completes; every pin ok on a clean run.
-   <span class="timestamp-wrapper"><span class="timestamp">[2026-06-11 Thu] </span></span> v1.0.0. First tagged release: a biopipe-standard-conformant cfDNA fragmentomics pipeline.  
    -   Restructured to the biopipe module pattern; migrated the fragment-length (GC normalization chain) and end-motif pipelines from the archived `cfdna` repo into `workflows/frag.smk`.
    -   Migrated feature extraction from `nf1_fragmentomics`: fragment-length histograms, DELFI arm z-scores, ratio row-normalization, and QC plots (all R).
    -   Migrated from `candetect/delfi`: NMF fragment-length features, F-profiles (NMF + NNLS on motifs), motif-diversity score (Shannon entropy).
    -   Added GitHub Actions CI (test-data, smk-dry, smk-run) and the `frag_check_ids` identifier-compatibility gate.
    -   Brought the repo structure in line with the biopipe standard: CLAUDE.md / AGENTS.md / CHANGELOG.md, test sample sheet at `config/samples.tsv`, DAG figure at `docs/frag.dag.png`; removed legacy `test/`, the org LaTeX-preview cache, the `frag.md` stub, the orphaned `count_scale.R`, and the tangle preamble headers.
    -   Added the `frag_checks.{R,py}` data-integrity library and conformed all 13 R analysis scripts plus the 2 Python scripts to the code style guide, guarding read inputs.
    -   Seeded `sample_frags_by_gc` (`set.seed(42)`) for reproducibility; made `figure_schema.R` optional for CI with relative test-config paths.
    -   Full-genome validation: 6-sample run (3 healthy, 3 cancer) on ncbi\_hg38; all steps pass.
    -   Added the test data set (WGS cfDNA, PRJNA326698, chr22).

