# Frag Repository Restructure — Design Spec

## Goal

Restructure the frag repo and frag.org to follow the emseq biopipe module pattern, migrating fragment length and end motif content from the cfdna repo. The result is a clean, self-contained Snakemake module for cfDNA fragmentomics analysis.

## Scope

### In Scope
- Rebuild frag.org pipeline sections using emseq.org as structural template
- Migrate and refactor fragment length pipeline from cfdna.org
- Migrate and refactor end motif pipeline from cfdna.org
- Refactor read processing rules (merge best of frag + cfdna)
- Rebuild `workflows/frag.smk` as a proper module snakefile
- Rebuild `workflows/test.smk` as a test wrapper
- Refactor test harness (`get_test_data.sh`, `config/test.yaml`)
- Search for open-access cfDNA WGS test data; fall back to EM-seq if needed
- Clean up legacy `workflow/` directory
- All scripts as org-babel tangle outputs

### Out of Scope (Do Not Touch)
- `** Fragment Length Window Optimization :bioinfo-dev:` section in frag.org
- `dev/` directory
- `data/` analysis outputs (windows_*.tsv, enrichment files, etc.)
- `resources/` reference materials (except adding end motif deps if needed)
- `.dir-locals.el`
- CNA pipeline (separate repo)
- Other cfdna.org content not related to fragment length or end motifs

## Structural Template: emseq.org

The emseq repo defines the biopipe module pattern:

1. **Single org file as source of truth** — all configs, workflows, scripts, tests tangled from org-babel blocks
2. **Section-level `:header-args:snakemake:`** on the main processing section tangles all nested snakemake blocks to `workflows/<module>.smk`
3. **Module snakefile** — prefixed rules, no `rule all`, no constants, references wrapper-provided variables
4. **Test wrapper** — `test.smk` defines `rule all`, loads config, includes the module
5. **Conda envs** — each as a separate org-babel YAML block, tangled to `config/`
6. **Scripts** — each in its own org-babel block with `:tangle ./scripts/<name>`
7. **Test harness** — `get_test_data.sh` fetches minimal test inputs, `config/test.yaml` defines test parameters

## Org File Structure

```
* Cell-free DNA Fragmentomics Analysis                              :biopipe:
:PROPERTIES:
:ID: <existing-or-new-uuid>
:END:

  ** README
  :PROPERTIES:
  :ID: 339b69f6-6c09-4e9d-a2ec-27bdf5747163
  :END:

  ** Repository setup and administration
     *** Org update
         - org-babel block → tools/shell/org_update.sh
     *** Get in-repo test data
         - org-babel block → tools/get_test_data.sh
     *** README
         - Export machinery (emacs_export_header_to_markdown.py reference)
     *** Frag conda environment
         - org-babel YAML block → config/frag-conda-env.yaml
     *** End motif conda environment
         - org-babel YAML block → config/frag-end-motif-conda-env.yaml
         - (only if deps differ from main env)

  ** Fragmentomics sequence processing
  :PROPERTIES:
  :header-args:snakemake: :tangle ./workflows/frag.smk
  :ID: <uuid>
  :END:

     *** Preamble
         - Modular snakefile comment block
         - No constants (provided by wrapper)

     *** Read processing
         **** FASTQ QC
              - snakemake: rule frag_fastp
              - bash block → scripts/fastp_wrapper.sh (or inline if simple)
         **** BWA index
              - snakemake: rule frag_bwa_index
         **** Alignment
              - snakemake: rule frag_align
              - bash block → scripts/bwa_mem_markdup_stream.sh
         **** Filter alignments
              - snakemake: rule frag_filter_alignments
              - bash block → scripts/filter_alignments.sh

     *** Fragment length
         **** BAM to fragment BED
              - snakemake: rule frag_bam_to_frag_bed
              - bash block → scripts/bam_to_frag_bed.sh
         **** GC distribution
              - snakemake: rule frag_gc_distro
              - R block → scripts/gc_distro.R
         **** Healthy GC summary
              - snakemake: rule frag_healthy_gc
              - R block → scripts/make_healthy_gc_summary.R
         **** Sample by GC
              - snakemake: rule frag_sample_by_gc
              - R block → scripts/sample_frags_by_gc.R
         **** Window sum
              - snakemake: rule frag_window_sum
              - bash block → scripts/frag_window_sum.sh
         **** Window count
              - snakemake: rule frag_window_count
              - bash block → scripts/frag_window_int.sh
         **** Count merge
              - snakemake: rule frag_count_merge
              - bash block → scripts/count_merge.sh
         **** Ratio normalization
              - snakemake: rule frag_ratio_normalize
              - R block → scripts/make_ratios.R

     *** End motifs
         **** Sample motifs
              - snakemake: rule frag_sample_motifs
              - script block → scripts/sample_motifs.sh (or .py)
         **** Motif matrix
              - snakemake: rule frag_motif_matrix
              - R block → scripts/end_motif_mat.R

  ** Testing
     *** Full test
         - snakemake block → workflows/test.smk
         - Defines rule all, loads config, includes frag.smk
         - resolve_config_paths() preamble
         - Constants: D_FRAG, D_LOGS, CONDA_FRAG, etc.

  ** Fragment Length Window Optimization :bioinfo-dev:
     ← UNTOUCHED — existing content preserved exactly

  ** Ideas
```

## Snakemake Architecture

### Module: `workflows/frag.smk`

- All rules prefixed `frag_`
- Style guide directive order: message, conda, input, log, benchmark, params, threads, output
- Each rule references external scripts, not inline shell logic
- Conda env via variable (e.g., `CONDA_FRAG`), defined by wrapper
- Paths use f-strings with `D_FRAG` directory variable
- No `rule all`, no constant definitions
- Modular comment header block

### Wrapper: `workflows/test.smk`

- `resolve_config_paths()` function
- Load config, define constants (`D_FRAG`, `D_LOGS`, `D_BENCHMARK`, `CONDA_FRAG`)
- Load sample TSV (`data/test-samples.tsv`)
- `rule all:` with explicit expand patterns covering all test outputs
- `include: "frag.smk"`

### Rule Summary

| Rule | Category | Input | Output |
|------|----------|-------|--------|
| `frag_fastp` | Read processing | Raw FASTQs | Trimmed FASTQs + reports |
| `frag_bwa_index` | Read processing | Reference FASTA | BWA index files |
| `frag_align` | Read processing | Trimmed FASTQs + index | Sorted, deduped BAM |
| `frag_filter_alignments` | Read processing | BAM + keep/exclude BEDs | Filtered BAM |
| `frag_bam_to_frag_bed` | Fragment length | Filtered BAM | Fragment BED (chr, start, end, gc, length) |
| `frag_gc_distro` | Fragment length | Fragment BED | Per-sample GC distribution |
| `frag_healthy_gc` | Fragment length | Multiple GC distros | Healthy reference GC profile |
| `frag_sample_by_gc` | Fragment length | Fragment BED + healthy GC | GC-resampled fragment BED |
| `frag_window_sum` | Fragment length | Fragment BED | Short/long partitioned counts |
| `frag_window_count` | Fragment length | Partitioned BED + genomic windows | Per-window intersection counts |
| `frag_count_merge` | Fragment length | Per-library counts | Merged count matrix |
| `frag_ratio_normalize` | Fragment length | Merged counts | Zero-centered, unit-SD ratios |
| `frag_sample_motifs` | End motifs | Filtered BAM | Per-sample 5' motif frequencies |
| `frag_motif_matrix` | End motifs | Multiple motif files | Motif frequency matrix |

## Configuration

### `config/test.yaml`

Nested YAML following style guide:

```yaml
conda:
  frag: "${HOME}/repos/frag/config/frag-conda-env.yaml"
  end_motif: "${HOME}/repos/frag/config/frag-end-motif-conda-env.yaml"

directories:
  inputs: tests/full/inputs
  frag: tests/full
  logs: tests/logs
  benchmark: tests/benchmark

sample-tsv-path: ${HOME}/repos/frag/data/test-samples.tsv
threads: 4

frag_ref_assemblies:
  chr22:
    url: <ncbi-url>
    input: chr22-test.fa.gz
    name: chr22

blklist: hg38-blacklist.v2.bed.gz
```

### Conda Environments

**`config/frag-conda-env.yaml`** — core tools:
- bwa, samtools, bedtools, fastp, fastqc, snakemake, r-tidyverse

**`config/frag-end-motif-conda-env.yaml`** — end motif tools (if needed):
- python, pysam, numpy, r-tidyverse

## Test Harness

### `tools/get_test_data.sh`
- Best-effort search for open-access cfDNA WGS on SRA/ENA (small, paired-end, hg38)
- **Default fallback: EM-seq test data** — mechanically identical (FASTQ → BAM → fragment BED) even though library prep differs. EM-seq reads exercise all pipeline rules; fragment length distributions will differ from production cfDNA but this is acceptable for CI validation of rule connectivity and script correctness.
- Download chr22 reference subset
- Download hg38 blacklist
- Generate keep/exclude BEDs

### Test Data
- 4 test libraries (2 healthy + 2 cancer preferred; EM-seq fallback: use emseq test libs)
- chr22 reference subset
- Small enough for CI (~60k reads per sample)

### Success Criterion
- `snakemake -s workflows/test.smk --dry-run` passes with no errors

## Cleanup

- Delete `workflow/` directory entirely (legacy)
- Remove any orphaned scripts not referenced by org-babel blocks
- Ensure `scripts/` contains only tangle outputs

## Migration Strategy

The agent has broad latitude to restructure scripts and rules as needed. Key principles:
- Org-babel blocks are the source of truth; scripts are tangle outputs
- **Fragment length pipeline**: cfdna.org is the authoritative source for GC normalization logic (gc_distro, healthy_gc, sample_by_gc, window_sum, window_count, count_merge, ratio_normalize). These rules implement a specific statistical pipeline: per-fragment GC annotation → healthy reference GC profile → GC-stratified resampling → short/long partitioning (100-150 vs 151-220 bp) → 5Mb window intersection → count merge → zero-centered unit-SD ratio normalization. Preserve this logic faithfully.
- **Read processing**: frag's current implementation (frag_fastp, frag_align with streaming markdup, frag_filter_alignments) is more recent and already refactored — prefer frag's versions, augmenting from cfdna where cfdna has additional functionality.
- **End motifs**: cfdna.org is the sole source. Refactor for style guide compliance.
- Merge the best of both where they overlap
- Follow the Snakemake style guide (basecamp.org line 9457) rigorously

## Configuration Requirements

All keys below are **required** in any config YAML that includes frag.smk:

```yaml
# Required keys
conda:
  frag: <path>            # Path to frag-conda-env.yaml
directories:
  inputs: <path>          # Input data directory
  frag: <path>            # Main output directory
  logs: <path>            # Log directory
  benchmark: <path>       # Benchmark directory
sample-tsv-path: <path>   # TSV with columns: library_id, r1, r2, cohort
threads: <int>            # Default thread count
frag_ref_assemblies:
  <name>:
    url: <url>            # Reference FASTA download URL
    input: <filename>     # Local filename for reference
    name: <name>          # Assembly name
blklist: <path>           # ENCODE blacklist BED

# Required for end motifs (optional if end motif rules not in rule all)
conda:
  end_motif: <path>       # Path to end-motif-conda-env.yaml (if separate)
```
