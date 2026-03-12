# Frag Repository Restructure Implementation Plan

> **For agentic workers:** REQUIRED: Use superpowers:subagent-driven-development (if subagents available) or superpowers:executing-plans to implement this plan. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Restructure frag.org and the frag repo to follow the emseq biopipe module pattern, migrating fragment length and end motif content from cfdna.

**Architecture:** Single org file (frag.org) as source of truth. All configs, scripts, and workflows tangled from org-babel blocks. Module snakefile (`workflows/frag.smk`) with `frag_` prefixed rules included by a test wrapper (`workflows/test.smk`). Fragment length pipeline from cfdna + read processing from frag + end motifs from cfdna, all refactored to style guide conventions.

**Tech Stack:** Snakemake, org-babel, bash, R, conda, bedtools, samtools, bwa

**Spec:** `docs/superpowers/specs/2026-03-12-frag-restructure-design.md`

**Key reference files (READ THESE FIRST):**
- Structural template: `/home/jeszyman/repos/emseq/emseq.org`, `/home/jeszyman/repos/emseq/workflows/test.smk`
- Content source (fragment length): `/home/jeszyman/repos/cfdna/workflows/frag.smk`, `/home/jeszyman/repos/cfdna/cfdna.org`
- Content source (end motifs): `/home/jeszyman/repos/cfdna/scripts/sample_motifs.sh`, `/home/jeszyman/repos/cfdna/scripts/end_motif_mat.R`
- Current frag state: `/home/jeszyman/repos/frag/frag.org`, `/home/jeszyman/repos/frag/workflows/frag.smk`, `/home/jeszyman/repos/frag/workflows/test.smk`
- Snakemake style guide: `/home/jeszyman/repos/basecamp/basecamp.org` lines 9457–9646

**CRITICAL CONSTRAINTS:**
- DO NOT touch the `** Fragment Length Window Optimization` section in frag.org (bioinfo-dev content)
- DO NOT touch `dev/`, `data/` analysis outputs, `.dir-locals.el`
- All edits to scripts/configs/workflows MUST be made in org-babel blocks in frag.org — the files are tangle outputs
- After editing frag.org, tangle it: `emacsclient --socket-name ~/.emacs.d/server/server -e '(progn (find-file "/home/jeszyman/repos/frag/frag.org") (org-babel-tangle))'`
- Rule prefix: `frag_` for ALL rules
- Style guide directive order: message, conda, input, log, benchmark, params, threads, output

---

## Chunk 1: Repository Setup and Administration

### Task 1: Restructure frag.org top-level headings

**Files:**
- Modify: `/home/jeszyman/repos/frag/frag.org`

**Context:** The current frag.org has the window optimization section near the top. We need to reorganize to match the emseq pattern: README → Repository setup → Conda envs → Processing → Testing → (bioinfo-dev stays) → Ideas.

Read the current frag.org to understand full structure before making changes. Preserve ALL existing content — just reorganize headings.

- [ ] **Step 1: Read frag.org completely**

Read the entire file to identify all sections and their boundaries.

- [ ] **Step 2: Reorganize top-level headings**

The target structure (each `**` is a second-level heading under the root `*`):

```org
* Cell-free DNA Fragmentomics Analysis                              :biopipe:
:PROPERTIES:
:ID: <keep-existing-or-generate>
:END:

** README
** Repository setup and administration
** Conda environments
** Fragmentomics sequence processing
** Testing
** Fragment Length Window Optimization :bioinfo-dev:  ← EXISTING, UNTOUCHED
** Development  ← EXISTING, if present
** Ideas  ← EXISTING, if present
```

Move existing content to match. The `** Fragment Length Window Optimization` section and everything below it stays exactly as-is.

- [ ] **Step 3: Commit**

```bash
git add frag.org
git commit -m "refactor: reorganize frag.org top-level headings to match biopipe pattern"
```

### Task 2: Repository setup section — org update and get test data

**Files:**
- Modify: `/home/jeszyman/repos/frag/frag.org` (Repository setup section)

**Context:** The current frag.org already has org_update.sh and get_test_data.sh content. Keep these — they're already good. Ensure they're properly placed under `** Repository setup and administration` with the right subsection structure.

- [ ] **Step 1: Verify org update block**

Ensure `*** Org update` subsection exists under `** Repository setup and administration` with a bash tangle block targeting `./tools/shell/org_update.sh`. The current content is fine — just verify placement.

- [ ] **Step 2: Verify get test data block**

Ensure `*** Get in-repo test data` subsection exists with a bash tangle block targeting `./tools/get_test_data.sh`. Current content is fine.

- [ ] **Step 3: Ensure README subsection exists**

Ensure `*** README` subsection exists with the export machinery reference (emacs_export_header_to_markdown.py). Current content is fine if present.

- [ ] **Step 4: Commit**

```bash
git add frag.org
git commit -m "refactor: verify repository setup section structure"
```

### Task 3: Conda environment blocks

**Files:**
- Modify: `/home/jeszyman/repos/frag/frag.org`

**Context:** Create a `** Conda environments` section with YAML tangle blocks. The main env needs: fastp, fastqc, bedtools, samtools, bwa, snakemake, r-tidyverse. End motif env may need: python, pysam if not already in main.

- [ ] **Step 1: Create conda environments section**

Under `** Conda environments`, add:

```org
*** Frag
#+begin_src yaml :tangle ./config/frag-conda-env.yaml
name: frag

channels:
  - conda-forge
  - bioconda

dependencies:
  - bedtools
  - bwa
  - fastp
  - fastqc
  - r-tidyverse
  - samtools
  - snakemake
#+end_src
```

Note: added `r-tidyverse` (needed for GC distro, ratios, end motif R scripts). Deps sorted alphabetically.

- [ ] **Step 2: Add end motif conda env if needed**

Check if pysam/python are needed beyond what the main env provides. If end motif scripts only use bash+samtools+bedtools (which `sample_motifs.sh` does) and R (which `end_motif_mat.R` does), the main env suffices. If so, skip this step.

If a separate env is needed:
```org
*** End motif
#+begin_src yaml :tangle ./config/frag-end-motif-conda-env.yaml
name: frag-end-motif

channels:
  - conda-forge
  - bioconda

dependencies:
  - bedtools
  - r-tidyverse
  - samtools
  - snakemake
#+end_src
```

- [ ] **Step 3: Tangle and commit**

```bash
emacsclient --socket-name ~/.emacs.d/server/server -e '(progn (find-file "/home/jeszyman/repos/frag/frag.org") (org-babel-tangle))'
git add frag.org config/frag-conda-env.yaml
git commit -m "refactor: add conda environment blocks to frag.org"
```

---

## Chunk 2: Fragmentomics Sequence Processing — Read Processing Rules

### Task 4: Create processing section with preamble and read processing rules

**Files:**
- Modify: `/home/jeszyman/repos/frag/frag.org`

**Context:** Create `** Fragmentomics sequence processing` with section-level `:header-args:snakemake: :tangle ./workflows/frag.smk`. This makes all nested snakemake blocks tangle to `workflows/frag.smk`.

Read processing rules come from frag's EXISTING implementations (frag_fastp, frag_bwa_index, frag_align, frag_filter_alignments, frag_bam_to_frag_bed). These are already refactored and good — keep them.

- [ ] **Step 1: Create the processing section header**

```org
** Fragmentomics sequence processing
:PROPERTIES:
:header-args:snakemake: :tangle ./workflows/frag.smk
:ID: <generate-uuid>
:END:
```

- [ ] **Step 2: Add preamble subsection**

```org
*** Preamble

#+begin_src snakemake
#########1#########2#########3#########4#########5#########6#########7#########8
#
# This is a modular snakefile, intended to be incorporated into a larger
# workflow using the "include:" directive. (See
# https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html)
#
#########1#########2#########3#########4#########5#########6#########7#########8
#+end_src

- [ ] **Step 3: Add read processing subsections**

Under `*** Read processing`, create four subsections — each containing a snakemake rule block and (where applicable) a script tangle block.

Keep the existing frag rules: `frag_fastp`, `frag_bwa_index`, `frag_align`, `frag_filter_alignments`, `frag_bam_to_frag_bed`.

Refactor each rule to match style guide directive order: message, conda, input, log, benchmark, params, threads, output.

The current frag.smk rules have some ordering issues (e.g., `frag_fastp` puts threads before params, output before log). Fix these.

Example corrected rule for frag_fastp:

```snakemake
rule frag_fastp:
    message:
        "Fragmentomics fastp FASTQ processing"
    conda:
        CONDA_FRAG
    input:
        r1 = f"{D_FRAG}/fastqs/{{library_id}}.raw_R1.fastq.gz",
        r2 = f"{D_FRAG}/fastqs/{{library_id}}.raw_R2.fastq.gz",
    log:
        cmd = f"{D_LOGS}/{{library_id}}_frag_fastp.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_frag_fastp.tsv"
    params:
        extra = config.get("fastp", {}).get("extra", ""),
    threads:
        8
    output:
        failed = f"{D_FRAG}/fastqs/{{library_id}}.failed.fastq.gz",
        html = f"{D_FRAG}/qc/{{library_id}}_frag_fastp.html",
        json = f"{D_FRAG}/qc/{{library_id}}_frag_fastp.json",
        r1 = f"{D_FRAG}/fastqs/{{library_id}}.processed_R1.fastq.gz",
        r2 = f"{D_FRAG}/fastqs/{{library_id}}.processed_R2.fastq.gz",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[fastp] $(date) lib={wildcards.library_id} threads={threads}"

        fastp \
          --detect_adapter_for_pe \
          --disable_quality_filtering \
          --in1 "{input.r1}" --in2 "{input.r2}" \
          --out1 "{output.r1}" --out2 "{output.r2}" \
          --failed_out "{output.failed}" \
          --json "{output.json}" --html "{output.html}" \
          --thread {threads} \
          {params.extra}
        """
```

Apply same directive reordering to frag_bwa_index, frag_align, frag_filter_alignments, frag_bam_to_frag_bed. Keep all the existing shell logic and script references.

Each rule's script block (bash tangle) should also be in the same subsection:
- `**** Alignment` contains both the snakemake rule block AND the `#+begin_src bash :tangle ./scripts/bwa_mem_markdup_stream.sh` block
- `**** Filter alignments` contains rule + `#+begin_src bash :tangle ./scripts/filter_alignments.sh`
- `**** BAM to fragment BED` contains rule + `#+begin_src bash :tangle ./scripts/bam_to_frag_bed.sh`

- [ ] **Step 4: Add frag_fastqc rule**

Keep the existing `frag_fastqc` rule, fix directive ordering.

- [ ] **Step 5: Tangle and verify**

```bash
emacsclient --socket-name ~/.emacs.d/server/server -e '(progn (find-file "/home/jeszyman/repos/frag/frag.org") (org-babel-tangle))'
```

Verify `workflows/frag.smk` was generated and contains all read processing rules.

- [ ] **Step 6: Commit**

```bash
git add frag.org workflows/frag.smk scripts/bwa_mem_markdup_stream.sh scripts/filter_alignments.sh scripts/bam_to_frag_bed.sh
git commit -m "refactor: add read processing rules to frag.org with style guide compliance"
```

---

## Chunk 3: Fragment Length Pipeline Rules

### Task 5: Add fragment length rules from cfdna

**Files:**
- Modify: `/home/jeszyman/repos/frag/frag.org`
- Reference: `/home/jeszyman/repos/cfdna/workflows/frag.smk` (authoritative source for fragment length pipeline logic)
- Reference: `/home/jeszyman/repos/cfdna/scripts/gc_distro.R`, `healthy_gc.R`, `gc_sample.R`, `frag_window_sum.sh`, `frag_window_int.sh`, `count_merge.sh`, `make_ratios.R`, `make_gc_map_bind.sh`

**Context:** The fragment length pipeline from cfdna implements: GC-filtered genomic bins → fragment BED → GC distribution → healthy GC reference → GC-stratified resampling → short/long partitioning → window intersection → count merge → ratio normalization. Migrate ALL of these, refactoring each to style guide conventions and `frag_` prefix.

The current frag.org already has `gc_distro` (unprefixed) — replace it with the full pipeline.

- [ ] **Step 1: Add GC/mappability bin rule**

Under `*** Fragment length`, add `**** GC and mappability restricted bins`:

Snakemake rule `frag_gc_map_bins` (refactored from cfdna's `make_gc_map_bind`):
- message, conda (CONDA_FRAG), input (gc5mb from config, blklist from config), log, benchmark, params (script path), threads (1), output (keep_5mb.bed)
- Style guide compliant
- Script tangle block: `#+begin_src bash :tangle ./scripts/make_gc_map_bind.sh`

Script content (from cfdna, add proper header/arg parsing):
```bash
#!/usr/bin/env bash
set -euo pipefail

gc5mb="$1"
blklist="$2"
keep="$3"

bedtools intersect -a "$gc5mb" -b "$blklist" -v -wa \
  | grep -v _ \
  | awk '{ if ($4 >= 0.3) print $0 }' > "$keep"
```

- [ ] **Step 2: Refactor gc_distro rule**

Replace current unprefixed `gc_distro` with `frag_gc_distro`:
- Per-library (not per-group) to match cfdna pattern
- Input: single fragment BED
- Output: per-library gc_distro.csv
- Script: `scripts/gc_distro.R` (from cfdna, refactored)

Script content (from cfdna's gc_distro.R — works on single BED file):
```r
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
bed_file <- args[1]
distro_file <- args[2]

library(tidyverse)

bed <- read.table(bed_file, sep = "\t")
names(bed) <- c("chr", "start", "end", "gc_raw", "len")

distro <- bed %>%
  mutate(gc_strata = round(gc_raw, 2)) %>%
  count(gc_strata) %>%
  mutate(fract_frags = n / sum(n)) %>%
  mutate(library_id = gsub("_frag.bed", "", gsub("^.*lib", "lib", bed_file))) %>%
  select(library_id, gc_strata, fract_frags)

write.csv(distro, file = distro_file, row.names = FALSE)
```

- [ ] **Step 3: Add healthy GC summary rule**

`frag_healthy_gc`:
- Input: expand gc_distro CSVs for healthy libraries
- Output: healthy_med.rds
- Script: `scripts/make_healthy_gc_summary.R` (from cfdna's healthy_gc.R)

- [ ] **Step 4: Add GC sample rule**

`frag_gc_sample`:
- Input: fragment BED + healthy_med.rds
- Output: sampled fragment BED
- Script: `scripts/sample_frags_by_gc.R` (from cfdna's gc_sample.R)

- [ ] **Step 5: Add window sum rule**

`frag_window_sum`:
- Input: sampled fragment BED
- Output: short BED + long BED (100-150 bp vs 151-220 bp)
- Script: `scripts/frag_window_sum.sh` (from cfdna)

Refactored script:
```bash
#!/usr/bin/env bash
set -euo pipefail

input_frag="$1"
output_short="$2"
output_long="$3"

awk '{ if ($4 >= 100 && $4 <= 150) print $0 }' "$input_frag" > "$output_short"
awk '{ if ($4 >= 151 && $4 <= 220) print $0 }' "$input_frag" > "$output_long"
```

Note: cfdna version uses `$4` and `$5` but fragment BED has length in column 5 (1-indexed: chr, start, end, gc, len). Verify column indexing against bam_to_frag_bed.sh output: `$1=chr, $2=start, $3=end, $4=gc_fraction, $5=fragment_length`. So length filtering should use `$5` not `$4`. Correct accordingly.

- [ ] **Step 6: Add window count rule**

`frag_window_count`:
- Input: short BED, long BED, keep_5mb.bed
- Output: short count tmp, long count tmp
- Script: `scripts/frag_window_int.sh` (from cfdna)

- [ ] **Step 7: Add count merge rule**

`frag_count_merge`:
- Input: expand count tmps across all libraries × {short, long}
- Output: frag_counts.tsv
- Script: `scripts/count_merge.sh` (from cfdna)

- [ ] **Step 8: Add ratio normalization rule**

`frag_ratio_normalize`:
- Input: frag_counts.tsv
- Output: ratios.tsv
- Script: `scripts/make_ratios.R` (from cfdna)

- [ ] **Step 9: Tangle and verify**

```bash
emacsclient --socket-name ~/.emacs.d/server/server -e '(progn (find-file "/home/jeszyman/repos/frag/frag.org") (org-babel-tangle))'
```

Verify `workflows/frag.smk` contains all fragment length rules, and all scripts exist in `scripts/`.

- [ ] **Step 10: Commit**

```bash
git add frag.org workflows/frag.smk scripts/make_gc_map_bind.sh scripts/gc_distro.R scripts/make_healthy_gc_summary.R scripts/sample_frags_by_gc.R scripts/frag_window_sum.sh scripts/frag_window_int.sh scripts/count_merge.sh scripts/make_ratios.R
git commit -m "feat: add fragment length pipeline rules from cfdna (GC normalization through ratios)"
```

---

## Chunk 4: End Motif Rules

### Task 6: Add end motif rules from cfdna

**Files:**
- Modify: `/home/jeszyman/repos/frag/frag.org`
- Reference: `/home/jeszyman/repos/cfdna/scripts/sample_motifs.sh`, `/home/jeszyman/repos/cfdna/scripts/end_motif_mat.R`
- Reference: `/home/jeszyman/repos/cfdna/cfdna.org` (end motifs section for snakemake rules)

**Context:** End motif analysis from cfdna: (1) sample 5' motifs from filtered BAM per library, (2) build motif frequency matrix across libraries. Refactor to frag_ prefix and style guide.

- [ ] **Step 1: Add sample motifs rule and script**

Under `*** End motifs`, add `**** Sample motifs`:

Snakemake rule `frag_sample_motifs`:
- message: "Sample 5-prime end motifs from filtered BAM"
- conda: CONDA_FRAG
- input: filtered BAM + reference FASTA
- log, benchmark
- params: script path, n_motif (4), n_reads (from config or default), seed, threads
- threads: 4
- output: motif txt file
- shell: calls script with positional args

Script tangle block for `scripts/sample_motifs.sh` — refactored from cfdna version with proper arg parsing, set -euo pipefail, quoting:

```bash
#!/usr/bin/env bash
set -euo pipefail

in_bam="$1"
in_fasta="$2"
n_motif="$3"
n_reads="$4"
seed="$5"
threads="$6"
out_merged="$7"

forward_motif() {
    local in_bam="$1" seed="$2" threads="$3" n_reads="$4" in_fasta="$5" n_motif="$6"
    local f_reads=$(( 3 * n_reads ))
    local factor
    factor=$(samtools idxstats "$in_bam" \
        | cut -f3 \
        | awk -v nreads="$f_reads" 'BEGIN {total=0} {total += $1} END {print nreads/total}')

    samtools view \
        --with-header \
        --min-MQ 60 \
        --require-flags 65 \
        --subsample "$factor" \
        --subsample-seed "$seed" \
        --threads "$threads" "$in_bam" \
      | bedtools bamtobed -i stdin \
      | head -n "$n_reads" \
      | bedtools getfasta -bed stdin -fi "$in_fasta" \
      | sed "1d; n; d" \
      | sed -E "s/(.{$n_motif}).*/\1/"
}

reverse_motif() {
    local in_bam="$1" seed="$2" threads="$3" n_reads="$4" in_fasta="$5" n_motif="$6"
    local f_reads=$(( 3 * n_reads ))
    local factor
    factor=$(samtools idxstats "$in_bam" \
        | cut -f3 \
        | awk -v nreads="$f_reads" 'BEGIN {total=0} {total += $1} END {print nreads/total}')

    samtools view \
        --with-header \
        --min-MQ 60 \
        --require-flags 129 \
        --subsample "$factor" \
        --subsample-seed "$seed" \
        --threads "$threads" "$in_bam" \
      | bedtools bamtobed -i stdin \
      | head -n "$n_reads" \
      | bedtools getfasta -bed stdin -fi "$in_fasta" \
      | sed "1d; n; d" \
      | sed -E "s/.*(.{$n_motif})/\1/" \
      | tr ACGT TGCA \
      | rev
}

forward_motif "$in_bam" "$seed" "$threads" "$n_reads" "$in_fasta" "$n_motif" > "$out_merged"
reverse_motif "$in_bam" "$seed" "$threads" "$n_reads" "$in_fasta" "$n_motif" >> "$out_merged"
```

- [ ] **Step 2: Add motif matrix rule and script**

Under `**** Motif matrix`:

Snakemake rule `frag_motif_matrix`:
- Input: lambda expanding motif files across libraries
- Output: motif frequency matrix TSV
- Script: `scripts/end_motif_mat.R`

Script from cfdna's end_motif_mat.R, cleaned up (remove hardcoded test paths):

```r
#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
motif_str <- args[1]
motif_tsv <- args[2]

library(tidyverse)

possible_motifs <- expand.grid(rep(list(c("A", "G", "T", "C")), 4)) %>%
  as_tibble() %>%
  mutate(motif = paste0(Var1, Var2, Var3, Var4)) %>%
  select(motif) %>%
  arrange(motif)

motif_files <- strsplit(motif_str, " ")[[1]]
names(motif_files) <- substr(gsub("^.*lib", "lib", motif_files), 1, 6)

ingest_motif <- function(motif_file) {
  read_tsv(motif_file, col_names = c("motif")) %>%
    group_by(motif) %>%
    summarise(count = n()) %>%
    mutate(fract = count / sum(count)) %>%
    select(motif, fract)
}

motif_tibs <- lapply(motif_files, ingest_motif)

motifs <- bind_rows(motif_tibs, .id = "library") %>%
  pivot_wider(names_from = library, values_from = fract) %>%
  filter(motif %in% possible_motifs$motif)

write_tsv(motifs, motif_tsv)
```

- [ ] **Step 3: Tangle and verify**

```bash
emacsclient --socket-name ~/.emacs.d/server/server -e '(progn (find-file "/home/jeszyman/repos/frag/frag.org") (org-babel-tangle))'
```

Verify scripts/sample_motifs.sh and scripts/end_motif_mat.R exist.

- [ ] **Step 4: Commit**

```bash
git add frag.org workflows/frag.smk scripts/sample_motifs.sh scripts/end_motif_mat.R
git commit -m "feat: add end motif rules from cfdna (sample motifs + motif matrix)"
```

---

## Chunk 5: Test Wrapper and Config

### Task 7: Rebuild test.smk wrapper

**Files:**
- Modify: `/home/jeszyman/repos/frag/frag.org` (Testing section)
- Reference: `/home/jeszyman/repos/emseq/workflows/test.smk` (structural pattern)

**Context:** The test wrapper defines rule all, loads config, sets up constants, includes frag.smk. Frag already has a working test.smk — extend it to include fragment length and end motif outputs in rule all.

- [ ] **Step 1: Update test.smk tangle block**

Under `** Testing > *** Full test`, update the snakemake tangle block (`#+begin_src snakemake :tangle ./workflows/test.smk`).

Key additions to the existing test.smk:
1. Add `CONDA_FRAG_END_MOTIF` if separate env needed (otherwise just CONDA_FRAG)
2. Add fragment length outputs to `rule all`:
   - GC-filtered bins: `f"{D_FRAG}/ref/keep_5mb.bed"`
   - GC distros: expand `f"{D_FRAG}/frags/{{library_id}}.{{ref_name}}.gc_distro.csv"`
   - Healthy GC: `f"{D_FRAG}/frags/{{group}}.healthy_med.rds"`
   - Sampled beds: expand `f"{D_FRAG}/beds/{{library_id}}.{{ref_name}}.sampled_frag.bed"`
   - Short/long beds: expand
   - Count tmps: expand
   - Merged counts: `f"{D_FRAG}/frags/{{group}}.frag_counts.tsv"`
   - Ratios: `f"{D_FRAG}/frags/{{group}}.ratios.tsv"`
3. Add end motif outputs to `rule all`:
   - Per-library motif files: expand
   - Motif matrix: expand per group

Note: The exact path patterns must match what the rules in frag.smk produce. Align carefully with the rule output paths defined in Task 5 and Task 6.

The wrapper also needs config keys for:
- `gc5mb`: path to GC 5MB BED file
- `blklist`: path to blacklist
- End motif params (n_motif, n_reads, seed) — can be defaults in the rule

- [ ] **Step 2: Update config/test.yaml**

Add required keys:

```yaml
conda:
  frag: "${HOME}/repos/frag/config/frag-conda-env.yaml"

directories:
  inputs: tests/full/inputs
  frag: tests/full
  logs: tests/logs
  benchmark: tests/benchmark

repos:
  frag: "${HOME}/repos/frag"

sample-tsv-path: ${HOME}/repos/frag/data/test-samples.tsv
threads: 4

frag_ref_assemblies:
  chr22:
    url: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
    name: chr22
    input: chr22-test.fa.gz

gc5mb: tests/full/inputs/gc5mb.bed
blklist: tests/full/inputs/hg38-blacklist.v2.bed.gz
```

Note: gc5mb test data needs to be in `tests/full/inputs/`. Either get_test_data.sh should generate it, or a minimal version for chr22 needs to be created.

- [ ] **Step 3: Tangle and run dry-run**

```bash
emacsclient --socket-name ~/.emacs.d/server/server -e '(progn (find-file "/home/jeszyman/repos/frag/frag.org") (org-babel-tangle))'
conda run -n basecamp snakemake -s workflows/test.smk --dry-run
```

Fix any errors iteratively until dry-run passes. Common issues:
- Missing config keys
- Path mismatches between rule outputs and rule all expand patterns
- Wildcard resolution failures
- Missing input files in expand patterns

- [ ] **Step 4: Commit**

```bash
git add frag.org workflows/test.smk config/test.yaml
git commit -m "refactor: update test wrapper with full pipeline outputs and config"
```

---

## Chunk 6: Cleanup

### Task 8: Remove legacy workflow/ directory and orphaned files

**Files:**
- Delete: `workflow/` directory (frag_bed.smk, frag_counts.smk, int_test.smk, int_test.yaml, scripts/)
- Delete: any orphaned scripts in `scripts/` not produced by tangle
- Modify: `/home/jeszyman/repos/frag/frag.org` (remove any org blocks that tangled to workflow/scripts/)

- [ ] **Step 1: Identify files to remove**

List contents of `workflow/` directory. These are all legacy files that have been superseded by the restructured `workflows/frag.smk`.

- [ ] **Step 2: Remove legacy workflow/ directory**

```bash
rm -rf workflow/
```

- [ ] **Step 3: Remove any org blocks targeting workflow/scripts/**

Search frag.org for any `:tangle ./workflow/` references and remove those blocks (the content should already have been migrated to the new structure).

- [ ] **Step 4: Remove orphaned scripts**

Compare `scripts/` contents against tangle targets in frag.org. Remove any files not produced by current tangle blocks (e.g., `scripts/<RULE WORKFLOW NAME>.R` placeholder).

- [ ] **Step 5: Tangle to verify clean state**

```bash
emacsclient --socket-name ~/.emacs.d/server/server -e '(progn (find-file "/home/jeszyman/repos/frag/frag.org") (org-babel-tangle))'
```

- [ ] **Step 6: Run dry-run to verify nothing broke**

```bash
conda run -n basecamp snakemake -s workflows/test.smk --dry-run
```

- [ ] **Step 7: Commit**

```bash
git add -A
git commit -m "cleanup: remove legacy workflow/ directory and orphaned scripts"
```

### Task 9: Final verification

- [ ] **Step 1: Full dry-run**

```bash
conda run -n basecamp snakemake -s workflows/test.smk --dry-run
```

Expected: passes with all rules listed in the DAG.

- [ ] **Step 2: Verify org structure**

Read frag.org and verify:
- Top-level structure matches spec
- All snakemake blocks under `** Fragmentomics sequence processing` have correct header-args inheritance
- All script blocks have correct tangle targets
- `** Fragment Length Window Optimization` section is untouched
- No references to `workflow/` directory remain

- [ ] **Step 3: Verify file tree**

```bash
ls -la scripts/
ls -la workflows/
ls -la config/
```

Expected:
- `scripts/` contains only tangle outputs
- `workflows/` contains only frag.smk and test.smk
- `config/` contains frag-conda-env.yaml, test.yaml (and optionally end motif env)
- No `workflow/` directory

- [ ] **Step 4: Final commit if needed**

```bash
git status
# If any remaining changes:
git add -A
git commit -m "cleanup: final verification and file tree cleanup"
```
