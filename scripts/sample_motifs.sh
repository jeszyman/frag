#!/usr/bin/env bash
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/frag/frag.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-12 15:04:13
# ============================================================

set -euo pipefail

# sample_motifs.sh
# Sample 5-prime end motifs from a filtered BAM file.
# Extracts forward and reverse read motifs, writing to a single output.
# Usage: sample_motifs.sh <bam> <fasta> <n_motif> <n_reads> <seed> <threads> <out.txt>

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
