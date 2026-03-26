#!/usr/bin/env bash
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

    # head causes SIGPIPE (141) upstream under pipefail; write to temp then process
    local tmpbed
    tmpbed=$(mktemp)
    samtools view \
        --with-header \
        --min-MQ 60 \
        --require-flags 65 \
        --subsample "$factor" \
        --subsample-seed "$seed" \
        --threads "$threads" "$in_bam" \
      | bedtools bamtobed -i stdin \
      | head -n "$n_reads" > "$tmpbed" || true

    bedtools getfasta -bed "$tmpbed" -fi "$in_fasta" \
      | sed "1d; n; d" \
      | sed -E "s/(.{$n_motif}).*/\1/"
    rm -f "$tmpbed"
}

reverse_motif() {
    local in_bam="$1" seed="$2" threads="$3" n_reads="$4" in_fasta="$5" n_motif="$6"
    local f_reads=$(( 3 * n_reads ))
    local factor
    factor=$(samtools idxstats "$in_bam" \
        | cut -f3 \
        | awk -v nreads="$f_reads" 'BEGIN {total=0} {total += $1} END {print nreads/total}')

    local tmpbed
    tmpbed=$(mktemp)
    samtools view \
        --with-header \
        --min-MQ 60 \
        --require-flags 129 \
        --subsample "$factor" \
        --subsample-seed "$seed" \
        --threads "$threads" "$in_bam" \
      | bedtools bamtobed -i stdin \
      | head -n "$n_reads" > "$tmpbed" || true

    bedtools getfasta -bed "$tmpbed" -fi "$in_fasta" \
      | sed "1d; n; d" \
      | sed -E "s/.*(.{$n_motif})/\1/" \
      | tr ACGT TGCA \
      | rev
    rm -f "$tmpbed"
}

forward_motif "$in_bam" "$seed" "$threads" "$n_reads" "$in_fasta" "$n_motif" > "$out_merged"
reverse_motif "$in_bam" "$seed" "$threads" "$n_reads" "$in_fasta" "$n_motif" >> "$out_merged"
