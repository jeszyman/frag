#!/usr/bin/env bash
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/frag/frag.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-24 13:15:21
# ============================================================

set -euo pipefail

# make_gc_map_bins.sh
# Create 5Mb genomic windows, compute GC content, filter by blacklist and GC.
# Usage: make_gc_map_bins.sh <regions.bed> <ref.fa> <blacklist.bed.gz> <output_keep.bed>

regions="$1"
fasta="$2"
blklist="$3"
keep="$4"

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

# Create genome file from FASTA index
samtools faidx "$fasta" 2>/dev/null || true
awk 'BEGIN{OFS="\t"} {print $1,$2}' "${fasta}.fai" > "$tmpdir/genome.txt"

# Window into 5Mb bins
bedtools makewindows -b "$regions" -w 5000000 > "$tmpdir/windows.bed"

# Compute GC content per window (column 5 = %GC in bedtools nuc output)
bedtools nuc -fi "$fasta" -bed "$tmpdir/windows.bed" \
  | tail -n +2 \
  | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$5}' > "$tmpdir/gc_windows.bed"

# Remove blacklist-overlapping bins, alt contigs, and low-GC bins
bedtools intersect -a "$tmpdir/gc_windows.bed" -b "$blklist" -v -wa \
  | grep -v '_' \
  | awk '{ if ($4 >= 0.3) print $0 }' > "$keep"
