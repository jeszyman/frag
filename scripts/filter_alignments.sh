#!/usr/bin/env bash
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/frag/frag.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-13 16:35:03
# ============================================================

set -euo pipefail

# filter_alignments.sh
# Filter BAM to MAPQ>=30 reads within keep-bed regions, fixmate, and re-sort.
# Usage: filter_alignments.sh <in.bam> <keep.bed> <threads> <out.bam>

in_bam="$1"
keep_bed="$2"
threads="$3"
out_bam="$4"

samtools view -@ "$threads" -b -h -L "$keep_bed" -q 30 "$in_bam" \
  | samtools sort -@ "$threads" -n -o - - \
  | samtools fixmate -@ "$threads" - - \
  | samtools sort -@ "$threads" -o "$out_bam" -

samtools index -@ "$threads" "$out_bam"
