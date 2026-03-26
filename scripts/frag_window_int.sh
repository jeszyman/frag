#!/usr/bin/env bash
set -euo pipefail

# frag_window_int.sh
# Count intersections of fragment BED with genomic bins.
# Usage: frag_window_int.sh <fragments.bed> <bins.bed> <output.tmp>

input="$1"
keep_bed="$2"
output="$3"

bedtools intersect -c \
  -a "$keep_bed" \
  -b "$input" > "$output"
