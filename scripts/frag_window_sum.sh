#!/usr/bin/env bash
set -euo pipefail

# frag_window_sum.sh
# Partition sampled fragment BED into short (100-150bp) and long (151-220bp).
# Sampled BED columns: chr, start, end, len, gc_strata
# Length is in column 4 ($4).
# Usage: frag_window_sum.sh <input.bed> <short.bed> <long.bed>

input_frag="$1"
output_short="$2"
output_long="$3"

awk '{ if ($4 >= 100 && $4 <= 150) print $0 }' "$input_frag" > "$output_short"
awk '{ if ($4 >= 151 && $4 <= 220) print $0 }' "$input_frag" > "$output_long"
