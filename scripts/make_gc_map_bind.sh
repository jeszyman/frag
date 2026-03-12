#!/usr/bin/env bash
set -euo pipefail

# make_gc_map_bind.sh
# Filter 5Mb GC bins by blacklist and minimum GC content.
# Usage: make_gc_map_bind.sh <gc5mb.bed> <blacklist.bed> <output_keep.bed>

gc5mb="$1"
blklist="$2"
keep="$3"

bedtools intersect -a "$gc5mb" -b "$blklist" -v -wa \
  | grep -v _ \
  | awk '{ if ($4 >= 0.3) print $0 }' > "$keep"
