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

# frag_window_int.sh
# Count intersections of fragment BED with genomic bins.
# Usage: frag_window_int.sh <fragments.bed> <bins.bed> <output.tmp>

input="$1"
keep_bed="$2"
output="$3"

bedtools intersect -c \
  -a "$keep_bed" \
  -b "$input" > "$output"
