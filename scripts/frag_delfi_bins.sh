#!/usr/bin/env bash
set -euo pipefail

# frag_delfi_bins.sh
# DELFI bins on one reference: the rows of the bin table (header
# chr start end arm gc map blacklisted_bases, 0-based half-open) whose
# chromosome is in the reference index and whose end lies within it, written
# as a BED of chr, start, end, gc.
# Usage: frag_delfi_bins.sh <delfi_bins.tsv> <ref.fa.fai> <out.bed>

bins="$1"
fai="$2"
out_bed="$3"

header=$(head -n 1 "$bins")
if [[ "$header" != $'chr\tstart\tend\tarm\tgc\tmap\tblacklisted_bases' ]]; then
  echo "ERROR: frag_delfi_bins: $bins header is '$header', expected chr start end arm gc map blacklisted_bases" >&2
  exit 1
fi

awk -F'\t' 'BEGIN {OFS="\t"}
  NR == FNR {len[$1] = $2; next}
  FNR > 1 && ($1 in len) && $3 <= len[$1] {print $1, $2, $3, $5}' "$fai" "$bins" \
  | sort -k1,1 -k2,2n > "$out_bed"

# Postcondition: at least one bin lies on the reference.
n=$(wc -l < "$out_bed")
if [[ "$n" -eq 0 ]]; then
  echo "ERROR: frag_delfi_bins: no bin of $bins lies on a contig of $fai" >&2
  echo "  Likely cause: 'chr' prefix mismatch or a bin table for another genome." >&2
  exit 1
fi
echo "[frag_delfi_bins] OK: $n of $(( $(wc -l < "$bins") - 1 )) bins on the reference in $out_bed" >&2
