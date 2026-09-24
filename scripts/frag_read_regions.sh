#!/usr/bin/env bash
set -euo pipefail

# frag_read_regions.sh
# Regions the filtered BAM keeps: the autosomes (index names matching
# ^chr[0-9]+$) minus every blacklisted base, sorted and merged. Windows that
# contain blacklisted bases are kept; only the blacklisted bases are removed.
# Usage: frag_read_regions.sh <ref.fa.fai> <blacklist.bed[.gz]> <out.bed>

fai="$1"
blklist="$2"
out_bed="$3"

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

awk 'BEGIN {OFS="\t"} $1 ~ /^chr[0-9]+$/ {print $1, 0, $2}' "$fai" \
  | sort -k1,1 -k2,2n > "$tmpdir/autosomes.bed"
if [[ ! -s "$tmpdir/autosomes.bed" ]]; then
  echo "ERROR: frag_read_regions: $fai has no autosome (names matching ^chr[0-9]+\$)" >&2
  exit 1
fi

zcat -f "$blklist" | cut -f1-3 | sort -k1,1 -k2,2n > "$tmpdir/blacklist.bed"
if [[ ! -s "$tmpdir/blacklist.bed" ]]; then
  echo "ERROR: frag_read_regions: blacklist $blklist is empty" >&2
  exit 1
fi

# Precondition: the blacklist and the autosomes share chromosome names.
# bedtools subtract removes nothing on a pure naming mismatch ('chr' prefix).
shared=$(comm -12 <(cut -f1 "$tmpdir/autosomes.bed" | sort -u) <(cut -f1 "$tmpdir/blacklist.bed" | sort -u))
if [[ -z "$shared" ]]; then
  echo "ERROR: frag_read_regions: $blklist and the autosomes of $fai share no chromosome name" >&2
  echo "  Likely cause: 'chr' prefix mismatch." >&2
  exit 1
fi

bedtools subtract -a "$tmpdir/autosomes.bed" -b "$tmpdir/blacklist.bed" \
  | sort -k1,1 -k2,2n \
  | bedtools merge -i stdin > "$out_bed"

total_bp=$(awk '{s += $3 - $2} END {print s+0}' "$tmpdir/autosomes.bed")
kept_bp=$(awk '{s += $3 - $2} END {print s+0}' "$out_bed")
if [[ "$kept_bp" -eq 0 ]]; then
  echo "ERROR: frag_read_regions: no autosome bases left after the blacklist" >&2
  exit 1
fi
echo "[frag_read_regions] OK: $kept_bp of $total_bp autosome bp kept in $(wc -l < "$out_bed") intervals" >&2
