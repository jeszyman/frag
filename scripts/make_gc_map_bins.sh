#!/usr/bin/env bash
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

# Ensure FASTA is indexed
samtools faidx "$fasta" 2>/dev/null || true

# Precondition: regions BED and FASTA must share at least one chromosome name.
# bedtools nuc silently skips windows whose chrom is not in the FASTA, producing
# an empty keep.bed with only "Feature beyond the length of X size (0 bp)" warnings.
regions_chroms=$(cut -f1 "$regions" | sort -u)
fai_chroms=$(cut -f1 "${fasta}.fai" | sort -u)
if [[ -z "$(comm -12 <(echo "$regions_chroms") <(echo "$fai_chroms"))" ]]; then
  echo "ERROR: make_gc_map_bins: $regions and ${fasta}.fai share zero chromosome names" >&2
  echo "  regions sample: $(echo "$regions_chroms" | head -3 | tr '\n' ' ')" >&2
  echo "  FASTA sample:   $(echo "$fai_chroms"     | head -3 | tr '\n' ' ')" >&2
  echo "  Likely cause: 'chr' prefix mismatch. Strip or add prefix on one side to match." >&2
  exit 1
fi

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

# Postcondition: keep.bed must be non-empty.
lines=$(wc -l < "$keep")
if [[ "$lines" -eq 0 ]]; then
  echo "ERROR: make_gc_map_bins produced empty $keep" >&2
  echo "  Possible causes:" >&2
  echo "    - blacklist removed all bins (verify $blklist chrom naming matches $regions)" >&2
  echo "    - all bins had GC < 0.3" >&2
  echo "    - bedtools nuc silently dropped windows (identifier mismatch against FASTA)" >&2
  exit 1
fi
echo "[make_gc_map_bins] OK: $lines bins kept in $keep" >&2
