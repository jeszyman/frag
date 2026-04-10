#!/usr/bin/env bash
set -euo pipefail

# check_ids.sh
# Verify that all pipeline inputs share a chromosome-naming convention with the
# reference FASTA. Catches the classic "chr1" vs "1" silent-failure bug where
# bedtools nuc, samtools view -L, and bedtools intersect emit nothing on a pure
# identifier mismatch — producing empty outputs and opaque downstream errors.
#
# The check is non-empty intersection: for each file, its chrom set must share
# at least one name with the reference's .fai. Subsets are fine (regions BED
# with only autosomes against a reference with autosomes + X + Y + decoys).
# Zero overlap is a hard halt.
#
# Usage:
#   check_ids.sh <ref.fa> <regions.bed> <blacklist.bed.gz> <cytoband.txt> \
#                <sentinel.ok> <bam1> [bam2 ...]

ref="$1"; regions="$2"; blklist="$3"; cyto="$4"; sentinel="$5"
shift 5

[[ -f "${ref}.fai" ]] || samtools faidx "$ref"
ref_chroms=$(cut -f1 "${ref}.fai" | sort -u)

check() {
  local label="$1" file_chroms="$2"
  local overlap
  overlap=$(comm -12 <(echo "$ref_chroms") <(echo "$file_chroms" | sort -u) | wc -l)
  if [[ "$overlap" -eq 0 ]]; then
    echo "ERROR: check_ids: $label has zero chrom overlap with reference ${ref}.fai" >&2
    echo "  ref sample:    $(echo "$ref_chroms"  | head -3 | tr '\n' ' ')" >&2
    echo "  $label sample: $(echo "$file_chroms" | head -3 | tr '\n' ' ')" >&2
    echo "  Likely cause: 'chr' prefix mismatch. All inputs must use the reference's convention." >&2
    exit 1
  fi
  echo "[check_ids] OK: $label ($overlap chroms match reference)" >&2
}

check "regions"   "$(cut -f1 "$regions")"
check "blacklist" "$(zcat -f "$blklist" | cut -f1)"
check "cytoband"  "$(cut -f1 "$cyto")"

for bam in "$@"; do
  check "$(basename "$bam")" \
    "$(samtools view -H "$bam" | awk '$1=="@SQ" {sub("SN:","",$2); print $2}')"
done

mkdir -p "$(dirname "$sentinel")"
touch "$sentinel"
echo "[check_ids] ALL OK — sentinel written to $sentinel" >&2
