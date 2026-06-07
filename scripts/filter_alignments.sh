#!/usr/bin/env bash
set -euo pipefail

# filter_alignments.sh
# Filter BAM to MAPQ>=30 reads within keep-bed regions, fixmate, and re-sort.
# Usage: filter_alignments.sh <in.bam> <keep.bed> <threads> <out.bam>

in_bam="$1"
keep_bed="$2"
threads="$3"
out_bam="$4"

# Precondition: BAM @SQ names and keep_bed chromosomes must share at least one
# name. samtools view -L silently drops every read on a pure identifier mismatch
# (typically 'chr' prefix), producing a header-only BAM with no error.
bam_sn=$(samtools view -H "$in_bam" | awk '$1=="@SQ" {sub("SN:","",$2); print $2}' | sort -u)
bed_chroms=$(cut -f1 "$keep_bed" | sort -u)
if [[ -z "$(comm -12 <(echo "$bam_sn") <(echo "$bed_chroms"))" ]]; then
  echo "ERROR: filter_alignments: $in_bam @SQ and $keep_bed chromosomes share zero overlap" >&2
  echo "  BAM sample: $(echo "$bam_sn"     | head -3 | tr '\n' ' ')" >&2
  echo "  BED sample: $(echo "$bed_chroms" | head -3 | tr '\n' ' ')" >&2
  echo "  Likely cause: 'chr' prefix mismatch. Re-stage $keep_bed or reheader the BAM to match." >&2
  exit 1
fi

# Samtools sort temp files go to CWD by default, which may be a small/home
# filesystem. On large BAMs the name-sort spill can exhaust disk and surface
# as bogus "Illegal seek" write errors on the tmp file. Force tmp to a known
# large scratch dir when TMPDIR_SAMTOOLS is set, else fall back to $TMPDIR, else /tmp.
sort_tmp="${TMPDIR_SAMTOOLS:-${TMPDIR:-/tmp}}/samtools_sort.$$"
mkdir -p "$(dirname "$sort_tmp")"

samtools view -@ "$threads" -b -h -L "$keep_bed" -q 30 "$in_bam" \
  | samtools sort -@ "$threads" -n -T "${sort_tmp}.name" -o - - \
  | samtools fixmate -@ "$threads" - - \
  | samtools sort -@ "$threads" -T "${sort_tmp}.coord" -o "$out_bam" -

samtools index -@ "$threads" "$out_bam"

# Postcondition: output BAM must contain reads (not just a header).
# Uses idxstats for O(1)-per-chrom speed rather than a full scan.
reads=$(samtools idxstats "$out_bam" | awk '{s += $3} END {print s+0}')
if [[ "$reads" -eq 0 ]]; then
  echo "ERROR: filter_alignments: $out_bam contains zero reads after filtering" >&2
  echo "  Input BAM: $in_bam" >&2
  echo "  Keep BED:  $keep_bed" >&2
  echo "  Possible causes: MAPQ threshold too strict, keep regions do not overlap reads, or upstream filter was already empty." >&2
  exit 1
fi
echo "[filter_alignments] OK: $reads reads in $out_bam" >&2
