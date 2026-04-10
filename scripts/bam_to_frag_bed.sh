#!/usr/bin/env bash
set -euo pipefail

# bam_to_frag_bed.sh
# Convert a name-sorted/fixmated BAM to a fragment-level BED with GC and length.
# Output columns: chr, start, end, gc_fraction, fragment_length
# Usage: bam_to_frag_bed.sh <in.bam> <ref.fa> <out.bed>

input_bam="$1"
params_fasta="$2"
output_frag_bed="$3"

# Precondition: input BAM must contain reads. A header-only BAM propagates
# silently through bedtools bamtobed and produces an empty frag.bed that
# breaks downstream R scripts with opaque "NULL data.table" errors.
in_reads=$(samtools idxstats "$input_bam" 2>/dev/null | awk '{s += $3} END {print s+0}')
if [[ "$in_reads" -eq 0 ]]; then
  echo "ERROR: bam_to_frag_bed: $input_bam contains zero reads" >&2
  exit 1
fi

# Precondition: BAM @SQ names must overlap FASTA index chromosomes. bedtools
# nuc silently emits bogus rows when a read's chrom is not in the FASTA.
if [[ -f "${params_fasta}.fai" ]]; then
  bam_sn=$(samtools view -H "$input_bam" | awk '$1=="@SQ" {sub("SN:","",$2); print $2}' | sort -u)
  fai_chroms=$(cut -f1 "${params_fasta}.fai" | sort -u)
  if [[ -z "$(comm -12 <(echo "$bam_sn") <(echo "$fai_chroms"))" ]]; then
    echo "ERROR: bam_to_frag_bed: $input_bam @SQ and ${params_fasta}.fai share zero overlap" >&2
    echo "  BAM sample:   $(echo "$bam_sn"     | head -3 | tr '\n' ' ')" >&2
    echo "  FASTA sample: $(echo "$fai_chroms" | head -3 | tr '\n' ' ')" >&2
    echo "  Likely cause: 'chr' prefix mismatch between BAM header and reference FASTA." >&2
    exit 1
  fi
fi

bedtools bamtobed -bedpe -i "$input_bam" \
  | awk '$1==$4 {print $0}' \
  | awk '$2 < $6 {print $0}' \
  | awk -v OFS='\t' '{print $1,$2,$6}' \
  | bedtools nuc -fi "$params_fasta" -bed stdin \
  | awk -v OFS='\t' '{print $1,$2,$3,$5,$12}' \
  | sed '1d' > "$output_frag_bed"

# Postcondition: output BED must be non-empty.
lines=$(wc -l < "$output_frag_bed")
if [[ "$lines" -eq 0 ]]; then
  echo "ERROR: bam_to_frag_bed: $output_frag_bed is empty" >&2
  echo "  Input BAM read count: $in_reads" >&2
  echo "  Likely cause: no intra-chromosome paired reads, or bedtools nuc dropped rows due to identifier mismatch." >&2
  exit 1
fi
echo "[bam_to_frag_bed] OK: $lines fragments in $output_frag_bed" >&2
