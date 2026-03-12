#!/usr/bin/env bash
set -euo pipefail

# bam_to_frag_bed.sh
# Convert a name-sorted/fixmated BAM to a fragment-level BED with GC and length.
# Output columns: chr, start, end, gc_fraction, fragment_length
# Usage: bam_to_frag_bed.sh <in.bam> <ref.fa> <out.bed>

input_bam="$1"
params_fasta="$2"
output_frag_bed="$3"

bedtools bamtobed -bedpe -i "$input_bam" \
  | awk '$1==$4 {print $0}' \
  | awk '$2 < $6 {print $0}' \
  | awk -v OFS='\t' '{print $1,$2,$6}' \
  | bedtools nuc -fi "$params_fasta" -bed stdin \
  | awk -v OFS='\t' '{print $1,$2,$3,$5,$12}' \
  | sed '1d' > "$output_frag_bed"
