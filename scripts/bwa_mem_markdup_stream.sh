#!/usr/bin/env bash
set -euo pipefail

# bwa_mem_markdup_stream.sh
# Align paired FASTQs with BWA MEM, coordinate sort, and index.
# Usage: bwa_mem_markdup_stream.sh <ref.fa> <r1.fq.gz> <r2.fq.gz> <out.bam> <threads>

ref="$1"
r1="$2"
r2="$3"
out_bam="$4"
threads="$5"

bwa mem -M -t "$threads" "$ref" "$r1" "$r2" \
  | samtools view -@ 4 -Sb - \
  | samtools sort -@ 4 -o "$out_bam" -

samtools index -@ 4 "$out_bam"
