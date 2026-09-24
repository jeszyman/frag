#!/usr/bin/env bash
set -euo pipefail

# bwa_mem_markdup_stream.sh
# Align paired FASTQs with BWA MEM, add mate tags (fixmate -m), coordinate
# sort, mark duplicates (samtools markdup), and index. Duplicates are flagged
# (0x400), not removed; the filtered BAM drops them.
# Usage: bwa_mem_markdup_stream.sh <ref.fa> <r1.fq.gz> <r2.fq.gz> <out.bam> <threads>

ref="$1"
r1="$2"
r2="$3"
out_bam="$4"
threads="$5"

# Sort temp files go next to the output, not to the working directory.
tmp_prefix="${out_bam%.bam}.tmp.$$"

bwa mem -M -t "$threads" "$ref" "$r1" "$r2" \
  | samtools fixmate -@ 4 -m -u - - \
  | samtools sort -@ 4 -u -T "${tmp_prefix}.sort" - \
  | samtools markdup -@ 4 -T "${tmp_prefix}.markdup" - "$out_bam"

samtools index -@ 4 "$out_bam"

dups=$(samtools view -c -f 1024 "$out_bam")
echo "[bwa_mem_markdup_stream] OK: $out_bam holds $dups reads marked duplicate" >&2
samtools flagstat "$out_bam" >&2
