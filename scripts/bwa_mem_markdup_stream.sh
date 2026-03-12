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
