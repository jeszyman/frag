#!/usr/bin/env bash
set -euo pipefail

# frag_end_motifs.sh
# Count the 5' fragment-end 4-mer of every qualifying read, strand-aware and
# read from the reference. Adapted from extract_5p_motifs.sh by
# Axel Hidalgo, commit 7a5767ab1 (2026-09-18).
#
# Qualifying read: paired, primary, mapped, mate mapped, non-duplicate,
# MAPQ >= 30 (samtools -q 30 -f 1 -F 3340). Both mates count, so each
# fragment gives two ends. Plus-strand read: reference [start, start+4).
# Minus-strand read: reference [end-4, end), reverse-complemented by
# bedtools getfasta -s. 4-mers with a non-ACGT base are counted as OTHER.
#
# max_ends > 0 and fewer than the qualifying reads: samtools --subsample keeps
# a fraction max_ends/qualifying reads, chosen by a hash of the read name, so
# mates stay together and every chromosome is sampled. max_ends 0 or none
# counts all.
#
# Output: 256 rows "motif<TAB>count" (A<C<G<T order, zeros included), then
# "OTHER<TAB>count".
# Usage: frag_end_motifs.sh <bam> <fasta> <out.tsv> <threads> <max_ends|none> <seed>

in_bam="$1"
in_fasta="$2"
out_tsv="$3"
threads="$4"
max_ends="$5"
seed="$6"

# max_ends "none" (or null/None/empty) means every end, the same as 0.
case "$max_ends" in none|None|null|"") max_ends=0 ;; esac
if ! [[ "$max_ends" =~ ^[0-9]+$ ]]; then
  echo "ERROR: frag_end_motifs: max_ends must be a non-negative integer or none, got '$max_ends'" >&2
  exit 1
fi

read_filter=(-q 30 -f 1 -F 3340)

# Precondition: BAM holds reads.
in_reads=$(samtools idxstats "$in_bam" 2>/dev/null | awk '{s += $3} END {print s+0}')
if [[ "$in_reads" -eq 0 ]]; then
  echo "ERROR: frag_end_motifs: $in_bam contains zero reads" >&2
  exit 1
fi

# Precondition: BAM @SQ names overlap the FASTA index; getfasta otherwise
# drops every interval.
if [[ -f "${in_fasta}.fai" ]]; then
  bam_sn=$(samtools view -H "$in_bam" | awk '$1=="@SQ" {sub("SN:","",$2); print $2}' | sort -u)
  fai_chroms=$(cut -f1 "${in_fasta}.fai" | sort -u)
  if [[ -z "$(comm -12 <(echo "$bam_sn") <(echo "$fai_chroms"))" ]]; then
    echo "ERROR: frag_end_motifs: $in_bam @SQ and ${in_fasta}.fai share zero names" >&2
    exit 1
  fi
fi

subsample=()
if (( max_ends > 0 )); then
  n_qual=$(samtools view -c "${read_filter[@]}" --threads "$threads" "$in_bam")
  if (( n_qual > max_ends )); then
    frac=$(awk -v m="$max_ends" -v n="$n_qual" 'BEGIN {printf "%.8f", m / n}')
    subsample=(--subsample "$frac" --subsample-seed "$seed")
    echo "[frag_end_motifs] subsampling $n_qual qualifying reads at fraction $frac" >&2
  fi
fi

samtools view -b "${read_filter[@]}" ${subsample[@]+"${subsample[@]}"} --threads "$threads" "$in_bam" \
  | bedtools bamtobed -i stdin \
  | awk 'BEGIN {OFS="\t"}
         $6=="+" {print $1, $2, $2+4, ".", ".", $6}
         $6=="-" {if ($3 >= 4) print $1, $3-4, $3, ".", ".", $6}' \
  | bedtools getfasta -fi "$in_fasta" -bed stdin -s -tab \
  | awk -F'\t' '
      {m = toupper($2); if (m ~ /^[ACGT][ACGT][ACGT][ACGT]$/) c[m]++; else other++}
      END {
        split("A C G T", b, " ")
        for (i = 1; i <= 4; i++) for (j = 1; j <= 4; j++)
          for (k = 1; k <= 4; k++) for (l = 1; l <= 4; l++) {
            m = b[i] b[j] b[k] b[l]; print m "\t" c[m] + 0
          }
        print "OTHER\t" other + 0
      }' > "$out_tsv"

# Postcondition: 257 rows and at least one counted end.
rows=$(wc -l < "$out_tsv")
total=$(awk -F'\t' '{s += $2} END {print s+0}' "$out_tsv")
if [[ "$rows" -ne 257 || "$total" -eq 0 ]]; then
  echo "ERROR: frag_end_motifs: $out_tsv has $rows rows and $total ends" >&2
  exit 1
fi
echo "[frag_end_motifs] OK: $total ends in $out_tsv" >&2
