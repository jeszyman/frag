#!/usr/bin/env bash
# 5' fragment-end 4-mer per read, strand-aware, reference-based. One pass over all paired primary reads.
#   extract_5p_motifs.sh <bam> <fasta> <out_counts.tsv> [threads] [region]
# samtools + bedtools from PATH, or TOOLS_BIN=/path/to/bin. Output: motif<TAB>count, plus an OTHER row
# (non-ACGT 4-mers, dropped). MAPQ>=30, paired, primary, mapped, non-duplicate; both mates counted.
set -euo pipefail
bam=$1; fasta=$2; out=$3; threads=${4:-3}; region=${5:-}
BIN=${TOOLS_BIN:-$(dirname "$(command -v samtools)")}
"$BIN/samtools" view -b -q 30 -f 1 -F 3340 --threads "$threads" "$bam" $region \
 | "$BIN/bedtools" bamtobed -i stdin \
 | awk 'BEGIN{OFS="\t"} $6=="+"{print $1,$2,$2+4,".",".",$6} $6=="-"{if($3>=4)print $1,$3-4,$3,".",".",$6}' \
 | "$BIN/bedtools" getfasta -fi "$fasta" -bed stdin -s -tab \
 | awk -F'\t' '{m=toupper($2); if(m ~ /^[ACGT][ACGT][ACGT][ACGT]$/) c[m]++; else o++} END{for(k in c) print k"\t"c[k]; print "OTHER\t"o+0}' \
 | sort > "$out"
