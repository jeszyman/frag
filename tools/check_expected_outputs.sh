#!/usr/bin/env bash
# check_expected_outputs.sh: compare a fixture run's deterministic outputs to
# the committed checksums in tests/expected_outputs.md5. Paths in that file are
# relative to the data dir; a .gz path is checksummed on its decompressed bytes.
# Usage: check_expected_outputs.sh <data-dir> [--write]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXPECTED="$(dirname "$SCRIPT_DIR")/tests/expected_outputs.md5"
DATA="${1:?usage: check_expected_outputs.sh <data-dir> [--write]}"
MODE="${2:-check}"
FILES=(
  counts/lib001.chr22.cnt_long.tmp
  counts/lib001.chr22.cnt_short.tmp
  counts/lib002.chr22.cnt_long.tmp
  counts/lib002.chr22.cnt_short.tmp
  features/chr22.arm_zscores.csv
  features/chr22.fprofiles.tsv
  features/chr22.motif_diversity.csv
  features/chr22.motif_per_fprofile.tsv
  features/chr22.nmf_2.H.csv
  features/chr22.nmf_2.W.csv
  features/chr22.ratios_normalized.csv
  frags/chr22.frag_counts.tsv
  frags/chr22.ratios.tsv
  frags/lib001.chr22.frag.bed
  frags/lib001.chr22.gc_distro.csv
  frags/lib001.chr22.norm_long.bed
  frags/lib001.chr22.norm_short.bed
  frags/lib001.chr22.sampled_frag.bed
  frags/lib002.chr22.frag.bed
  frags/lib002.chr22.gc_distro.csv
  frags/lib002.chr22.norm_long.bed
  frags/lib002.chr22.norm_short.bed
  frags/lib002.chr22.sampled_frag.bed
  histograms/chr22.count_histogram.csv
  histograms/chr22.freq_histogram.csv
  histograms/lib001.chr22.length_hist.tsv
  histograms/lib002.chr22.length_hist.tsv
  motifs/chr22.all_motifs.tsv
  motifs/chr22.motif_counts.tsv
  motifs/lib001.chr22.motif_counts.tsv
  motifs/lib002.chr22.motif_counts.tsv
  ref/keep_5mb.bed
)

sum_one() {
  case "$1" in
    *.gz) zcat "$DATA/$1" | md5sum | cut -d' ' -f1 ;;
    *)    md5sum "$DATA/$1" | cut -d' ' -f1 ;;
  esac
}

if [ "$MODE" = "--write" ]; then
  for f in "${FILES[@]}"; do echo "$(sum_one "$f")  $f"; done > "$EXPECTED"
  cat "$EXPECTED"
  exit 0
fi

status=0
while read -r want path; do
  if [ ! -f "$DATA/$path" ]; then
    echo "DIFF  $path (missing)"
    status=1
    continue
  fi
  got=$(sum_one "$path")
  if [ "$got" = "$want" ]; then
    echo "ok    $path"
  else
    echo "DIFF  $path (expected $want, got $got)"
    status=1
  fi
done < "$EXPECTED"
exit $status
