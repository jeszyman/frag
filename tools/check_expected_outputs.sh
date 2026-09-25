#!/usr/bin/env bash
# check_expected_outputs.sh: compare a fixture run's deterministic outputs to
# the committed expectations. Byte-stable files are compared by md5 against
# tests/expected_outputs.md5 (a .gz path is summed on its decompressed bytes).
# Floating-point tables whose last digits depend on the BLAS kernel are
# compared numerically against copies in tests/expected/ (same shape, row and
# column labels; numpy.isclose rtol 1e-6, atol 1e-9). Paths are relative to the
# data dir. The check fails if the md5 pin file does not list exactly FILES.
# Usage: check_expected_outputs.sh <data-dir> [--write]
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(dirname "$SCRIPT_DIR")"
EXPECTED="$REPO/tests/expected_outputs.md5"
EXPECTED_DIR="$REPO/tests/expected"
DATA="${1:?usage: check_expected_outputs.sh <data-dir> [--write]}"
MODE="${2:-check}"
FILES=(
  counts/lib001.chr22.cnt_long.tmp
  counts/lib001.chr22.cnt_short.tmp
  counts/lib002.chr22.cnt_long.tmp
  counts/lib002.chr22.cnt_short.tmp
  features/chr22.motif_diversity.csv
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
  ref/chr22.delfi_bins.bed
  ref/chr22.read_regions.bed
)
FLOAT_FILES=(
  features/chr22.arm_zscores.csv
  features/chr22.fprofiles.tsv
  features/chr22.motif_per_fprofile.tsv
  features/chr22.nmf_2.H.csv
  features/chr22.nmf_2.W.csv
)

sum_one() {
  case "$1" in
    *.gz) zcat "$DATA/$1" | md5sum | cut -d' ' -f1 ;;
    *)    md5sum "$DATA/$1" | cut -d' ' -f1 ;;
  esac
}

# Numeric comparison of two tables; prints nothing and exits 0 when they match.
compare_float() {
  python3 - "$1" "$2" <<'EOF'
import sys

import numpy as np
import pandas as pd

want_path, got_path = sys.argv[1:3]
sep = "\t" if want_path.endswith(".tsv") else ","
want = pd.read_csv(want_path, sep=sep, index_col=0)
got = pd.read_csv(got_path, sep=sep, index_col=0)
if want.shape != got.shape:
    sys.exit(f"shape {got.shape}, expected {want.shape}")
if list(want.index) != list(got.index) or list(want.columns) != list(got.columns):
    sys.exit("row or column labels differ")
w, g = want.to_numpy(dtype=float), got.to_numpy(dtype=float)
close = np.isclose(g, w, rtol=1e-6, atol=1e-9, equal_nan=True)
if not close.all():
    sys.exit(f"{(~close).sum()} values differ, largest by {np.nanmax(np.abs(g - w)):.3g}")
EOF
}

if [ "$MODE" = "--write" ]; then
  for f in "${FILES[@]}"; do echo "$(sum_one "$f")  $f"; done > "$EXPECTED"
  for f in "${FLOAT_FILES[@]}"; do
    mkdir -p "$(dirname "$EXPECTED_DIR/$f")"
    cp "$DATA/$f" "$EXPECTED_DIR/$f"
  done
  cat "$EXPECTED"
  exit 0
fi

status=0

# The pin file must list exactly FILES: an empty or truncated file fails.
if ! diff <(printf '%s\n' "${FILES[@]}" | sort) <(awk '{print $2}' "$EXPECTED" | sort) > /dev/null; then
  echo "DIFF  $EXPECTED does not list exactly the ${#FILES[@]} pinned files"
  status=1
fi

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

for f in "${FLOAT_FILES[@]}"; do
  if [ ! -f "$EXPECTED_DIR/$f" ] || [ ! -f "$DATA/$f" ]; then
    echo "DIFF  $f (missing expected copy or output)"
    status=1
  elif msg=$(compare_float "$EXPECTED_DIR/$f" "$DATA/$f" 2>&1); then
    echo "ok    $f (numeric)"
  else
    echo "DIFF  $f ($msg)"
    status=1
  fi
done
exit $status
