#!/usr/bin/env bash
set -euo pipefail

# count_merge.sh
# Merge fragment count files into a single TSV.
# Expects filenames: {library_id}.{ref_name}.cnt_{length}.tmp
# Usage: count_merge.sh <counts_dir> <output.tsv>

counts_dir="${1}"
out_tsv="${2}"

echo -e "library\tlen_class\tchr\tstart\tend\tgc\tcount" > "$out_tsv"

for file in "${counts_dir}"/*.tmp; do
    base=$(basename "$file" .tmp)
    # Parse library_id and len_class from filename: lib001.chr22.cnt_short
    library_id="${base%%.*}"
    len_class="${base##*.cnt_}"
    awk -v lib="$library_id" -v lc="$len_class" \
        'BEGIN{OFS="\t"} {print lib, lc, $0}' "$file" >> "$out_tsv"
done
