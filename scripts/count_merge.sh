#!/usr/bin/env bash
set -euo pipefail

# count_merge.sh
# Merge fragment count files into a single TSV.
# Usage: count_merge.sh <counts_dir> <output.tsv>

counts_dir="${1}"
out_tsv="${2}"

if [ -f "$out_tsv" ]; then rm "$out_tsv"; fi

for file in "${counts_dir}"/*; do
    awk '{{print FILENAME (NF?"\t":"") $0}}' "$file" |
        sed 's/^.*lib/lib/g' |
        sed 's/_.*_/\t/g' |
        sed 's/.tmp//g' |
        sed 's/\.bed//g' >> "$out_tsv"
done

sed -i '1 i\library\tlen_class\tchr\tstart\tend\tgc\tcount' "$out_tsv"
