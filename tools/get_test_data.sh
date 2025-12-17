#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# Script to collect and pre-process input data for a minimal working example
# of the repo functionality. The input data will persist in the repo. This
# script just documents how they were obtained. The other contents of
# tests/full are NOT committed and are ignored by git.
# -----------------------------------------------------------------------------

parse_args() {
    declare -g script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    declare -g repo_dir="$(cd "${script_dir}/.." && pwd)"
    declare -g out_dir="${repo_dir}/tests/full/inputs"

    declare -g blk_url="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz"
    declare -g ref_url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr22.fa.gz"


    declare -g fa_head_lines=4000000
    declare -g nreads=60000

    declare -g keep_bed="${out_dir}/chr22.keep.bed"
    declare -g out_blk="${out_dir}/hg38-blacklist.v2.bed.gz"
    declare -g out_fa="${out_dir}/chr22-test.fa.gz"
    declare -g excl_bed="${out_dir}/chr22.exclude.blacklist.bed.gz"

    declare -g RUNS=( SRR3819940 SRR3819939 SRR3819938 SRR3819937 )

}

main(){
    parse_args "$@"
    clean_inputs_dir
    ensure_gitignore
    for acc in "${RUNS[@]}"; do
	fetch_cfdna_reads
    done
    fetch_chr22_ref_subset
    fetch_blacklist
    make_exclude_keep_beds
}

clean_inputs_dir() {
  if [[ -d "$out_dir" ]]; then
    echo "[clean] removing ${out_dir}"
    rm -rf "$out_dir"
  fi
  mkdir -p "$out_dir"
}

ensure_gitignore() {
  gi="${repo_dir}/.gitignore"

  req1='tests/full/*'
  req2='!tests/full/inputs/'
  req3='!tests/full/inputs/**'

  touch "$gi"

  have1=0; have2=0; have3=0
  grep -Fxq "$req1" "$gi" && have1=1 || true
  grep -Fxq "$req2" "$gi" && have2=1 || true
  grep -Fxq "$req3" "$gi" && have3=1 || true

  if (( ! have1 || ! have2 || ! have3 )); then
    {
      echo ''
      echo '# --- test data (auto-managed by get_test_data.sh) ---'
      (( have1 )) || echo "$req1"
      (( have2 )) || echo "$req2"
      (( have3 )) || echo "$req3"
    } >> "$gi"
    echo "[gitignore] ensured tests/full ignores with inputs preserved"
  fi
}

fetch_cfdna_reads() {
  mkdir -p "$out_dir"

  ( cd "$out_dir" && fastq-dump --split-files --gzip -X "$nreads" "$acc" )

  echo "R1 reads: $(zcat "${out_dir}/${acc}_1.fastq.gz" | wc -l | awk '{print $1/4}')"
  echo "R2 reads: $(zcat "${out_dir}/${acc}_2.fastq.gz" | wc -l | awk '{print $1/4}')"
}

fetch_chr22_ref_subset() {
  mkdir -p "$out_dir"

  echo "[ref] chr22 subset → ${out_fa}"
  ( set +o pipefail
    wget -qO- "$ref_url" | zcat 2>/dev/null | head -n "$fa_head_lines" | gzip > "$out_fa"
  ) || true

  test -s "$out_fa"
}

fetch_blacklist() {
  mkdir -p "$out_dir"

  echo "[ref] hg38 blacklist → ${out_blk}"
  wget -qO "$out_blk" "$blk_url"
  [[ -s "$out_blk" ]] || { echo "ERR: failed to write ${out_blk}" >&2; exit 1; }
}

make_exclude_keep_beds() {
  mkdir -p "$out_dir"

  echo "[beds] building keep/exclude from ${out_fa}"

  fa_ungz="${out_dir}/chr22.test.fa"
  zcat "$out_fa" > "$fa_ungz"
  test -s "$fa_ungz"

  samtools faidx "$fa_ungz"
  test -s "${fa_ungz}.fai"

  # KEEP: spans for all contigs present in the .fai (chr22-only here)
  awk 'BEGIN{OFS="\t"} {print $1,0,$2}' "${fa_ungz}.fai" \
    | sort -k1,1 -k2,2n > "$keep_bed"
  test -s "$keep_bed"

  # EXCLUDE: clip blacklist to contigs present in .fai and bgzip
  awk 'NR==FNR{ok[$1]=1; next} ok[$1]' <(cut -f1 "${fa_ungz}.fai") <(zcat "$out_blk") \
    | sort -k1,1 -k2,2n | bgzip -c > "$excl_bed"
  test -s "$excl_bed"

  tabix -p bed "$excl_bed" || true

  rm -f "$fa_ungz" "${fa_ungz}.fai"

  echo "[beds] keep=${keep_bed}  exclude=${excl_bed}"
}

main "$@"
