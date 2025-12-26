#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# org_update.sh
#
# Orchestrate a small repo maintenance workflow:
#   1) Save open Emacs buffers via emacsclient (if a server is available)
#   2) Tangle selected repo Org files
#   3) Export a README.md from a specific Org node (via a Python helper)
#   4) Save buffers again (best-effort)
#   5) Run a git update workflow in each repo (via basecamp_functions.sh)
#
# Diagnostics are written to stderr. Stdout is reserved for primary outputs.
# -----------------------------------------------------------------------------

# =============================================================================
# SECTION: CONFIGURATION
# =============================================================================

BASECAMP_LIB="${BASECAMP_LIB:-${HOME}/repos/basecamp/lib/basecamp_functions.sh}"

REPOS=(
  frag
)

REPO_DIR="${HOME}/repos/frag"
README_ORG="${REPO_DIR}/frag.org"
README_NODE="339b69f6-6c09-4e9d-a2ec-27bdf5747163"
README_EXPORT="${HOME}/repos/emacs/scripts/emacs_export_header_to_markdown.py"

# =============================================================================
# SECTION: USAGE AND INPUT PARSING
# =============================================================================

print_usage() {
  cat <<EOF
Usage: ${0##*/} [OPTIONS]

Run a repo maintenance workflow: save Emacs buffers, tangle Org files, export
README.md, and run a git update workflow.

Options:
  -h, --help  Show this help message and exit

Examples:
  ${0##*/}

EOF
}

parse_args() {
  if [[ $# -gt 0 ]]; then
    case "${1:-}" in
      -h|--help)
        print_usage
        exit 0
        ;;
      *)
        echo "Error: Unrecognized argument: ${1:-}" >&2
        echo >&2
        print_usage >&2
        exit 2
        ;;
    esac
  fi

  if [[ ! -r "$BASECAMP_LIB" ]]; then
    echo "Error: BASECAMP_LIB not readable: $BASECAMP_LIB" >&2
    exit 1
  fi

  if [[ ! -d "$REPO_DIR" ]]; then
    echo "Error: REPO_DIR not found: $REPO_DIR" >&2
    exit 1
  fi

  if [[ ! -f "$README_ORG" ]]; then
    echo "Error: README_ORG not found: $README_ORG" >&2
    exit 1
  fi

  if [[ ! -f "$README_EXPORT" ]]; then
    echo "Error: README_EXPORT script not found: $README_EXPORT" >&2
    exit 1
  fi
}


################
###   Main   ###
################


main() {
  parse_args "$@"

  source "$BASECAMP_LIB"

  save_in_emacs

  for repo in "${REPOS[@]}"; do
    tangle_repo_org "$repo" || true
  done

  update_readme
  save_in_emacs

  for repo in "${REPOS[@]}"; do
    git_update_repo "$repo" || true
  done
}

# =============================================================================
# SECTION: WORK FUNCTIONS
# =============================================================================

emacs_server_available() {
  command -v emacsclient >/dev/null 2>&1 && emacsclient -e "(progn t)" >/dev/null 2>&1
}

save_in_emacs() {
  if emacs_server_available; then
    emacsclient -e "(save-some-buffers t)" >/dev/null
    echo "[INFO] Saved buffers via emacsclient" >&2
  else
    echo "[INFO] No Emacs server; skipping save" >&2
  fi
}

update_readme() {
  echo "[INFO] Exporting README.md from: $README_ORG" >&2

  pushd "$REPO_DIR" >/dev/null
  python3 "$README_EXPORT" --org_file "$README_ORG" --node_id "$README_NODE"
  popd >/dev/null
}

tangle_repo_org() {
  local repo="$1"
  local org_file="${HOME}/repos/${repo}/${repo}.org"

  echo "[INFO] Repo: $repo" >&2

  if [[ -f "$org_file" ]]; then
    echo "[INFO] Tangling: $org_file" >&2
    tangle "$org_file"
  else
    echo "[WARN] Missing Org file; skipping: $org_file" >&2
  fi
}

git_update_repo() {
  local repo="$1"
  local repo_dir="${HOME}/repos/${repo}"
  local output=""

  echo "[INFO] Updating git workflow: $repo" >&2

  if [[ ! -d "$repo_dir" ]]; then
    echo "[WARN] Repo directory missing; skipping: $repo_dir" >&2
    return 0
  fi

  pushd "$repo_dir" >/dev/null

  if ! output="$(git_wkflow_up 2>&1)"; then
    echo "[ERROR] git_wkflow_up failed in: $repo" >&2
    echo "$output" >&2
    popd >/dev/null
    return 1
  fi

  printf '%s\n%s\n' "$(date)" "$output" >&2
  popd >/dev/null
}

# =============================================================================
# SECTION: ENTRY POINT
# =============================================================================

main "$@"
