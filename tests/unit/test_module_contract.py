"""workflows/frag.smk: the preamble lists exactly the wrapper variables the
module uses, and every script runs through the repository root R_FRAG, so the
module works when included from another directory."""

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SMK = REPO / "workflows" / "frag.smk"
NAME = re.compile(
    r"\b(R_[A-Z_]+|D_[A-Z_]+|FRAG_[A-Z_]+|CONDA_[A-Z_]+|NMF_[A-Z_]+|frag_ref_names)\b"
)
LINES = SMK.read_text().splitlines()
HEADER = "\n".join(l for l in LINES if l.startswith("#"))
BODY = "\n".join(l for l in LINES if not l.lstrip().startswith("#"))
LISTED = set(re.findall(r"^#\s+([A-Za-z_]+)\s{2,}", HEADER, flags=re.M))


def test_preamble_lists_every_wrapper_variable():
    used = set(NAME.findall(BODY))
    assert used, "no wrapper variables found"
    assert sorted(used - LISTED) == []


def test_preamble_lists_no_unused_variable():
    used = set(NAME.findall(BODY))
    assert sorted(LISTED - used) == []


def test_scripts_run_through_repo_root():
    calls = re.findall(r"\b(?:bash|Rscript|python3?)\s+(\S*scripts/\S+)", BODY)
    assert calls, "no script calls found"
    assert [c for c in calls if not c.startswith("{R_FRAG}/scripts/")] == []


def test_r_scripts_source_helpers_from_their_own_directory():
    bad = [
        p.name
        for p in (REPO / "scripts").glob("*.R")
        if 'source("scripts/' in p.read_text()
    ]
    assert bad == []
