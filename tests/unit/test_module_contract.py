"""workflows/frag.smk lists in its preamble every wrapper variable it uses."""

import re
from pathlib import Path

SMK = Path(__file__).resolve().parents[2] / "workflows" / "frag.smk"
NAME = re.compile(
    r"\b(D_[A-Z_]+|FRAG_[A-Z_]+|CONDA_[A-Z_]+|NMF_[A-Z_]+|frag_ref_names)\b"
)


def test_preamble_lists_every_wrapper_variable():
    lines = SMK.read_text().splitlines()
    header = "\n".join(l for l in lines if l.startswith("#"))
    body = "\n".join(l for l in lines if not l.lstrip().startswith("#"))
    used = set(NAME.findall(body))
    listed = set(re.findall(r"^#\s+(\w+)\s", header, flags=re.M))
    assert used, "no wrapper variables found"
    assert sorted(used - listed) == []
