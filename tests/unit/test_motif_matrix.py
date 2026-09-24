"""Tests for scripts/frag_motif_matrix.py: layout, orientation, values."""

import subprocess
import sys
from itertools import product
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "frag_motif_matrix.py"
MOTIFS = ["".join(p) for p in product("ACGT", repeat=4)]


def write_counts(path, counts, other=0):
    lines = [f"{m}\t{counts.get(m, 0)}" for m in MOTIFS] + [f"OTHER\t{other}"]
    path.write_text("\n".join(lines) + "\n")


def run(tmp_path, files):
    c_out, f_out = tmp_path / "counts.tsv", tmp_path / "fractions.tsv"
    subprocess.run(
        [sys.executable, str(SCRIPT), str(c_out), str(f_out), *map(str, files)],
        check=True,
    )
    return pd.read_csv(c_out, sep="\t"), pd.read_csv(f_out, sep="\t")


def two_libraries(tmp_path):
    a = tmp_path / "libA.chrT.motif_counts.tsv"
    b = tmp_path / "libB.chrT.motif_counts.tsv"
    write_counts(a, {"CCCA": 3, "AAAA": 1}, other=5)
    write_counts(b, {"CCTG": 2, "TTTT": 2})
    return [a, b]


def test_orientation(tmp_path):
    counts, fracs = run(tmp_path, two_libraries(tmp_path))
    for df in (counts, fracs):
        assert list(df.columns) == ["motif", "libA", "libB"]
        assert list(df["motif"]) == MOTIFS
        assert len(df) == 256


def test_values(tmp_path):
    counts, fracs = run(tmp_path, two_libraries(tmp_path))
    c = counts.set_index("motif")
    f = fracs.set_index("motif")
    assert c.loc["CCCA", "libA"] == 3 and c.loc["AAAA", "libA"] == 1
    assert c["libA"].sum() == 4  # OTHER is not a motif row
    assert f.loc["CCCA", "libA"] == 0.75
    assert f.loc["CCTG", "libB"] == 0.5
    assert abs(f.sum() - 1).max() < 1e-12


def test_library_id_with_dot(tmp_path):
    p = tmp_path / "kh_01.rep.2.ncbi_hg38.motif_counts.tsv"
    write_counts(p, {"CCCA": 1})
    counts, _ = run(tmp_path, [p])
    assert list(counts.columns) == ["motif", "kh_01.rep.2"]


def test_library_with_no_ends_fails(tmp_path):
    p = tmp_path / "empty.chrT.motif_counts.tsv"
    write_counts(p, {}, other=3)
    res = subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            str(tmp_path / "c.tsv"),
            str(tmp_path / "f.tsv"),
            str(p),
        ],
        capture_output=True,
        text=True,
    )
    assert res.returncode != 0
    assert "empty" in res.stderr
