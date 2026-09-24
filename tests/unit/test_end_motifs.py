"""Known-answer tests for scripts/frag_end_motifs.sh.

A 400 bp synthetic contig and a hand-written SAM give exact expected counts.
The expected motifs are computed here by slicing the contig directly, which
is independent of the script's samtools/bedtools route.
"""

import random
import shutil
import subprocess
from collections import Counter
from itertools import product
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "frag_end_motifs.sh"
MOTIFS = ["".join(p) for p in product("ACGT", repeat=4)]
READ_LEN = 30

pytestmark = pytest.mark.skipif(
    not (shutil.which("samtools") and shutil.which("bedtools")),
    reason="samtools and bedtools must be on PATH",
)


def make_ref():
    rng = random.Random(7)
    seq = [rng.choice("ACGT") for _ in range(400)]
    seq[300:304] = list("NNNN")  # the 5' 4-mer of a read starting at 301 contains N
    return "".join(seq)


REF = make_ref()
COMP = str.maketrans("ACGTN", "TGCAN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def plus_motif(pos1):
    """5' 4-mer of a plus-strand read whose leftmost aligned base is 1-based pos1."""
    return REF[pos1 - 1 : pos1 + 3]


def minus_motif(pos1, aligned_len=READ_LEN):
    """5' 4-mer of a minus-strand read: last 4 aligned bases, reverse-complemented."""
    end = pos1 - 1 + aligned_len
    return revcomp(REF[end - 4 : end])


def rec(name, flag, pos, mapq=60, cigar="30M", pnext=0):
    return f"{name}\t{flag}\tchrT\t{pos}\t{mapq}\t{cigar}\t=\t{pnext}\t0\t*\t*"


def pair(name, f1, p1, f2, p2, mapq=60, cigar1="30M"):
    return [rec(name, f1, p1, mapq, cigar1, p2), rec(name, f2, p2, mapq, "30M", p1)]


# Each entry: SAM lines, and the ends the definition counts from them.
CASES = {
    "r1_plus": (pair("r1_plus", 99, 11, 147, 151), [plus_motif(11), minus_motif(151)]),
    # Read 1 on the minus strand: the case read-number keying gets wrong.
    "r1_minus": (
        pair("r1_minus", 83, 201, 163, 61),
        [minus_motif(201), plus_motif(61)],
    ),
    "duplicate": (pair("dup", 99 + 1024, 21, 147 + 1024, 121), []),
    "secondary": ([rec("sec", 99 + 256, 31, pnext=131)], []),
    "supplementary": ([rec("supp", 99 + 2048, 41, pnext=141)], []),
    "mate_unmapped": (
        [rec("mu", 73, 51, pnext=51), rec("mu", 133, 51, mapq=0, cigar="*", pnext=51)],
        [],
    ),
    "mapq29": (pair("lowq", 99, 71, 147, 171, mapq=29), []),
    # QC-fail (0x200) reads are kept, as in the reference implementation.
    "qcfail": (
        pair("qcf", 99 + 512, 131, 147 + 512, 260),
        [plus_motif(131), minus_motif(260)],
    ),
    # Soft clip: the end is the first aligned base, not the clipped read start.
    "softclip": (
        pair("sc", 99, 101, 147, 250, cigar1="5S25M"),
        [plus_motif(101), minus_motif(250)],
    ),
    "n_end": (pair("nend", 99, 301, 147, 351), ["OTHER", minus_motif(351)]),
}


def write_fixture(tmp_path, sam_lines):
    fa = tmp_path / "ref.fa"
    fa.write_text(">chrT\n" + REF + "\n")
    subprocess.run(["samtools", "faidx", str(fa)], check=True)
    sam = tmp_path / "in.sam"
    sam.write_text(
        "@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:chrT\tLN:400\n" + "\n".join(sam_lines) + "\n"
    )
    bam = tmp_path / "in.bam"
    subprocess.run(["samtools", "sort", "-o", str(bam), str(sam)], check=True)
    subprocess.run(["samtools", "index", str(bam)], check=True)
    return bam, fa


def run_script(tmp_path, bam, fa, max_ends=0, seed=42, name="out.tsv"):
    out = tmp_path / name
    subprocess.run(
        [
            "bash",
            str(SCRIPT),
            str(bam),
            str(fa),
            str(out),
            "1",
            str(max_ends),
            str(seed),
        ],
        check=True,
    )
    rows = [line.split("\t") for line in out.read_text().splitlines()]
    return rows, {m: int(c) for m, c in rows}


@pytest.fixture
def all_cases(tmp_path):
    lines, expected = [], Counter()
    for sam_lines, ends in CASES.values():
        lines += sam_lines
        expected.update(ends)
    bam, fa = write_fixture(tmp_path, lines)
    return tmp_path, bam, fa, expected


def test_output_layout(all_cases):
    tmp_path, bam, fa, _ = all_cases
    rows, _ = run_script(tmp_path, bam, fa)
    assert [r[0] for r in rows] == MOTIFS + ["OTHER"]


def test_known_answer_counts(all_cases):
    tmp_path, bam, fa, expected = all_cases
    _, got = run_script(tmp_path, bam, fa)
    want = {m: expected.get(m, 0) for m in MOTIFS + ["OTHER"]}
    assert got == want


def test_read1_minus_takes_fragment_5prime_end(tmp_path):
    """Read-number keying takes read 1's leftmost 4 bases and read 2's rightmost 4
    bases reverse-complemented; for a read-1-minus pair both are read 3' ends,
    interior to the fragment. Neither may appear in the output."""
    bam, fa = write_fixture(tmp_path, CASES["r1_minus"][0])
    _, got = run_script(tmp_path, bam, fa)
    true_ends = Counter([minus_motif(201), plus_motif(61)])
    interior = {REF[200:204], revcomp(REF[86:90])}
    assert not interior & set(
        true_ends
    ), "fixture seed makes an interior motif equal a true end"
    assert {m: c for m, c in got.items() if c} == dict(true_ends)
    assert all(got[m] == 0 for m in interior)


def many_pairs(n, seed=3):
    rng = random.Random(seed)
    lines = []
    for i in range(n):
        p1 = rng.randint(1, 170)
        p2 = rng.randint(p1, 340)
        lines += pair(f"q{i}", 99, p1, 147, p2)
    return lines


def test_max_ends_zero_counts_all(tmp_path):
    bam, fa = write_fixture(tmp_path, many_pairs(500))
    _, got = run_script(tmp_path, bam, fa, max_ends=0)
    assert sum(got.values()) == 1000


def test_max_ends_above_total_equals_all(tmp_path):
    bam, fa = write_fixture(tmp_path, many_pairs(500))
    _, all_ends = run_script(tmp_path, bam, fa, max_ends=0, name="a.tsv")
    _, capped = run_script(tmp_path, bam, fa, max_ends=10**9, name="b.tsv")
    assert capped == all_ends


def test_subsample_keeps_mates_together_and_is_seeded(tmp_path):
    bam, fa = write_fixture(tmp_path, many_pairs(2000))
    _, a = run_script(tmp_path, bam, fa, max_ends=1000, seed=42, name="a.tsv")
    _, b = run_script(tmp_path, bam, fa, max_ends=1000, seed=42, name="b.tsv")
    total = sum(a.values())
    assert 0 < total < 4000
    assert total % 2 == 0  # every counted pair contributes both mates
    assert a == b
