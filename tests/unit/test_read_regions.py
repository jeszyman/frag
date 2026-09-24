"""Known-answer test for scripts/frag_read_regions.sh."""

import shutil
import subprocess
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "frag_read_regions.sh"

pytestmark = pytest.mark.skipif(
    not shutil.which("bedtools"), reason="bedtools must be on PATH"
)


def test_autosomes_minus_blacklisted_bases(tmp_path):
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t1000\t6\t60\t61\nchrX\t500\t1030\t60\t61\n")
    blk = tmp_path / "blacklist.bed"
    blk.write_text(
        "chr1\t100\t200\tHigh Signal Region\nchrX\t10\t20\tLow Mappability\n"
    )
    out = tmp_path / "regions.bed"
    subprocess.run(["bash", str(SCRIPT), str(fai), str(blk), str(out)], check=True)
    assert out.read_text() == "chr1\t0\t100\nchr1\t200\t1000\n"


def test_no_shared_chromosome_fails(tmp_path):
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t1000\t6\t60\t61\n")
    blk = tmp_path / "blacklist.bed"
    blk.write_text("1\t100\t200\tHigh Signal Region\n")
    res = subprocess.run(
        ["bash", str(SCRIPT), str(fai), str(blk), str(tmp_path / "o.bed")],
        capture_output=True,
        text=True,
    )
    assert res.returncode != 0
    assert "share no chromosome" in res.stderr
