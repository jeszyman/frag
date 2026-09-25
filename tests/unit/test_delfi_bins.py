"""Known-answer test for scripts/frag_delfi_bins.sh."""

import subprocess
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "frag_delfi_bins.sh"
HEADER = "chr\tstart\tend\tarm\tgc\tmap\tblacklisted_bases\n"


def run(tmp_path, table):
    bins = tmp_path / "bins.tsv"
    bins.write_text(table)
    fai = tmp_path / "ref.fa.fai"
    fai.write_text("chr1\t12000000\t6\t60\t61\n")
    out = tmp_path / "bins.bed"
    res = subprocess.run(
        ["bash", str(SCRIPT), str(bins), str(fai), str(out)],
        capture_output=True,
        text=True,
    )
    return res, out


def test_bins_on_reference_contigs(tmp_path):
    res, out = run(
        tmp_path,
        HEADER
        + "chr1\t0\t5000000\t1p\t0.41\t0.95\t0\n"
        + "chr1\t5000000\t10000000\t1p\t0.42\t0.96\t12\n"
        + "chr1\t10000000\t15000000\t1q\t0.43\t0.97\t0\n"
        + "chr2\t0\t5000000\t2p\t0.44\t0.98\t0\n",
    )
    assert res.returncode == 0, res.stderr
    assert out.read_text() == "chr1\t0\t5000000\t0.41\nchr1\t5000000\t10000000\t0.42\n"


def test_wrong_header_fails(tmp_path):
    res, _ = run(tmp_path, "chr\tstart\tend\tgc\nchr1\t0\t5000000\t0.41\n")
    assert res.returncode != 0
    assert "expected chr start end arm gc map blacklisted_bases" in res.stderr
