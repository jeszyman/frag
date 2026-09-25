"""tools/get_test_data.sh cleanup removes only the files the script rebuilds."""

import shutil
import subprocess
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[2] / "tools" / "get_test_data.sh"
REBUILT = [
    "SRR3819937_1.fastq.gz",
    "chr22-test.fa.gz",
    "hg38-blacklist.v2.bed.gz",
    "chr22.exclude.blacklist.bed.gz",
    "chr22.exclude.blacklist.bed.gz.tbi",
]
KEPT = ["cytoBand.chr22.txt", "delfi_bins.chr22.tsv"]


def test_cleanup_keeps_fixtures_it_does_not_rebuild(tmp_path):
    (tmp_path / "tools").mkdir()
    shutil.copy(SCRIPT, tmp_path / "tools" / "get_test_data.sh")
    inputs = tmp_path / "tests" / "full" / "inputs"
    inputs.mkdir(parents=True)
    for name in REBUILT + KEPT:
        (inputs / name).write_text("x\n")
    subprocess.run(
        [
            "bash",
            "-c",
            "source tools/get_test_data.sh && parse_args && clean_inputs_dir",
        ],
        cwd=tmp_path,
        check=True,
        timeout=30,
    )
    assert sorted(p.name for p in inputs.iterdir()) == sorted(KEPT)
