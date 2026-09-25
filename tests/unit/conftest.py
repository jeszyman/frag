"""pytest setup: put scripts/ on sys.path so tests can import frag's Python scripts.

Tests that need samtools and bedtools skip when those are not on PATH. With
FRAG_REQUIRE_TOOLS=1 (set in CI) a missing tool stops the run instead, so a
broken environment cannot pass by skipping.
"""

import os
import shutil
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))

REQUIRED_TOOLS = ("samtools", "bedtools")


def pytest_configure(config):
    if os.environ.get("FRAG_REQUIRE_TOOLS") == "1":
        missing = [t for t in REQUIRED_TOOLS if not shutil.which(t)]
        if missing:
            raise pytest.UsageError(
                f"FRAG_REQUIRE_TOOLS=1 but not on PATH: {', '.join(missing)}"
            )
