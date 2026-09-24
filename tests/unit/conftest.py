"""pytest setup: put scripts/ on sys.path so tests can import frag's Python scripts."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
