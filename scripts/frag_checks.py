#!/usr/bin/env python3
"""frag_checks.py — data-integrity helpers for frag analysis scripts.

Import with a sibling-path shim (rules run scripts as `python3 scripts/<x>.py`
from the repo root):

    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from frag_checks import check_nonempty, log_n
"""
import sys


def log_n(df, label="rows"):
    """Log row count to stderr; return df for chaining."""
    print(f"[check] {label}: {len(df)}", file=sys.stderr)
    return df


def check_n(df, n, label="rows"):
    """Raise unless df has exactly n rows."""
    if len(df) != n:
        raise ValueError(f"[check] {label}: expected {n}, got {len(df)}")
    return df


def check_nonempty(df, label="object"):
    """Raise if df has zero rows."""
    if len(df) == 0:
        raise ValueError(f"[check] {label} is empty")
    return df


def assert_cols(df, cols, label="object"):
    """Raise unless all `cols` are present in df.columns."""
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise ValueError(f"[check] {label} missing columns: {', '.join(missing)}")
    return df
