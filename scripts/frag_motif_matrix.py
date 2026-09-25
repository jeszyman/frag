#!/usr/bin/env python3
"""frag_motif_matrix.py: per-library end-motif count files to motif x library matrices.

Adapted from aggregate_motif_counts.py by Axel Hidalgo, commit 7a5767ab1
(2026-09-18). Unlike that
script, both outputs keep motifs as rows and libraries as columns, the layout
frag_motif_diversity.R and frag_fprofiles.py read.

Usage:
    frag_motif_matrix.py <counts_out.tsv> <fractions_out.tsv> <count_file>...

Each count file is frag_end_motifs.sh output named {library_id}.{ref_name}.motif_counts.tsv.
counts_out: integer counts, 256 ACGT 4-mers. fractions_out: counts divided by each
library's total over those 256 motifs (the OTHER row is excluded).
"""
import argparse
import os
import re
import sys
from itertools import product

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from frag_checks import check_n, log_n  # noqa: E402

MOTIFS = ["".join(p) for p in product("ACGT", repeat=4)]


def library_id(path):
    """Library id: file name without the trailing .{ref_name}.motif_counts.tsv."""
    return re.sub(r"\.[^.]+\.motif_counts\.tsv$", "", os.path.basename(path))


def read_counts(path):
    """One library's counts over the 256 ACGT motifs, in MOTIFS order."""
    d = pd.read_csv(
        path, sep="\t", header=None, names=["motif", "count"], index_col="motif"
    )
    check_n(d, 257, f"{path} rows")
    counts = d["count"].reindex(MOTIFS)
    if counts.isna().any() or counts.sum() == 0:
        raise ValueError(f"[check] {path}: empty or missing motif rows")
    return counts.astype(int)


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("counts_out", help="motif x library count matrix (TSV)")
    p.add_argument("fractions_out", help="motif x library fraction matrix (TSV)")
    p.add_argument("count_files", nargs="+", help="per-library count files")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    counts = pd.DataFrame({library_id(f): read_counts(f) for f in args.count_files})
    counts.index.name = "motif"
    log_n(counts, "motif count matrix rows")
    counts.to_csv(args.counts_out, sep="\t")
    (counts / counts.sum()).to_csv(args.fractions_out, sep="\t")


if __name__ == "__main__":
    main()
