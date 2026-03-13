#!/usr/bin/env python3
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/frag/frag.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-13 16:35:03
# ============================================================

"""NMF decomposition of fragment length frequency matrix.

Produces W (samples x components) and H (components x lengths) matrices.
Uses KL divergence with multiplicative update solver (Renaud 2022).
"""
import sys
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF

freq_csv = sys.argv[1]
n_components = int(sys.argv[2])
out_w = sys.argv[3]
out_h = sys.argv[4]

freq_mat = pd.read_csv(freq_csv, index_col=0)

nmf = NMF(
    n_components=n_components,
    init="random",
    random_state=42,
    solver="mu",
    beta_loss="kullback-leibler",
    max_iter=5000,
)
W = nmf.fit_transform(freq_mat.values)
H = nmf.components_

W_df = pd.DataFrame(
    W,
    index=freq_mat.index,
    columns=[f"nmf_{i+1}" for i in range(n_components)],
)
H_df = pd.DataFrame(
    H,
    index=[f"nmf_{i+1}" for i in range(n_components)],
    columns=freq_mat.columns,
)

W_df.to_csv(out_w)
H_df.to_csv(out_h)
