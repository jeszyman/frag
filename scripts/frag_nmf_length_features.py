#!/usr/bin/env python3
"""NMF decomposition of fragment length frequency matrix.

Produces W (samples x components) and H (components x lengths) matrices.
Uses KL divergence with multiplicative update solver (Renaud 2022).
"""
import os
import sys
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from frag_checks import check_nonempty, log_n

freq_csv = sys.argv[1]
n_components = int(sys.argv[2])
out_w = sys.argv[3]
out_h = sys.argv[4]

freq_mat = pd.read_csv(freq_csv, index_col=0)
check_nonempty(freq_mat, "freq_mat")
log_n(freq_mat, "freq matrix")

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
