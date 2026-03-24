#!/usr/bin/env python3
# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/frag/frag.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-24 13:15:21
# ============================================================

"""F-profile decomposition of end motif frequencies.

NMF learns latent motif signatures, then NNLS deconvolves each sample's
motif profile into F-profile contributions (percentages summing to 100).
"""
import sys
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from scipy.optimize import nnls

motif_tsv = sys.argv[1]
n_components = int(sys.argv[2])
out_fprof = sys.argv[3]
out_motif_per_fprof = sys.argv[4]

motifs_df = pd.read_csv(motif_tsv, sep="\t").fillna(0)
M = motifs_df.iloc[:, 1:].astype(float).T  # samples x motifs

n_components = min(n_components, M.shape[0], M.shape[1])

nmf = NMF(
    n_components=n_components,
    init="random",
    random_state=42,
    max_iter=100000,
)
W = nmf.fit_transform(M)
F = nmf.components_  # components x motifs

# NNLS deconvolution per sample
nnls_df = pd.DataFrame()
for lib in motifs_df.columns[1:]:
    sample_freqs = motifs_df[lib].values.astype(float)
    P, _ = nnls(F.T, sample_freqs)
    total = np.sum(P)
    nnls_df[lib] = 100 * P / total if total > 0 else P

nnls_df.index = [f"fprof{i+1}" for i in range(F.shape[0])]
nnls_df.to_csv(out_fprof, sep="\t")

# Normalized motif profiles per F-profile
F_norm = F / F.sum(axis=1, keepdims=True)
F_df = pd.DataFrame(
    F_norm,
    columns=motifs_df.iloc[:, 0],
    index=[f"fprof{i+1}" for i in range(F.shape[0])],
)
F_df.to_csv(out_motif_per_fprof, sep="\t")
