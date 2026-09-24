#!/usr/bin/env python3
"""end_motif_summary.py: compare end-motif count files (frag_end_motifs.sh format).

Usage: end_motif_summary.py <count_file>...
Prints one TSV row per file: total ends, OTHER count, CC-end share (sum of the
16 CC.. motifs over all ACGT motifs), five most frequent motifs, and the
Pearson r of its 256 motif fractions against the first file.
"""
import sys

import pandas as pd

files = sys.argv[1:]
fracs = {}
print("file\ttotal_ends\tother\tcc_share\ttop5\tpearson_r_vs_first")
for f in files:
    d = pd.read_csv(f, sep="\t", header=None, index_col=0)[1]
    m = d.drop("OTHER")
    fr = m / m.sum()
    fracs[f] = fr
    cc = fr[fr.index.str.startswith("CC")].sum()
    top5 = ",".join(fr.sort_values(ascending=False).index[:5])
    r = fr.corr(fracs[files[0]])
    print(f"{f}\t{int(d.sum())}\t{int(d['OTHER'])}\t{cc:.4f}\t{top5}\t{r:.4f}")
