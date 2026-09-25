#!/usr/bin/env python3
"""end_motif_summary.py: compare end-motif count files (frag_end_motifs.sh format).

Usage: end_motif_summary.py <count_file>...
Prints one TSV row per file: total ends, OTHER count, CC-end share (sum of the
16 CC.. motifs over all ACGT motifs), the share of the six CC motifs Chan et
al. 2020 report together, five most frequent motifs, and the Pearson r of its
256 motif fractions against the first file.
"""
import argparse

import pandas as pd

# The six CC motifs Chan et al. 2020 (Am J Hum Genet 107:882) report together
# for healthy plasma: 10.64 % (range 9.76-11.69 %) of end motifs.
CHAN6 = ["CCCA", "CCAG", "CCTG", "CCAA", "CCCT", "CCAT"]


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("count_files", nargs="+", help="frag_end_motifs.sh count files")
    return p.parse_args(argv)


def summarize(files):
    """One row per file; Pearson r is against the first file's fractions."""
    fracs = {}
    rows = []
    for f in files:
        d = pd.read_csv(f, sep="\t", header=None, index_col=0)[1]
        m = d.drop("OTHER")
        fr = m / m.sum()
        fracs[f] = fr
        cc = fr[fr.index.str.startswith("CC")].sum()
        chan6 = fr[CHAN6].sum()
        top5 = ",".join(fr.sort_values(ascending=False).index[:5])
        r = fr.corr(fracs[files[0]])
        rows.append(
            f"{f}\t{int(d.sum())}\t{int(d['OTHER'])}\t{cc:.4f}\t{chan6:.4f}\t{top5}\t{r:.4f}"
        )
    return rows


def main(argv=None):
    args = parse_args(argv)
    print("file\ttotal_ends\tother\tcc_share\tchan6_share\ttop5\tpearson_r_vs_first")
    for row in summarize(args.count_files):
        print(row)


if __name__ == "__main__":
    main()
