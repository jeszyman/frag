"""Per-library count files from extract_5p_motifs.sh -> motif_counts.tsv (motifs x libraries) and
motifs_rel_freq_wide.tsv (libraries x motifs, rows sum to 1); same layout as the frag pipeline's
count_motifs rule, so the two are directly diffable. Library id = file stem.

  python aggregate_motif_counts.py <counts_dir> <out_dir>
"""
import sys
from itertools import product
from pathlib import Path

import pandas as pd

MOTIFS = ["".join(p) for p in product("ACGT", repeat=4)]
src, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
cols, other = {}, {}
for f in sorted(src.glob("*.tsv")):
    d = pd.read_csv(f, sep="\t", header=None, index_col=0)[1]
    cols[f.stem] = d.reindex(MOTIFS).fillna(0).astype(int)
    other[f.stem] = int(d.get("OTHER", 0))
if not cols:
    sys.exit(f"[FAIL] no *.tsv under {src}")
counts = pd.DataFrame(cols).rename_axis("motif")
counts.to_csv(out / "motif_counts.tsv", sep="\t")
(counts / counts.sum()).T.rename_axis("library").sort_index().to_csv(out / "motifs_rel_freq_wide.tsv", sep="\t")
tot = counts.sum()
print(f"{len(cols)} libraries; ends per library median {tot.median()/1e6:.1f}M [{tot.min()/1e6:.1f}, {tot.max()/1e6:.1f}]; "
      f"OTHER max {max(other.values())}", file=sys.stderr)
