#!/usr/bin/env python3
"""make_delfi_bins.py: DELFI 5 Mb bins on an hg38 reference.

Applies the bin rule of Mathios et al. 2021 to hg38: chr1-chr22 are tiled
into non-overlapping bins from position 0 (the last, partial bin of each
chromosome is dropped); each bin gets its GC fraction over non-N bases, its
mean mappability over all bases (bases with no bigWig value count as 0), its
blacklisted bases and its chromosome arm; bins with GC < min-gc or
mappability < min-map are dropped.

Output TSV with a header, 0-based half-open coordinates:
    chr start end arm gc map blacklisted_bases

Usage:
    make_delfi_bins.py --fasta hg38.fa --blacklist hg38-blacklist.v2.bed.gz \
        --mappability k100.Umap.MultiTrackMappability.bw \
        --cytoband cytoBand.txt --out delfi_bins_hg38.tsv
"""
import argparse
import gzip
import io
import re
import subprocess
import sys
import tempfile

import pandas as pd
import pyBigWig

# =============================================================================
# ARGUMENTS
# =============================================================================


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--fasta", required=True, help="reference FASTA with a .fai index")
    p.add_argument("--blacklist", required=True, help="blacklist BED (.bed or .bed.gz)")
    p.add_argument("--mappability", required=True, help="mappability bigWig")
    p.add_argument("--cytoband", required=True, help="UCSC cytoBand table")
    p.add_argument("--out", required=True, help="output TSV")
    p.add_argument("--bin-size", type=int, default=5_000_000)
    p.add_argument("--min-gc", type=float, default=0.3)
    p.add_argument("--min-map", type=float, default=0.9)
    return p.parse_args(argv)


# =============================================================================
# WORK FUNCTIONS
# =============================================================================


def tile_autosomes(fai, bin_size):
    """Full-size bins over every chr1-chr22 contig of the index, from position 0."""
    rows = []
    with open(fai) as fh:
        for line in fh:
            chrom, length = line.split("\t")[:2]
            if re.fullmatch(r"chr[0-9]+", chrom):
                for i in range(int(length) // bin_size):
                    rows.append((chrom, i * bin_size, (i + 1) * bin_size))
    if not rows:
        sys.exit(f"[FAIL] no chr1-chr22 contig of at least one bin in {fai}")
    return pd.DataFrame(rows, columns=["chr", "start", "end"])


def gc_non_n(fasta, bins, bed_path):
    """GC fraction over A, C, G, T bases of each bin (bedtools nuc counts).

    A bin with no A, C, G or T base (all N) gets NaN and fails the GC filter.
    """
    out = subprocess.run(
        ["bedtools", "nuc", "-fi", fasta, "-bed", bed_path],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    nuc = pd.read_csv(io.StringIO(out), sep="\t")
    a, c, g, t = (nuc.iloc[:, i] for i in (5, 6, 7, 8))
    return ((c + g) / (a + c + g + t)).to_numpy()


def mean_mappability(bigwig, bins):
    """Mean bigWig value over all bases of each bin; bases without a value count as 0."""
    bw = pyBigWig.open(bigwig)
    values = []
    for chrom, start, end in bins[["chr", "start", "end"]].itertuples(index=False):
        mean = bw.stats(chrom, start, end, type="mean", exact=True)[0]
        covered = bw.stats(chrom, start, end, type="coverage", exact=True)[0]
        values.append(0.0 if mean is None else mean * covered)
    bw.close()
    return values


def blacklisted_bases(blacklist, bed_path, tmpdir):
    """Bases of each bin covered by the merged blacklist."""
    opener = gzip.open if blacklist.endswith(".gz") else open
    with opener(blacklist, "rt") as fh:
        blk = pd.read_csv(fh, sep="\t", header=None, usecols=[0, 1, 2])
    blk_path = f"{tmpdir}/blacklist.bed"
    blk.sort_values([0, 1]).to_csv(blk_path, sep="\t", header=False, index=False)
    merged = subprocess.run(
        ["bedtools", "merge", "-i", blk_path],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    merged_path = f"{tmpdir}/blacklist.merged.bed"
    with open(merged_path, "w") as fh:
        fh.write(merged)
    out = subprocess.run(
        ["bedtools", "intersect", "-a", bed_path, "-b", merged_path, "-wao"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout
    hits = pd.read_csv(io.StringIO(out), sep="\t", header=None)
    return hits.groupby([0, 1, 2], sort=False)[hits.columns[-1]].sum().to_numpy()


def arms(cytoband, bins):
    """p when the bin midpoint lies before the end of the chromosome's last p band."""
    cyto = pd.read_csv(cytoband, sep="\t", header=None, usecols=[0, 1, 2, 3])
    cyto.columns = ["chr", "start", "end", "band"]
    p_end = (
        cyto[cyto["band"].astype(str).str.startswith("p")].groupby("chr")["end"].max()
    )
    mid = (bins["start"] + bins["end"]) / 2
    arm = [
        chrom[3:] + ("p" if m < p_end.get(chrom, 0) else "q")
        for chrom, m in zip(bins["chr"], mid)
    ]
    return arm


def main(argv=None):
    args = parse_args(argv)
    bins = tile_autosomes(args.fasta + ".fai", args.bin_size)
    print(f"[make_delfi_bins] {len(bins)} full bins on chr1-chr22", file=sys.stderr)
    with tempfile.TemporaryDirectory() as tmpdir:
        bed_path = f"{tmpdir}/bins.bed"
        bins.to_csv(bed_path, sep="\t", header=False, index=False)
        bins["arm"] = arms(args.cytoband, bins)
        bins["gc"] = gc_non_n(args.fasta, bins, bed_path)
        bins["map"] = mean_mappability(args.mappability, bins)
        bins["blacklisted_bases"] = blacklisted_bases(args.blacklist, bed_path, tmpdir)
    keep = bins[(bins["gc"] >= args.min_gc) & (bins["map"] >= args.min_map)]
    keep.to_csv(args.out, sep="\t", index=False, float_format="%.6f")
    print(
        f"[make_delfi_bins] kept {len(keep)} of {len(bins)} bins "
        f"({len(keep) * args.bin_size / 1e9:.3f} Gb) in {args.out}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
