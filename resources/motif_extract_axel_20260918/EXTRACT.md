# End-motif extraction, strand-aware — for comparing against an existing motif matrix

Two files. Needs `samtools` and `bedtools` on PATH (or `TOOLS_BIN=/path/to/bin`), Python 3 + pandas
for the aggregator, and the same reference FASTA (+ `.fai`) the BAMs were aligned to.

```bash
# one library (3 threads); ~1 h per 30 GB BAM, embarrassingly parallel across libraries
./extract_5p_motifs.sh SAMPLE.bam ref.fa counts/SAMPLE.tsv 3

# all libraries, 16 at a time
ls bams/*.bam | xargs -P 16 -I{} bash -c 'b={}; ./extract_5p_motifs.sh $b ref.fa counts/$(basename ${b%.bam}).tsv 3'

# -> motif_counts.tsv (motifs x libraries) and motifs_rel_freq_wide.tsv (libraries x motifs)
python aggregate_motif_counts.py counts/ out/
```

**Definition.** For every paired, primary, mapped, non-duplicate read with MAPQ ≥ 30 (both mates), the
4 reference bases at the read's **5′ end**: a plus-strand read gives `[start, start+4)`, a minus-strand
read gives `[end-4, end)` reverse-complemented (`bedtools getfasta -s` does the flip). Non-ACGT 4-mers
are dropped and counted in the `OTHER` row. Reference-based, so it is correct for bisulfite / EM-seq
alignments where the read sequence is converted.

**What it is not.** Keying the end on read 1 vs read 2 instead of strand (take the left edge of R1, the
right edge of R2) returns the read's 3′ end for the half of fragments where R1 maps to the minus strand;
that end lies inside the fragment whenever the fragment is longer than the read. Measured on one EM-seq
library: 50 % of R1 minus-strand, 80 % of fragments longer than the read, 40 % of motifs interior. The
signature of that mix is AAAA / TTTT at the top of the frequency table instead of CCCA / CCTG / CCAG.

**Comparing to an existing matrix.** Same layout, so per library: Pearson r over the 256 motifs, total
ends (should match to within the OTHER count if the read filter is the same), and the CC-end share
(sum of the 16 `CC..` motifs; ~0.11–0.12 in plasma cfDNA under this definition).
