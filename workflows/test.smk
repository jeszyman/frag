# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/frag/frag.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-13 16:35:03
# ============================================================

# CFDNA FRAGMENTOMICS FULL PIPELINE TEST WRAPPER SNAKEFILE
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

import os
import pandas as pd

# ------------------------------------------------------------------------------
# Load YAML Configuration
# ------------------------------------------------------------------------------

configfile: "config/test.yaml"

def resolve_config_paths(config_dict):
    for k, v in config_dict.items():
        if isinstance(v, str):
            config_dict[k] = os.path.expandvars(os.path.expanduser(v))
        elif isinstance(v, dict):
            resolve_config_paths(v)
        elif isinstance(v, list):
            config_dict[k] = [os.path.expandvars(os.path.expanduser(i)) if isinstance(i, str) else i for i in v]

resolve_config_paths(config)

# ------------------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------------------

D_INPUTS    = config["directories"]["inputs"]
D_FRAG      = config["directories"]["frag"]
D_LOGS      = config["directories"]["logs"]
D_BENCHMARK = config["directories"]["benchmark"]

CONDA_FRAG = config["conda"]["frag"]

frag_ref_names = ["chr22"]

FRAG_HEALTHY_LIBRARIES = config.get("healthy_libraries", [])

NMF_N_COMPONENTS = config.get("nmf", {}).get("n_components", 2)

# ------------------------------------------------------------------------------
# Load Tabular Configuration
# ------------------------------------------------------------------------------

class SampleTable:
    def __init__(self, tsv_path, selected_ids):
        df = pd.read_csv(tsv_path, sep="\t")

        missing = sorted(set(selected_ids) - set(df["library_id"]))
        if missing:
            raise ValueError(f"library_id values not found in TSV: {missing}")

        df = df[df["library_id"].isin(selected_ids)].copy()
        self.df = df

    @property
    def frag_library_ids(self):
        return sorted(self.df["library_id"].unique())

    @property
    def r1_map(self):
        return dict(zip(self.df["library_id"], self.df["r1_basename"]))

    @property
    def r2_map(self):
        return dict(zip(self.df["library_id"], self.df["r2_basename"]))

samples = SampleTable(
    tsv_path=config["sample-tsv-path"],
    selected_ids=["lib001", "lib002"],
)

FRAG_LIBRARY_IDS = samples.frag_library_ids

# ------------------------------------------------------------------------------
# Rule all
# ------------------------------------------------------------------------------

rule all:
    input:
        # FASTQs (raw + processed)
        expand(
            f"{D_FRAG}/fastqs/{{library_id}}.{{processing}}_{{read}}.fastq.gz",
            library_id=FRAG_LIBRARY_IDS,
            processing=["raw", "processed"],
            read=["R1", "R2"],
        ),
        # BWA index
        expand(
            f"{D_FRAG}/ref/bwa/{{ref_name}}/{{ref_name}}.fa.sa",
            ref_name=frag_ref_names,
        ),
        # Alignments
        expand(
            f"{D_FRAG}/bams/{{library_id}}.bwa.{{ref_name}}.coorsort.bam",
            library_id=FRAG_LIBRARY_IDS,
            ref_name=frag_ref_names,
        ),
        # GC-filtered bins
        f"{D_FRAG}/ref/keep_5mb.bed",
        # Filtered BAMs
        expand(
            f"{D_FRAG}/bams/{{library_id}}.bwa.{{ref_name}}.filt.bam",
            library_id=FRAG_LIBRARY_IDS,
            ref_name=frag_ref_names,
        ),
        # Fragment BEDs
        expand(
            f"{D_FRAG}/frags/{{library_id}}.{{ref_name}}.frag.bed",
            library_id=FRAG_LIBRARY_IDS,
            ref_name=frag_ref_names,
        ),
        # GC distributions
        expand(
            f"{D_FRAG}/frags/{{library_id}}.{{ref_name}}.gc_distro.csv",
            library_id=FRAG_LIBRARY_IDS,
            ref_name=frag_ref_names,
        ),
        # Healthy GC summary
        expand(
            f"{D_FRAG}/frags/{{ref_name}}.healthy_med.rds",
            ref_name=frag_ref_names,
        ),
        # Sampled fragment BEDs
        expand(
            f"{D_FRAG}/frags/{{library_id}}.{{ref_name}}.sampled_frag.bed",
            library_id=FRAG_LIBRARY_IDS,
            ref_name=frag_ref_names,
        ),
        # Short/long BEDs
        expand(
            f"{D_FRAG}/frags/{{library_id}}.{{ref_name}}.norm_{{length}}.bed",
            library_id=FRAG_LIBRARY_IDS,
            ref_name=frag_ref_names,
            length=["short", "long"],
        ),
        # Fragment counts per bin
        expand(
            f"{D_FRAG}/counts/{{library_id}}.{{ref_name}}.cnt_{{length}}.tmp",
            library_id=FRAG_LIBRARY_IDS,
            ref_name=frag_ref_names,
            length=["short", "long"],
        ),
        # Merged counts
        expand(
            f"{D_FRAG}/frags/{{ref_name}}.frag_counts.tsv",
            ref_name=frag_ref_names,
        ),
        # Ratios
        expand(
            f"{D_FRAG}/frags/{{ref_name}}.ratios.tsv",
            ref_name=frag_ref_names,
        ),
        # End motifs per library
        expand(
            f"{D_FRAG}/motifs/{{library_id}}.{{ref_name}}.motifs.txt",
            library_id=FRAG_LIBRARY_IDS,
            ref_name=frag_ref_names,
        ),
        # Motif matrix
        expand(
            f"{D_FRAG}/motifs/{{ref_name}}.all_motifs.tsv",
            ref_name=frag_ref_names,
        ),
        # Fragment length histograms
        expand(
            f"{D_FRAG}/histograms/{{library_id}}.{{ref_name}}.length_hist.tsv",
            library_id=FRAG_LIBRARY_IDS,
            ref_name=frag_ref_names,
        ),
        # Frequency matrix
        expand(
            f"{D_FRAG}/histograms/{{ref_name}}.count_histogram.csv",
            ref_name=frag_ref_names,
        ),
        expand(
            f"{D_FRAG}/histograms/{{ref_name}}.freq_histogram.csv",
            ref_name=frag_ref_names,
        ),
        # Arm z-scores
        expand(
            f"{D_FRAG}/features/{{ref_name}}.arm_zscores.csv",
            ref_name=frag_ref_names,
        ),
        # Normalized ratios
        expand(
            f"{D_FRAG}/features/{{ref_name}}.ratios_normalized.csv",
            ref_name=frag_ref_names,
        ),
        # QC plots
        expand(
            f"{D_FRAG}/qc/{{ref_name}}.frag_length_overlay.pdf",
            ref_name=frag_ref_names,
        ),
        expand(
            f"{D_FRAG}/qc/{{ref_name}}.ratio_profile.pdf",
            ref_name=frag_ref_names,
        ),
        expand(
            f"{D_FRAG}/qc/{{ref_name}}.arm_zscore_heatmap.pdf",
            ref_name=frag_ref_names,
        ),
        # NMF fragment length features
        expand(
            f"{D_FRAG}/features/{{ref_name}}.nmf_{NMF_N_COMPONENTS}.W.csv",
            ref_name=frag_ref_names,
        ),
        # Motif diversity score
        expand(
            f"{D_FRAG}/features/{{ref_name}}.motif_diversity.csv",
            ref_name=frag_ref_names,
        ),
        # F-profiles
        expand(
            f"{D_FRAG}/features/{{ref_name}}.fprofiles.tsv",
            ref_name=frag_ref_names,
        ),

# ------------------------------------------------------------------------------
# Symlink input FASTQs
# ------------------------------------------------------------------------------

rule symlink_input_fastqs:
    message:
        "Create symlinks for raw input FASTQs"
    input:
        r1 = lambda wc: f"{D_INPUTS}/{samples.r1_map[wc.library_id]}",
        r2 = lambda wc: f"{D_INPUTS}/{samples.r2_map[wc.library_id]}",
    output:
        r1 = f"{D_FRAG}/fastqs/{{library_id}}.raw_R1.fastq.gz",
        r2 = f"{D_FRAG}/fastqs/{{library_id}}.raw_R2.fastq.gz",
    params:
        out_dir = f"{D_FRAG}/fastqs",
    shell:
        """
        mkdir -p "{params.out_dir}"
        ln -sfr "{input.r1}" "{output.r1}"
        ln -sfr "{input.r2}" "{output.r2}"
        """

# ------------------------------------------------------------------------------
# Include module
# ------------------------------------------------------------------------------

include: "frag.smk"
