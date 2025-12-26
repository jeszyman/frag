# ==============================================================================
# CFDNA FRAGMENTOMICS FULL PIPELINE TEST WRAPPER SNAKEFILE
# ==============================================================================

# ------------------------------------------------------------------------------
#
# Comment
#
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

# Load packages
import os
import pandas as pd

# --------------------------------------------------------------------------------------
# Load YAML Configuration
# --------------------------------------------------------------------------------------

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

D_INPUTS = config["directories"]["inputs"]
D_FRAG = config["directories"]["frag"]
D_LOGS = config["directories"]["logs"]
D_BENCHMARK = config["directories"]["benchmark"]

CONDA_FRAG = config["conda"]["frag"]

frag_ref_names = ["chr22"]

# --------------------------------------------------------------------------------------
# Load Tablular Configuration
# --------------------------------------------------------------------------------------

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
        # Map unmerged_library_id → R1 filename
        return dict(zip(self.df["library_id"], self.df["r1_basename"]))

    @property
    def r2_map(self):
        # Map unmerged_library_id → R1 filename
        return dict(zip(self.df["library_id"], self.df["r2_basename"]))

samples = SampleTable(
    tsv_path=config["sample-tsv-path"],
    selected_ids=["lib001","lib002"]
)

# ------------------------------------------------------------------------------
# Wrapper Rules
# ------------------------------------------------------------------------------

# All rule

rule all:
    input:
        # FASTQs
        expand(f"{D_FRAG}/fastqs/{{library_id}}.{{processing}}_{{read}}.fastq.gz",
               library_id = samples.frag_library_ids,
               processing = ["raw","processed"],
               read = ["R1", "R2"]),
        #
        # Alignments
        #
        # Index
        expand(f"{D_FRAG}/ref/bwa/{{ref_name}}/{{ref_name}}.fa.sa",
               ref_name = frag_ref_names),
        #
        # Align
        expand(f"{D_FRAG}/bams/{{library_id}}.bwa.{{ref_name}}.coorsort.bam",
               library_id = samples.frag_library_ids,
               ref_name = frag_ref_names),

# Symlink input FASTQ files

rule symlink_input_fastqs:
    input:
        r1 = lambda wc: f"{D_INPUTS}/{samples.r1_map[wc.library_id]}",
        r2 = lambda wc: f"{D_INPUTS}/{samples.r2_map[wc.library_id]}",
    output:
        r1=f"{D_FRAG}/fastqs/{{library_id}}.raw_R1.fastq.gz",
        r2=f"{D_FRAG}/fastqs/{{library_id}}.raw_R2.fastq.gz",
    params:
        out_dir=f"{D_FRAG}/fastqs",
    shell:
        r"""
        mkdir -p "{params.out_dir}"
        ln -sfr "{input.r1}" "{output.r1}"
        ln -sfr "{input.r2}" "{output.r2}"
        """

# --------------------------------------------------------------------------------------
# Snakemake Includes
# --------------------------------------------------------------------------------------

include: "frag.smk"
