# ============================================================
# AUTO-GENERATED — DO NOT EDIT DIRECTLY
# Edits will be overwritten on next org-babel tangle.
# 
# Source:  /home/jeszyman/repos/frag/frag.org
# Author:  Jeffrey Szymanski
# Tangled: 2026-03-12 15:04:13
# ============================================================

#########1#########2#########3#########4#########5#########6#########7#########8
#
# This is a modular snakefile, intended to be incorporated into a larger
# workflow using the "include:" directive. (See
# https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html)
#
#########1#########2#########3#########4#########5#########6#########7#########8
rule frag_fastp:
    message:
        "Fragmentomics fastp FASTQ processing"
    conda:
        CONDA_FRAG
    input:
        r1 = f"{D_FRAG}/fastqs/{{library_id}}.raw_R1.fastq.gz",
        r2 = f"{D_FRAG}/fastqs/{{library_id}}.raw_R2.fastq.gz",
    log:
        cmd = f"{D_LOGS}/{{library_id}}_frag_fastp.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_frag_fastp.tsv"
    params:
        extra = config.get("fastp", {}).get("extra", ""),
    threads:
        8
    output:
        failed = f"{D_FRAG}/fastqs/{{library_id}}.failed.fastq.gz",
        html = f"{D_FRAG}/qc/{{library_id}}_frag_fastp.html",
        json = f"{D_FRAG}/qc/{{library_id}}_frag_fastp.json",
        r1 = f"{D_FRAG}/fastqs/{{library_id}}.processed_R1.fastq.gz",
        r2 = f"{D_FRAG}/fastqs/{{library_id}}.processed_R2.fastq.gz",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[fastp] $(date) lib={wildcards.library_id} threads={threads}"

        fastp \
          --detect_adapter_for_pe \
          --disable_quality_filtering \
          --in1 "{input.r1}" --in2 "{input.r2}" \
          --out1 "{output.r1}" --out2 "{output.r2}" \
          --failed_out "{output.failed}" \
          --json "{output.json}" --html "{output.html}" \
          --thread {threads} \
          {params.extra}
        """
rule frag_bwa_index:
    message:
        "Index reference FASTA for BWA alignment"
    conda:
        CONDA_FRAG
    input:
        lambda wc: f"{D_INPUTS}/{config['frag_ref_assemblies'][wc.ref_name]['input']}"
    log:
        cmd = f"{D_LOGS}/{{ref_name}}_bwa_index.log",
    benchmark:
        f"{D_BENCHMARK}/{{ref_name}}_bwa_index.tsv"
    params:
        out_dir = lambda wc: f"{D_FRAG}/ref/bwa/{wc.ref_name}",
        fasta_target = lambda wc: f"{D_FRAG}/ref/bwa/{wc.ref_name}/{wc.ref_name}.fa",
    threads:
        50
    output:
        fa  = f"{D_FRAG}/ref/bwa/{{ref_name}}/{{ref_name}}.fa",
        fai = f"{D_FRAG}/ref/bwa/{{ref_name}}/{{ref_name}}.fa.fai",
        amb = f"{D_FRAG}/ref/bwa/{{ref_name}}/{{ref_name}}.fa.amb",
        ann = f"{D_FRAG}/ref/bwa/{{ref_name}}/{{ref_name}}.fa.ann",
        bwt = f"{D_FRAG}/ref/bwa/{{ref_name}}/{{ref_name}}.fa.bwt",
        pac = f"{D_FRAG}/ref/bwa/{{ref_name}}/{{ref_name}}.fa.pac",
        sa  = f"{D_FRAG}/ref/bwa/{{ref_name}}/{{ref_name}}.fa.sa",
    shell:
        """
        set -euo pipefail
        exec &>> "{log.cmd}"
        echo "[bwa index] $(date) ref_name={wildcards.ref_name} threads={threads}"

        mkdir -p "{params.out_dir}"

        if file -b "{input}" | grep -qi gzip; then
            zcat "{input}" > "{params.fasta_target}"
        else
            cat "{input}" > "{params.fasta_target}"
        fi

        samtools faidx "{params.fasta_target}"

        bwa index "{params.fasta_target}"
        """
rule frag_align:
    message:
        "Fragmentomics alignment with BWA MEM and streaming markdup"
    conda:
        CONDA_FRAG
    input:
        r1  = f"{D_FRAG}/fastqs/{{library_id}}.processed_R1.fastq.gz",
        r2  = f"{D_FRAG}/fastqs/{{library_id}}.processed_R2.fastq.gz",
        ref = f"{D_FRAG}/ref/bwa/{{ref_name}}/{{ref_name}}.fa",
    log:
        cmd = f"{D_LOGS}/{{library_id}}_{{ref_name}}_frag_align.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_{{ref_name}}_frag_align.benchmark.txt"
    threads:
        25
    output:
        bam = f"{D_FRAG}/bams/{{library_id}}.bwa.{{ref_name}}.coorsort.bam",
        bai = f"{D_FRAG}/bams/{{library_id}}.bwa.{{ref_name}}.coorsort.bam.bai",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[bwa mem] $(date) lib={wildcards.library_id} ref={wildcards.ref_name} threads={threads}"

        bash scripts/bwa_mem_markdup_stream.sh \
          "{input.ref}" "{input.r1}" "{input.r2}" \
          "{output.bam}" {threads}
        """
rule frag_filter_alignments:
    message:
        "Filter alignments by MAPQ and genomic region"
    conda:
        CONDA_FRAG
    input:
        bam      = f"{D_FRAG}/bams/{{library_id}}.bwa.{{ref_name}}.coorsort.bam",
        keep_bed = f"{D_FRAG}/ref/keep_5mb.bed",
    log:
        cmd = f"{D_LOGS}/{{library_id}}_{{ref_name}}_frag_filter_alignments.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_{{ref_name}}_frag_filter_alignments.benchmark.txt"
    threads:
        4
    output:
        bam = f"{D_FRAG}/bams/{{library_id}}.bwa.{{ref_name}}.filt.bam",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[filter] $(date) lib={wildcards.library_id} ref={wildcards.ref_name} threads={threads}"

        bash scripts/filter_alignments.sh \
          "{input.bam}" "{input.keep_bed}" {threads} "{output.bam}"
        """
rule frag_bam_to_frag_bed:
    message:
        "Convert filtered BAM to fragment BED with GC and length annotations"
    conda:
        CONDA_FRAG
    input:
        bam   = f"{D_FRAG}/bams/{{library_id}}.bwa.{{ref_name}}.filt.bam",
        fasta = f"{D_FRAG}/ref/bwa/{{ref_name}}/{{ref_name}}.fa",
    log:
        cmd = f"{D_LOGS}/{{library_id}}_{{ref_name}}_frag_bam_to_frag_bed.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_{{ref_name}}_frag_bam_to_frag_bed.benchmark.txt"
    threads:
        1
    output:
        bed = f"{D_FRAG}/frags/{{library_id}}.{{ref_name}}.frag.bed",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[bam2bed] $(date) lib={wildcards.library_id} ref={wildcards.ref_name}"

        bash scripts/bam_to_frag_bed.sh \
          "{input.bam}" "{input.fasta}" "{output.bed}"
        """
rule frag_gc_map_bins:
    message:
        "Create GC and mappability restricted 5Mb bins"
    conda:
        CONDA_FRAG
    input:
        regions = config["gc5mb"],
        blklist = config["blklist"],
        fasta   = expand(f"{D_FRAG}/ref/bwa/{{ref_name}}/{{ref_name}}.fa", ref_name=frag_ref_names)[0],
    log:
        cmd = f"{D_LOGS}/frag_gc_map_bins.log",
    benchmark:
        f"{D_BENCHMARK}/frag_gc_map_bins.tsv"
    threads:
        1
    output:
        keep = f"{D_FRAG}/ref/keep_5mb.bed",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[gc_map_bins] $(date)"

        bash scripts/make_gc_map_bins.sh \
          "{input.regions}" "{input.fasta}" "{input.blklist}" "{output.keep}"
        """
rule frag_gc_distro:
    message:
        "Compute per-library GC distribution from fragment BED"
    conda:
        CONDA_FRAG
    input:
        bed = f"{D_FRAG}/frags/{{library_id}}.{{ref_name}}.frag.bed",
    log:
        cmd = f"{D_LOGS}/{{library_id}}_{{ref_name}}_frag_gc_distro.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_{{ref_name}}_frag_gc_distro.tsv"
    threads:
        1
    output:
        csv = f"{D_FRAG}/frags/{{library_id}}.{{ref_name}}.gc_distro.csv",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[gc_distro] $(date) lib={wildcards.library_id} ref={wildcards.ref_name}"

        Rscript scripts/gc_distro.R \
          "{input.bed}" "{output.csv}"
        """
rule frag_healthy_gc:
    message:
        "Compute median GC distribution from healthy libraries"
    conda:
        CONDA_FRAG
    input:
        csvs = lambda wc: expand(
            f"{D_FRAG}/frags/{{library_id}}.{{ref_name}}.gc_distro.csv",
            library_id=FRAG_HEALTHY_LIBRARIES,
            ref_name=wc.ref_name,
        ),
    log:
        cmd = f"{D_LOGS}/{{ref_name}}_frag_healthy_gc.log",
    benchmark:
        f"{D_BENCHMARK}/{{ref_name}}_frag_healthy_gc.tsv"
    threads:
        1
    output:
        rds = f"{D_FRAG}/frags/{{ref_name}}.healthy_med.rds",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[healthy_gc] $(date) ref={wildcards.ref_name}"

        Rscript scripts/make_healthy_gc_summary.R \
          "{input.csvs}" "{output.rds}"
        """
rule frag_gc_sample:
    message:
        "Resample fragments by healthy GC proportions"
    conda:
        CONDA_FRAG
    input:
        frag_bed    = f"{D_FRAG}/frags/{{library_id}}.{{ref_name}}.frag.bed",
        healthy_med = f"{D_FRAG}/frags/{{ref_name}}.healthy_med.rds",
    log:
        cmd = f"{D_LOGS}/{{library_id}}_{{ref_name}}_frag_gc_sample.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_{{ref_name}}_frag_gc_sample.tsv"
    threads:
        1
    output:
        bed = f"{D_FRAG}/frags/{{library_id}}.{{ref_name}}.sampled_frag.bed",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[gc_sample] $(date) lib={wildcards.library_id} ref={wildcards.ref_name}"

        Rscript scripts/sample_frags_by_gc.R \
          "{input.healthy_med}" "{input.frag_bed}" "{output.bed}"
        """
rule frag_window_sum:
    message:
        "Partition fragments into short (100-150bp) and long (151-220bp) groups"
    conda:
        CONDA_FRAG
    input:
        bed = f"{D_FRAG}/frags/{{library_id}}.{{ref_name}}.sampled_frag.bed",
    log:
        cmd = f"{D_LOGS}/{{library_id}}_{{ref_name}}_frag_window_sum.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_{{ref_name}}_frag_window_sum.tsv"
    threads:
        1
    output:
        short = f"{D_FRAG}/frags/{{library_id}}.{{ref_name}}.norm_short.bed",
        long  = f"{D_FRAG}/frags/{{library_id}}.{{ref_name}}.norm_long.bed",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[window_sum] $(date) lib={wildcards.library_id} ref={wildcards.ref_name}"

        bash scripts/frag_window_sum.sh \
          "{input.bed}" "{output.short}" "{output.long}"
        """
rule frag_window_count:
    message:
        "Count short and long fragments per 5Mb genomic bin"
    conda:
        CONDA_FRAG
    input:
        short  = f"{D_FRAG}/frags/{{library_id}}.{{ref_name}}.norm_short.bed",
        long   = f"{D_FRAG}/frags/{{library_id}}.{{ref_name}}.norm_long.bed",
        matbed = f"{D_FRAG}/ref/keep_5mb.bed",
    log:
        cmd = f"{D_LOGS}/{{library_id}}_{{ref_name}}_frag_window_count.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_{{ref_name}}_frag_window_count.tsv"
    threads:
        1
    output:
        short = f"{D_FRAG}/counts/{{library_id}}.{{ref_name}}.cnt_short.tmp",
        long  = f"{D_FRAG}/counts/{{library_id}}.{{ref_name}}.cnt_long.tmp",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[window_count] $(date) lib={wildcards.library_id} ref={wildcards.ref_name}"

        bash scripts/frag_window_int.sh \
          "{input.short}" "{input.matbed}" "{output.short}"
        bash scripts/frag_window_int.sh \
          "{input.long}" "{input.matbed}" "{output.long}"
        """
rule frag_count_merge:
    message:
        "Merge short and long fragment counts across libraries"
    conda:
        CONDA_FRAG
    input:
        counts = lambda wc: expand(
            f"{D_FRAG}/counts/{{library_id}}.{{ref_name}}.cnt_{{length}}.tmp",
            library_id=FRAG_LIBRARY_IDS,
            ref_name=wc.ref_name,
            length=["short", "long"],
        ),
    log:
        cmd = f"{D_LOGS}/{{ref_name}}_frag_count_merge.log",
    benchmark:
        f"{D_BENCHMARK}/{{ref_name}}_frag_count_merge.tsv"
    threads:
        1
    output:
        tsv = f"{D_FRAG}/frags/{{ref_name}}.frag_counts.tsv",
    params:
        counts_dir = f"{D_FRAG}/counts",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[count_merge] $(date) ref={wildcards.ref_name}"

        bash scripts/count_merge.sh \
          "{params.counts_dir}" "{output.tsv}"
        """
rule frag_ratio_normalize:
    message:
        "Compute zero-centered fragment length ratios per library"
    conda:
        CONDA_FRAG
    input:
        tsv = f"{D_FRAG}/frags/{{ref_name}}.frag_counts.tsv",
    log:
        cmd = f"{D_LOGS}/{{ref_name}}_frag_ratio_normalize.log",
    benchmark:
        f"{D_BENCHMARK}/{{ref_name}}_frag_ratio_normalize.tsv"
    threads:
        1
    output:
        tsv = f"{D_FRAG}/frags/{{ref_name}}.ratios.tsv",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[ratios] $(date) ref={wildcards.ref_name}"

        Rscript scripts/make_ratios.R \
          "{input.tsv}" "{output.tsv}"
        """
rule frag_sample_motifs:
    message:
        "Sample 5-prime end motifs from filtered BAM"
    conda:
        CONDA_FRAG
    input:
        bam   = f"{D_FRAG}/bams/{{library_id}}.bwa.{{ref_name}}.filt.bam",
        fasta = f"{D_FRAG}/ref/bwa/{{ref_name}}/{{ref_name}}.fa",
    log:
        cmd = f"{D_LOGS}/{{library_id}}_{{ref_name}}_frag_sample_motifs.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_{{ref_name}}_frag_sample_motifs.tsv"
    params:
        n_motif = config.get("end_motif", {}).get("n_motif", 4),
        n_reads = config.get("end_motif", {}).get("n_reads", 1000000),
        seed    = config.get("end_motif", {}).get("seed", 42),
    threads:
        4
    output:
        txt = f"{D_FRAG}/motifs/{{library_id}}.{{ref_name}}.motifs.txt",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[motifs] $(date) lib={wildcards.library_id} ref={wildcards.ref_name} threads={threads}"

        bash scripts/sample_motifs.sh \
          "{input.bam}" "{input.fasta}" \
          {params.n_motif} {params.n_reads} {params.seed} {threads} \
          "{output.txt}"
        """
rule frag_motif_matrix:
    message:
        "Build motif frequency matrix across libraries"
    conda:
        CONDA_FRAG
    input:
        txts = lambda wc: expand(
            f"{D_FRAG}/motifs/{{library_id}}.{{ref_name}}.motifs.txt",
            library_id=FRAG_LIBRARY_IDS,
            ref_name=wc.ref_name,
        ),
    log:
        cmd = f"{D_LOGS}/{{ref_name}}_frag_motif_matrix.log",
    benchmark:
        f"{D_BENCHMARK}/{{ref_name}}_frag_motif_matrix.tsv"
    threads:
        1
    output:
        tsv = f"{D_FRAG}/motifs/{{ref_name}}.all_motifs.tsv",
    shell:
        """
        exec &>> "{log.cmd}"
        echo "[motif_matrix] $(date) ref={wildcards.ref_name}"

        Rscript scripts/end_motif_mat.R \
          "{input.txts}" "{output.tsv}"
        """
