rule frag_fastp:
    message: "Fragmentomics fastp FASTQ processing"
    conda: CONDA_FRAG
    threads: 8
    params:
        extra = config.get("fastp", {}).get("extra", ""),
    input:
        r1 = f"{D_FRAG}/fastqs/{{library_id}}.raw_R1.fastq.gz",
        r2 = f"{D_FRAG}/fastqs/{{library_id}}.raw_R2.fastq.gz",
    output:
        failed = f"{D_FRAG}/fastqs/{{library_id}}.failed.fastq.gz",
        html = f"{D_FRAG}/qc/{{library_id}}_frag_fastp.html",
        json = f"{D_FRAG}/qc/{{library_id}}_frag_fastp.json",
        r1     = f"{D_FRAG}/fastqs/{{library_id}}.processed_R1.fastq.gz",
        r2     = f"{D_FRAG}/fastqs/{{library_id}}.processed_R2.fastq.gz",
    log:
        cmd  = f"{D_LOGS}/{{library_id}}_frag_fastp.log",
    benchmark:
        f"{D_BENCHMARK}/{{library_id}}_frag_fastp.tsv"
    shell:
        """
        # Logging and console output
        exec &>> "{log.cmd}"
        echo "[fastp] $(date) lib={wildcards.library_id} threads={threads}"

        # Main
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
    message: "Index reference fasta for BWA alignment",
    conda: CONDA_FRAG,
    input:
        lambda wc: f"{D_INPUTS}/{config['frag_ref_assemblies'][wc.ref_name]['input']}"
    output:
        fa   = f"{D_FRAG}/ref/bwa/{{ref_name}}/{{ref_name}}.fa",
        fai  = f"{D_FRAG}/ref/bwa/{{ref_name}}/{{ref_name}}.fa.fai",
        amb  = f"{D_FRAG}/ref/bwa/{{ref_name}}/{{ref_name}}.fa.amb",
        ann  = f"{D_FRAG}/ref/bwa/{{ref_name}}/{{ref_name}}.fa.ann",
        bwt  = f"{D_FRAG}/ref/bwa/{{ref_name}}/{{ref_name}}.fa.bwt",
        pac  = f"{D_FRAG}/ref/bwa/{{ref_name}}/{{ref_name}}.fa.pac",
        sa   = f"{D_FRAG}/ref/bwa/{{ref_name}}/{{ref_name}}.fa.sa",
    params:
        out_dir = lambda wc: f"{D_FRAG}/ref/bwa/{wc.ref_name}",
        fasta_target = lambda wc: f"{D_FRAG}/ref/bwa/{wc.ref_name}/{wc.ref_name}.fa",
    log:
        cmd = f"{D_LOGS}/{{ref_name}}_bwa_index.log",
    benchmark:
        f"{D_BENCHMARK}/{{ref_name}}_bwa_index.tsv"
    threads: 50
    shell:
        r"""
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
    message: "Fragmentomics alignment with BWA MEM",
    conda: CONDA_FRAG,
    input:
        r1     = f"{D_FRAG}/fastqs/{{library_id}}.processed_R1.fastq.gz",
        r2     = f"{D_FRAG}/fastqs/{{library_id}}.processed_R2.fastq.gz",
        ref = lambda wc: f"{D_FRAG}/ref/bwa/{{ref_name}}/{{ref_name}}.fa",
    output:
        bam = f"{D_FRAG}/bams/{{library_id}}.bwa.{{ref_name}}.coorsort.bam",
    threads: 25
    shell:
        """
        bwa mem -M -t {threads} \
        {input.ref} {input.r1} {input.r2} \
        | samtools view -@ 4 -Sb - -o - \
        | samtools sort -@ 4 - -o {output.bam}
        samtools index -@ 4 {output.bam}
        """
