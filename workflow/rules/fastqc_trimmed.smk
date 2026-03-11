# ============================================================
# Trimmed FASTQ QC
# ============================================================
# Input is resolved via get_trimmed_r1/r2:
#   start_from: fastq   → uses results/trim/ output
#   start_from: trimmed → uses path from samples.tsv
# ============================================================


rule fastqc_trimmed_sample:
    """
    Run FastQC on trimmed paired-end FASTQs for one sample.
    """
    input:
        r1 = get_trimmed_r1,
        r2 = get_trimmed_r2
    output:
        html_r1 = "results/fastqc_trimmed/{sample}/{sample}_R1_fastqc.html",
        html_r2 = "results/fastqc_trimmed/{sample}/{sample}_R2_fastqc.html",
        zip_r1  = "results/fastqc_trimmed/{sample}/{sample}_R1_fastqc.zip",
        zip_r2  = "results/fastqc_trimmed/{sample}/{sample}_R2_fastqc.zip"
    params:
        outdir = "results/fastqc_trimmed/{sample}"
    conda:
        "../envs/fastqc.yaml"
    threads: 2
    log:
        "logs/fastqc_trimmed/{sample}.log"
    shell:
        """
        mkdir -p {params.outdir}

        fastqc \
            --threads {threads} \
            --outdir {params.outdir} \
            {input.r1} {input.r2} \
        2> {log}

        r1_base=$(basename {input.r1} .fastq.gz)
        r2_base=$(basename {input.r2} .fastq.gz)

        [ -f "{params.outdir}/${{r1_base}}_fastqc.html" ] && \
            mv "{params.outdir}/${{r1_base}}_fastqc.html" {output.html_r1}
        [ -f "{params.outdir}/${{r1_base}}_fastqc.zip" ] && \
            mv "{params.outdir}/${{r1_base}}_fastqc.zip"  {output.zip_r1}
        [ -f "{params.outdir}/${{r2_base}}_fastqc.html" ] && \
            mv "{params.outdir}/${{r2_base}}_fastqc.html" {output.html_r2}
        [ -f "{params.outdir}/${{r2_base}}_fastqc.zip" ] && \
            mv "{params.outdir}/${{r2_base}}_fastqc.zip"  {output.zip_r2}
        """


rule multiqc_trimmed:
    """
    Aggregate trimmed FastQC + fastp JSON reports into MultiQC.
    fastp JSON is included only when start_from: fastq (trim ran).
    """
    input:
        fastqc = expand(
            "results/fastqc_trimmed/{sample}/{sample}_R1_fastqc.zip",
            sample=get_all_samples()
        ),
        fastp = (
            expand(
                "results/trim/{sample}/{sample}_fastp.json",
                sample=get_all_samples()
            )
            if config["start_from"] == "fastq" else []
        )
    output:
        html = "results/fastqc_trimmed/multiqc/multiqc_report.html"
    params:
        indirs = (
            "results/fastqc_trimmed results/trim"
            if config["start_from"] == "fastq"
            else "results/fastqc_trimmed"
        ),
        outdir = "results/fastqc_trimmed/multiqc"
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/fastqc_trimmed/multiqc.log"
    shell:
        """
        multiqc \
            --force \
            --outdir {params.outdir} \
            {params.indirs} \
        2> {log}
        """
