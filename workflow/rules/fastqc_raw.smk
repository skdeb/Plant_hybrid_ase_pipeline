# ============================================================
# Raw FASTQ QC
# ============================================================


rule fastqc_raw_sample:
    """
    Run FastQC on raw paired-end FASTQs for one sample.
    Outputs are renamed to sample-based names since FastQC
    names outputs after the input filename prefix.
    """
    input:
        r1 = get_raw_r1,
        r2 = get_raw_r2
    output:
        html_r1 = "results/fastqc_raw/{sample}/{sample}_R1_fastqc.html",
        html_r2 = "results/fastqc_raw/{sample}/{sample}_R2_fastqc.html",
        zip_r1  = "results/fastqc_raw/{sample}/{sample}_R1_fastqc.zip",
        zip_r2  = "results/fastqc_raw/{sample}/{sample}_R2_fastqc.zip"
    params:
        outdir = "results/fastqc_raw/{sample}"
    conda:
        "../envs/fastqc.yaml"
    threads: 2
    log:
        "logs/fastqc_raw/{sample}.log"
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


rule multiqc_raw:
    """
    Aggregate raw FastQC reports into a single MultiQC report.
    """
    input:
        expand(
            "results/fastqc_raw/{sample}/{sample}_R1_fastqc.zip",
            sample=get_all_samples()
        )
    output:
        html = "results/fastqc_raw/multiqc/multiqc_report.html"
    params:
        indir  = "results/fastqc_raw",
        outdir = "results/fastqc_raw/multiqc"
    conda:
        "../envs/fastqc.yaml"
    log:
        "logs/fastqc_raw/multiqc.log"
    shell:
        """
        multiqc \
            --force \
            --outdir {params.outdir} \
            {params.indir} \
        2> {log}
        """
