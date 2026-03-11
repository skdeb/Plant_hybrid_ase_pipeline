# ============================================================
# Read trimming with fastp
# ============================================================
# Parameters (min_length, quality, threads) are data-type
# specific and looked up from config by sample data_type.
# ============================================================


rule fastp_trim:
    """
    Trim paired-end reads with fastp.
    Adapter detection is automatic (--detect_adapter_for_pe).
    Parameters are looked up per data_type from config fastp section.
    """
    input:
        r1 = get_raw_r1,
        r2 = get_raw_r2
    output:
        r1   = "results/trim/{sample}/{sample}_R1.fastq.gz",
        r2   = "results/trim/{sample}/{sample}_R2.fastq.gz",
        html = "results/trim/{sample}/{sample}_fastp.html",
        json = "results/trim/{sample}/{sample}_fastp.json"
    params:
        min_length = lambda wc: config["fastp"][
            samples.loc[wc.sample, "data_type"]]["min_length"],
        quality    = lambda wc: config["fastp"][
            samples.loc[wc.sample, "data_type"]]["quality"],
        extra      = config.get("extra_fastp_params", "")
    conda:
        "../envs/fastp.yaml"
    threads:
        lambda wc: config["fastp"][
            samples.loc[wc.sample, "data_type"]]["threads"]
    log:
        "logs/trim/{sample}.log"
    shell:
        """
        mkdir -p results/trim/{wildcards.sample}

        fastp \
            --in1 {input.r1} \
            --in2 {input.r2} \
            --out1 {output.r1} \
            --out2 {output.r2} \
            --html {output.html} \
            --json {output.json} \
            --length_required {params.min_length} \
            --qualified_quality_phred {params.quality} \
            --thread {threads} \
            --detect_adapter_for_pe \
            {params.extra} \
        2> {log}
        """
