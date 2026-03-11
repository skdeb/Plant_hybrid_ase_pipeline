# ============================================================
# phASER per-sample allele counting
# ============================================================
# Runs phaser.py then phaser_gene_ae.py for each sample.
#
# MAPQ threshold is data-type specific:
#   RNA-seq        → phaser_mapq_rnaseq (255 = STAR unique mappers)
#   ATAC/ChIP/MOA  → phaser_mapq_other  (matches bowtie2_mapq)
#
# vcf_sample maps the pipeline sample ID to the correct VCF
# column — handles both single-sample and multi-sample VCFs.
#
# Review these outputs before running phaser_pop:
#   {sample}.allelic_counts.txt  per-SNP haplotypic read counts
#   {sample}.gene_ae.txt         per-gene allelic expression
# ============================================================


rule run_phaser:
    """
    Run phaser.py for one sample.
    Produces phased haplotypic read counts at het-SNP positions.
    """
    input:
        bam = get_wasp_bam,
        bai = get_wasp_bai,
        vcf = config["vcf_phaser"],
        tbi = config["vcf_phaser"] + ".tbi"
    output:
        counts     = "results/phaser_sample/{sample}/{sample}.allelic_counts.txt",
        haplotypes = "results/phaser_sample/{sample}/{sample}.haplotypes.txt",
        vcf_out    = "results/phaser_sample/{sample}/{sample}.vcf.gz"
    params:
        outprefix  = "results/phaser_sample/{sample}/{sample}",
        vcf_sample = lambda wc: get_vcf_sample(wc.sample),
        mapq       = lambda wc: get_phaser_mapq(wc),
        baseq      = config["phaser_baseq"],
        paired_end = config["phaser_paired_end"],
        unphased   = config["phaser_unphased_vars"],
        pass_only  = config["phaser_pass_only"],
        extra      = config.get("extra_phaser_params", "")
    conda:
        "../envs/phaser.yaml"
    threads: 4
    log:
        "logs/phaser_sample/{sample}_phaser.log"
    shell:
        """
        mkdir -p results/phaser_sample/{wildcards.sample}

        python $(which phaser.py) \
            --bam {input.bam} \
            --vcf {input.vcf} \
            --sample {params.vcf_sample} \
            --mapq {params.mapq} \
            --baseq {params.baseq} \
            --paired_end {params.paired_end} \
            --unphased_vars {params.unphased} \
            --pass_only {params.pass_only} \
            --threads {threads} \
            --o {params.outprefix} \
            {params.extra} \
        2> {log}
        """


rule run_phaser_gene_ae:
    """
    Summarize phASER haplotypic counts into per-gene allelic expression.
    Gene annotation BED is generated from the GTF on first run using
    phASER's own annotation script. A lock file prevents race conditions
    when multiple samples run in parallel.
    """
    input:
        counts = "results/phaser_sample/{sample}/{sample}.allelic_counts.txt",
        gtf    = config["genome_gtf"]
    output:
        gene_ae = "results/phaser_sample/{sample}/{sample}.gene_ae.txt"
    params:
        outprefix = "results/phaser_sample/{sample}/{sample}",
        bed       = "resources/gene_annotation.bed"
    conda:
        "../envs/phaser.yaml"
    log:
        "logs/phaser_sample/{sample}_gene_ae.log"
    shell:
        """
        if [ ! -f {params.bed} ]; then
            flock {params.bed}.lock \
            python $(which phaser_annotate.py) \
                --gtf {input.gtf} \
                --o {params.bed} \
            2>> {log}
        fi

        python $(which phaser_gene_ae.py) \
            --haplotypic_counts {input.counts} \
            --features {params.bed} \
            --o {params.outprefix} \
        2> {log}
        """
