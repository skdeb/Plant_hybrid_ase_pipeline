# ============================================================
# VCF preparation
# ============================================================
# Prepares three outputs from vcf_raw:
#
#   vcf_star    → uncompressed VCF, biallelic het SNPs only.
#                 Used by STAR --varVCFfile for WASP inline.
#                 STAR reads REF/ALT alleles only — phased (0|1)
#                 and unphased (0/1) genotypes are treated identically.
#                 STAR requires uncompressed input (not bgzipped).
#
#   vcf_phaser  → bgzipped + tabix indexed, biallelic het SNPs only.
#                 Used by phASER --vcf.
#                 phASER uses phase information (0|1) when present.
#                 If input is unphased (0/1), phASER treats each site
#                 independently — still valid for allele counting.
#
#   snp_dir     → per-chromosome SNP text files for WASP standalone.
#                 Used by ATAC-seq, ChIP-seq, MOA-seq samples only.
#                 Format per file (e.g. chr1.snps.txt):
#                   position  ref_allele  alt_allele
#
# Works with both single-sample and multi-sample VCFs.
# ============================================================


rule vcf_prep_star:
    """
    Prepare uncompressed VCF for STAR --varVCFfile.
    Keeps only biallelic het SNPs. No genotype format conversion needed.
    """
    input:
        vcf = config["vcf_raw"]
    output:
        vcf = config["vcf_star"]
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/vcf_prep/vcf_prep_star.log"
    shell:
        """
        bcftools view \
            --type snps \
            --max-alleles 2 \
            --min-alleles 2 \
            --genotype het \
            {input.vcf} \
        | bcftools annotate \
            --set-id '%CHROM:%POS:%REF:%ALT' \
        > {output.vcf} \
        2> {log}
        """


rule vcf_prep_phaser:
    """
    Prepare bgzipped + tabix indexed VCF for phASER.
    Keeps only biallelic het SNPs.
    Phase information (0|1) is preserved if present in input.
    """
    input:
        vcf = config["vcf_raw"]
    output:
        vcf = config["vcf_phaser"],
        tbi = config["vcf_phaser"] + ".tbi"
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/vcf_prep/vcf_prep_phaser.log"
    shell:
        """
        bcftools view \
            --type snps \
            --max-alleles 2 \
            --min-alleles 2 \
            --genotype het \
            {input.vcf} \
        | bcftools annotate \
            --set-id '%CHROM:%POS:%REF:%ALT' \
        | bgzip -c \
        > {output.vcf} \
        2> {log}

        tabix -p vcf {output.vcf} 2>> {log}
        """


rule vcf_prep_snp_dir:
    """
    Prepare per-chromosome SNP directory for WASP standalone.
    Used by ATAC-seq, ChIP-seq, and MOA-seq samples only.
    RNA-seq uses STAR+WASP inline which reads vcf_star directly.
    Touches a .done sentinel file on completion.
    """
    input:
        vcf = config["vcf_raw"]
    output:
        done = config["snp_dir"] + "/.done"
    params:
        snp_dir = config["snp_dir"]
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/vcf_prep/vcf_prep_snp_dir.log"
    shell:
        """
        mkdir -p {params.snp_dir}

        bcftools view \
            --type snps \
            --max-alleles 2 \
            --min-alleles 2 \
            --genotype het \
            {input.vcf} \
        | bcftools query \
            -f '%CHROM\t%POS\t%REF\t%ALT\n' \
        | awk '{{
            chrom = $1
            print $2"\t"$3"\t"$4 \
                >> "{params.snp_dir}/"chrom".snps.txt"
        }}' \
        2> {log}

        touch {output.done}
        """
