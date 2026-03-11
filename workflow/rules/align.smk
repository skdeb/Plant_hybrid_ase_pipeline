# ============================================================
# Alignment
# ============================================================
# RNA-seq    → STAR with hardcoded flags:
#                --twopassMode Basic      (improved junction detection;
#                                          ~2x runtime, ~30-50 Gb temp/sample)
#                --waspOutputMode SAMtag  (inline WASP bias correction)
#                --outSAMtype BAM SortedByCoordinate
#                --outSAMattributes NH HI AS NM MD vA vW
#              _STARtmp and _STARpass1 deleted after completion.
#
# ATAC/ChIP/MOA → bowtie2 paired-end, MAPQ filtered, sorted + indexed
#
# Index creation is handled separately in index.smk.
# If indices already exist, point star_index and bowtie2_index
# in config.yaml to them — index rules will be skipped.
#
# Routing:
#   Snakemake routes each sample to the correct rule via output
#   file pattern. RNA-seq produces:
#     results/align/{sample}/{sample}.Aligned.sortedByCoord.out.bam
#   ATAC/ChIP/MOA produces:
#     results/align/{sample}/{sample}.sorted.bam
#   get_bam() in the Snakefile returns the correct path per sample.
# ============================================================


rule star_align:
    """
    Align RNA-seq paired-end reads with STAR.
    Two-pass mode + WASP inline bias correction.
    BAM is coordinate-sorted and indexed with samtools.
    Temp directories (_STARtmp, _STARpass1) are removed after completion.

    Index input uses .done sentinel produced by build_star_index.
    If you provide a pre-built index, touch a .done file inside it:
      touch resources/star_index/.done
    """
    input:
        r1    = get_trimmed_r1,
        r2    = get_trimmed_r2,
        index = config["star_index"] + "/.done",
        vcf   = config["vcf_star"]
    output:
        bam = "results/align/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        bai = "results/align/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai",
        log = "results/align/{sample}/{sample}.Log.final.out",
        sj  = "results/align/{sample}/{sample}.SJ.out.tab"
    params:
        prefix           = "results/align/{sample}/{sample}.",
        index_dir        = config["star_index"],
        gtf              = config["genome_gtf"],
        max_multimappers = config["star_max_multimappers"],
        mapq_unique      = config["star_mapq_unique"],
        extra            = config.get("extra_star_params", "")
    conda:
        "../envs/star.yaml"
    threads:
        config["star_threads"]
    log:
        "logs/align/{sample}_star.log"
    shell:
        """
        mkdir -p results/align/{wildcards.sample}

        STAR \
            --runThreadN {threads} \
            --genomeDir {params.index_dir} \
            --sjdbGTFfile {params.gtf} \
            --readFilesIn {input.r1} {input.r2} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.prefix} \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes NH HI AS NM MD vA vW \
            --outFilterMultimapNmax {params.max_multimappers} \
            --outSAMmapqUnique {params.mapq_unique} \
            --twopassMode Basic \
            --waspOutputMode SAMtag \
            --varVCFfile {input.vcf} \
            --outSAMattrRGline ID:{wildcards.sample} SM:{wildcards.sample} \
            {params.extra} \
        2> {log}

        samtools index {output.bam} 2>> {log}

        # Remove large temp directories created by two-pass mode
        rm -rf {params.prefix}_STARtmp
        rm -rf {params.prefix}_STARpass1
        """


rule bowtie2_align:
    """
    Align ATAC-seq, ChIP-seq, or MOA-seq paired-end reads with bowtie2.
    --no-mixed and --no-discordant enforce proper pairs only.
    Reads below bowtie2_mapq are discarded before WASP filtering.
    BAM is coordinate-sorted and indexed with samtools.

    Index input uses .1.bt2 file produced by build_bowtie2_index.
    If you provide a pre-built index, the .1.bt2 file must exist
    at the bowtie2_index prefix path.
    """
    input:
        r1    = get_trimmed_r1,
        r2    = get_trimmed_r2,
        index = config["bowtie2_index"] + ".1.bt2"
    output:
        bam = "results/align/{sample}/{sample}.sorted.bam",
        bai = "results/align/{sample}/{sample}.sorted.bam.bai"
    params:
        index_prefix = config["bowtie2_index"],
        mapq         = config["bowtie2_mapq"],
        extra        = config.get("extra_bowtie2_params", "")
    conda:
        "../envs/bowtie2.yaml"
    threads:
        config["bowtie2_threads"]
    log:
        "logs/align/{sample}_bowtie2.log"
    shell:
        """
        mkdir -p results/align/{wildcards.sample}

        bowtie2 \
            --threads {threads} \
            -x {params.index_prefix} \
            -1 {input.r1} \
            -2 {input.r2} \
            --no-mixed \
            --no-discordant \
            {params.extra} \
            2> {log} \
        | samtools view -bS -q {params.mapq} \
        | samtools sort -@ {threads} -o {output.bam}

        samtools index {output.bam} 2>> {log}
        """
