# ============================================================
# WASP mapping bias correction
# ============================================================
# Routing is driven by samples.tsv data_type column via the
# output file pattern — Snakemake picks the correct rule
# automatically for each sample:
#
#   RNA-seq → wasp_filter_rnaseq
#     STAR+WASP inline already ran during alignment.
#     This rule filters the vW SAM tag:
#       no vW tag  → read did not overlap any SNP → keep
#       vW:i:1     → overlapped SNP, passed WASP remap → keep
#       vW:i:2+    → failed WASP remap (reference biased) → discard
#     Output: results/wasp/{sample}/{sample}.wasp.bam
#
#   ATAC/ChIP/MOA → 5-step WASP standalone pipeline
#     1. find_intersecting_snps.py  identify reads at het-SNPs
#     2. bowtie2 remap              remap allele-swapped reads
#     3. filter_remapped_reads.py   discard reads that shifted position
#     4. samtools merge             merge kept + filtered reads
#     5. rmdup_pe.py                unbiased random deduplication
#     Output: results/wasp/{sample}/{sample}.wasp.rmdup.bam
#
# IMPORTANT: rmdup_pe.py (WASP) must be used for step 5.
# samtools rmdup preferentially keeps the highest-quality read
# which correlates with the reference allele, reintroducing
# the mapping bias that WASP was designed to eliminate.
# WASP rmdup_pe.py selects a random read among duplicates.
# ============================================================


# ── RNA-seq: vW tag filter ───────────────────────────────────

rule wasp_filter_rnaseq:
    """
    Filter STAR+WASP BAM to keep only unbiased reads.
    """
    input:
        bam = "results/align/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
        bai = "results/align/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
    output:
        bam = "results/wasp/{sample}/{sample}.wasp.bam",
        bai = "results/wasp/{sample}/{sample}.wasp.bam.bai"
    conda:
        "../envs/star.yaml"
    threads: 4
    log:
        "logs/wasp/{sample}_filter_rnaseq.log"
    shell:
        """
        mkdir -p results/wasp/{wildcards.sample}

        samtools view -h -@ {threads} {input.bam} \
        | awk 'BEGIN{{OFS="\t"}}
               /^@/ {{print; next}}
               {{
                 vw = ""
                 for (i=12; i<=NF; i++) {{
                   if ($i ~ /^vW:i:/) {{ vw = substr($i, 6); break }}
                 }}
                 if (vw == "" || vw == "1") print
               }}' \
        | samtools sort -@ {threads} -o {output.bam} \
        2> {log}

        samtools index {output.bam} 2>> {log}
        """


# ── ATAC/ChIP/MOA: WASP standalone 5-step pipeline ──────────

rule wasp_find_intersecting_snps:
    """
    Step 1: Identify reads overlapping het-SNPs.
    Intermediate files are marked temp() and deleted after use.
    """
    input:
        bam  = "results/align/{sample}/{sample}.sorted.bam",
        bai  = "results/align/{sample}/{sample}.sorted.bam.bai",
        done = config["snp_dir"] + "/.done"
    output:
        keep_bam  = temp("results/wasp/{sample}/map1.keep.bam"),
        remap_bam = temp("results/wasp/{sample}/map1.to_remap.bam"),
        remap_fq1 = temp("results/wasp/{sample}/map1.remap.fq1.gz"),
        remap_fq2 = temp("results/wasp/{sample}/map1.remap.fq2.gz")
    params:
        snp_dir  = config["snp_dir"],
        wasp_dir = config["wasp_dir"],
        outdir   = "results/wasp/{sample}",
        bam_base = "{sample}.sorted"
    conda:
        "../envs/wasp.yaml"
    log:
        "logs/wasp/{sample}_find_intersecting.log"
    shell:
        """
        mkdir -p {params.outdir}

        python {params.wasp_dir}/mapping/find_intersecting_snps.py \
            --is_paired_end \
            --is_sorted \
            --output_dir {params.outdir} \
            --snp_dir {params.snp_dir} \
            {input.bam} \
        2> {log}

        mv {params.outdir}/{params.bam_base}.keep.bam     {output.keep_bam}
        mv {params.outdir}/{params.bam_base}.to.remap.bam {output.remap_bam}
        mv {params.outdir}/{params.bam_base}.remap.fq1.gz {output.remap_fq1}
        mv {params.outdir}/{params.bam_base}.remap.fq2.gz {output.remap_fq2}
        """


rule wasp_remap:
    """
    Step 2: Remap allele-swapped reads with bowtie2.
    No MAPQ filter here — filter_remapped_reads.py needs all
    alignments to compare positions correctly.
    """
    input:
        fq1   = "results/wasp/{sample}/map1.remap.fq1.gz",
        fq2   = "results/wasp/{sample}/map1.remap.fq2.gz",
        index = config["bowtie2_index"] + ".1.bt2"
    output:
        bam = temp("results/wasp/{sample}/map2.sorted.bam"),
        bai = temp("results/wasp/{sample}/map2.sorted.bam.bai")
    params:
        index_prefix = config["bowtie2_index"],
        extra        = config.get("extra_bowtie2_params", "")
    conda:
        "../envs/bowtie2.yaml"
    threads:
        config["bowtie2_threads"]
    log:
        "logs/wasp/{sample}_remap.log"
    shell:
        """
        bowtie2 \
            --threads {threads} \
            -x {params.index_prefix} \
            -1 {input.fq1} \
            -2 {input.fq2} \
            --no-mixed \
            --no-discordant \
            {params.extra} \
            2> {log} \
        | samtools sort -@ {threads} -o {output.bam}

        samtools index {output.bam} 2>> {log}
        """


rule wasp_filter_remapped:
    """
    Step 3: Discard reads that map to a different position
    after allele-swapping — these are reference-biased.
    """
    input:
        original = "results/wasp/{sample}/map1.to_remap.bam",
        remapped = "results/wasp/{sample}/map2.sorted.bam"
    output:
        bam = temp("results/wasp/{sample}/map2.keep.bam")
    params:
        wasp_dir = config["wasp_dir"]
    conda:
        "../envs/wasp.yaml"
    log:
        "logs/wasp/{sample}_filter_remapped.log"
    shell:
        """
        python {params.wasp_dir}/mapping/filter_remapped_reads.py \
            {input.original} \
            {input.remapped} \
            {output.bam} \
        2> {log}
        """


rule wasp_merge:
    """
    Step 4: Merge pass-through reads (no SNP overlap) with
    filtered remapped reads into one coordinate-sorted BAM.
    """
    input:
        keep_bam     = "results/wasp/{sample}/map1.keep.bam",
        remapped_bam = "results/wasp/{sample}/map2.keep.bam"
    output:
        bam = temp("results/wasp/{sample}/merged.bam"),
        bai = temp("results/wasp/{sample}/merged.bam.bai")
    conda:
        "../envs/wasp.yaml"
    threads: 4
    log:
        "logs/wasp/{sample}_merge.log"
    shell:
        """
        samtools merge -f -@ {threads} - \
            {input.keep_bam} \
            {input.remapped_bam} \
        | samtools sort -@ {threads} -o {output.bam} \
        2> {log}

        samtools index {output.bam} 2>> {log}
        """


rule wasp_rmdup:
    """
    Step 5: Unbiased random deduplication with WASP rmdup_pe.py.
    This is the final WASP-corrected BAM used by phASER.
    """
    input:
        bam = "results/wasp/{sample}/merged.bam",
        bai = "results/wasp/{sample}/merged.bam.bai"
    output:
        bam = "results/wasp/{sample}/{sample}.wasp.rmdup.bam",
        bai = "results/wasp/{sample}/{sample}.wasp.rmdup.bam.bai"
    params:
        wasp_dir = config["wasp_dir"],
        tmp_bam  = "results/wasp/{sample}/rmdup.unsorted.bam"
    conda:
        "../envs/wasp.yaml"
    threads: 4
    log:
        "logs/wasp/{sample}_rmdup.log"
    shell:
        """
        python {params.wasp_dir}/mapping/rmdup_pe.py \
            {input.bam} \
            {params.tmp_bam} \
        2> {log}

        samtools sort \
            -@ {threads} \
            -o {output.bam} \
            {params.tmp_bam} \
        2>> {log}

        samtools index {output.bam} 2>> {log}

        rm -f {params.tmp_bam}
        """
