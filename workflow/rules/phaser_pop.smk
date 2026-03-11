# ============================================================
# phASER population-level merge
# ============================================================
# Merges per-sample phASER outputs into a population matrix.
# Only runs if phaser_pop_samples is non-empty in config.yaml.
#
# Fill phaser_pop_samples ONLY after reviewing per-sample QC
# outputs. Include only samples that passed:
#   - biallelic rate >= min_biallelic_rate
#   - ref/alt ratio within expected range
#   - sufficient coverage at het-SNPs
# ============================================================


def get_pop_count_files():
    pop = config.get("phaser_pop_samples", [])
    return expand(
        "results/phaser_sample/{sample}/{sample}.allelic_counts.txt",
        sample=pop
    )

def get_pop_gene_ae_files():
    pop = config.get("phaser_pop_samples", [])
    return expand(
        "results/phaser_sample/{sample}/{sample}.gene_ae.txt",
        sample=pop
    )

def phaser_pop_ids():
    return ",".join(config.get("phaser_pop_samples", []))

def phaser_pop_count_paths():
    pop = config.get("phaser_pop_samples", [])
    return ",".join(
        f"results/phaser_sample/{s}/{s}.allelic_counts.txt"
        for s in pop
    )


rule run_phaser_pop:
    """
    Merge per-sample phASER haplotypic counts into a population matrix.
    Produces:
      phaser_pop.txt          population-level haplotypic counts
      phaser_pop_gene_ae.txt  population-level gene allelic expression
    """
    input:
        counts  = get_pop_count_files(),
        gene_ae = get_pop_gene_ae_files(),
        bed     = "resources/gene_annotation.bed"
    output:
        pop     = "results/phaser_pop/phaser_pop.txt",
        gene_ae = "results/phaser_pop/phaser_pop_gene_ae.txt"
    params:
        outprefix = "results/phaser_pop/phaser_pop",
        ids       = lambda wc: phaser_pop_ids(),
        counts    = lambda wc: phaser_pop_count_paths()
    conda:
        "../envs/phaser.yaml"
    log:
        "logs/phaser_pop/phaser_pop.log"
    shell:
        """
        mkdir -p results/phaser_pop

        python $(which phaser_pop.py) \
            --haplotypic_counts {params.counts} \
            --features {input.bed} \
            --id {params.ids} \
            --o {params.outprefix} \
        2> {log}

        mv {params.outprefix}.haplotypic_counts.txt {output.pop}
        mv {params.outprefix}.gene_ae.txt           {output.gene_ae}
        """
