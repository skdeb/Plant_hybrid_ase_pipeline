# ============================================================
# QC plots
# ============================================================
# Runs an R script that reads all per-sample phASER outputs
# and generates QC visualizations:
#   - biallelic rate per sample
#   - ref/alt ratio distribution
#   - coverage at het-SNPs
#   - ASE effect size distribution
#
# Thresholds from config are passed as arguments and used to
# flag samples on plots. They do NOT hard-filter output data.
# Review plots, then update phaser_pop_samples in config.yaml.
# ============================================================


rule run_qc_plots:
    """
    Generate QC plots from all per-sample phASER outputs.
    All samples in samples.tsv are included regardless of
    phaser_pop_samples — QC review precedes population filtering.
    Touches qc_plots.done on completion.
    """
    input:
        gene_ae = expand(
            "results/phaser_sample/{sample}/{sample}.gene_ae.txt",
            sample=get_all_samples()
        ),
        counts = expand(
            "results/phaser_sample/{sample}/{sample}.allelic_counts.txt",
            sample=get_all_samples()
        )
    output:
        done = "results/qc_plots/qc_plots.done"
    params:
        indir              = "results/phaser_sample",
        outdir             = "results/qc_plots",
        samples_tsv        = config["samples"],
        min_biallelic_rate = config["min_biallelic_rate"],
        max_ref_ratio      = config["max_ref_ratio"],
        min_ref_ratio      = config["min_ref_ratio"],
        effect_size_cutoff = config["effect_size_cutoff"],
        min_coverage       = config["min_coverage"],
        empirical_null     = config["empirical_null"]
    conda:
        "../envs/r_qc.yaml"
    log:
        "logs/qc_plots/qc_plots.log"
    shell:
        """
        mkdir -p {params.outdir}

        Rscript workflow/scripts/qc_plots.R \
            --indir              {params.indir} \
            --outdir             {params.outdir} \
            --samples            {params.samples_tsv} \
            --min_biallelic_rate {params.min_biallelic_rate} \
            --max_ref_ratio      {params.max_ref_ratio} \
            --min_ref_ratio      {params.min_ref_ratio} \
            --effect_size_cutoff {params.effect_size_cutoff} \
            --min_coverage       {params.min_coverage} \
            --empirical_null     {params.empirical_null} \
        2> {log}

        touch {output.done}
        """
