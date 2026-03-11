# ============================================================
# Genome index creation — one-time setup
# ============================================================
# Run these BEFORE the main pipeline if indices do not exist.
#
# Usage:
#   snakemake build_star_index    --cores 8  --use-conda
#   snakemake build_bowtie2_index --cores 8  --use-conda
#   snakemake build_indices       --cores 16 --use-conda
#
# If you already have indices, set the paths in config.yaml:
#   star_index:    path/to/existing/star_index/
#   bowtie2_index: path/to/existing/bowtie2_index/genome
# These rules will never run if the index already exists.
#
# Resource guidance for maize genome (2.4 Gb):
#   STAR index:    ~30-60 min, ~30 Gb RAM
#   bowtie2 index: ~15-30 min, ~8 Gb RAM
# ============================================================


rule build_star_index:
    """
    Build STAR genome index.

    star_sjdb_overhang: default 100, works well for most read lengths
    including 75 bp and 150 bp. Only reduce if reads are shorter than
    50 bp (set to read_length - 1 in that case).

    genomeSAindexNbases: controls suffix array index size. Default 14
    is appropriate for maize and most plant genomes. Reduce for genomes
    smaller than 100 Mb (e.g. --genomeSAindexNbases 11).
    """
    input:
        fasta = config["genome_fasta"],
        gtf   = config["genome_gtf"]
    output:
        done = config["star_index"] + "/.done"
    params:
        index_dir     = config["star_index"],
        sjdb_overhang = config.get("star_sjdb_overhang", 100),
        extra         = config.get("extra_star_index_params", "")
    conda:
        "../envs/star.yaml"
    threads:
        config["star_threads"]
    log:
        "logs/index/star_index.log"
    shell:
        """
        mkdir -p {params.index_dir}

        STAR \
            --runMode genomeGenerate \
            --runThreadN {threads} \
            --genomeDir {params.index_dir} \
            --genomeFastaFiles {input.fasta} \
            --sjdbGTFfile {input.gtf} \
            --sjdbOverhang {params.sjdb_overhang} \
            {params.extra} \
        2> {log}

        touch {output.done}
        """


rule build_bowtie2_index:
    """
    Build bowtie2 genome index.
    Output prefix is set by bowtie2_index in config.
    The index files produced will be:
      {bowtie2_index}.1.bt2  {bowtie2_index}.2.bt2
      {bowtie2_index}.3.bt2  {bowtie2_index}.4.bt2
      {bowtie2_index}.rev.1.bt2  {bowtie2_index}.rev.2.bt2
    bowtie2_align uses the .1.bt2 file as an existence check input.
    """
    input:
        fasta = config["genome_fasta"]
    output:
        bt2_1  = config["bowtie2_index"] + ".1.bt2",
        bt2_2  = config["bowtie2_index"] + ".2.bt2",
        bt2_3  = config["bowtie2_index"] + ".3.bt2",
        bt2_4  = config["bowtie2_index"] + ".4.bt2",
        bt2_r1 = config["bowtie2_index"] + ".rev.1.bt2",
        bt2_r2 = config["bowtie2_index"] + ".rev.2.bt2"
    params:
        index_prefix = config["bowtie2_index"],
        extra        = config.get("extra_bowtie2_index_params", "")
    conda:
        "../envs/bowtie2.yaml"
    threads:
        config["bowtie2_threads"]
    log:
        "logs/index/bowtie2_index.log"
    shell:
        """
        mkdir -p $(dirname {params.index_prefix})

        bowtie2-build \
            --threads {threads} \
            {params.extra} \
            {input.fasta} \
            {params.index_prefix} \
        2> {log}
        """


rule build_indices:
    """
    Convenience target to build both STAR and bowtie2 indices.
    Only builds what is missing — existing indices are skipped.
    """
    input:
        star    = config["star_index"] + "/.done",
        bowtie2 = config["bowtie2_index"] + ".1.bt2"
