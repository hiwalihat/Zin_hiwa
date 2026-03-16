# =============================================================================
# Rule module: Metagenomic Binning
# Tools: MetaBAT2 (binning) + CheckM (bin quality assessment)
# =============================================================================

rule metabat2_bin:
    """
    Bin assembled contigs using MetaBAT2.
    Requires per-contig coverage information derived from read mappings.
    """
    input:
        contigs = rules.megahit_assembly.output.contigs,
        bam     = rules.map_reads_to_contigs.output.bam,
        bai     = rules.map_reads_to_contigs.output.bai,
    output:
        depth   = f"{OUT}/binning/{{sample}}/depth.txt",
        bin_dir = directory(f"{OUT}/binning/{{sample}}/bins"),
    params:
        prefix         = f"{OUT}/binning/{{sample}}/bins/bin",
        min_contig     = config["metabat2"]["min_contig_size"],
        extra          = config["metabat2"]["extra_args"],
    threads: config["metabat2"]["threads"]
    conda: "../../envs/binning.yaml"
    log: f"{OUT}/logs/binning/{{sample}}_metabat2.log"
    shell:
        """
        # Compute per-contig depth
        jgi_summarize_bam_contig_depths \
            --outputDepth {output.depth} \
            {input.bam} \
            >> {log} 2>&1

        mkdir -p {output.bin_dir}

        # Run MetaBAT2
        metabat2 \
            -i {input.contigs} \
            -a {output.depth} \
            -o {params.prefix} \
            -m {params.min_contig} \
            --numThreads {threads} \
            {params.extra} \
            >> {log} 2>&1
        """


rule checkm_qa:
    """Assess bin completeness and contamination with CheckM."""
    input:
        bin_dir = rules.metabat2_bin.output.bin_dir,
    output:
        summary = f"{OUT}/binning/{{sample}}/checkm/checkm_summary.tsv",
    params:
        checkm_dir   = f"{OUT}/binning/{{sample}}/checkm",
        checkm_data  = config["databases"]["checkm_data"],
        completeness = config["checkm"]["completeness_threshold"],
        contamination= config["checkm"]["contamination_threshold"],
    threads: config["checkm"]["threads"]
    conda: "../../envs/binning.yaml"
    log: f"{OUT}/logs/binning/{{sample}}_checkm.log"
    shell:
        """
        export CHECKM_DATA_PATH={params.checkm_data}

        checkm lineage_wf \
            -t {threads} \
            -x fa \
            {input.bin_dir} \
            {params.checkm_dir} \
            > {log} 2>&1

        checkm qa \
            -o 2 \
            --tab_table \
            -f {output.summary} \
            {params.checkm_dir}/lineage.ms \
            {params.checkm_dir} \
            >> {log} 2>&1
        """


rule filter_bins:
    """
    Filter bins by completeness and contamination thresholds.
    Outputs a directory containing only high-quality bins.
    """
    input:
        summary = rules.checkm_qa.output.summary,
        bin_dir = rules.metabat2_bin.output.bin_dir,
    output:
        hq_dir  = directory(f"{OUT}/binning/{{sample}}/high_quality_bins"),
        hq_list = f"{OUT}/binning/{{sample}}/high_quality_bins.txt",
    params:
        completeness  = config["checkm"]["completeness_threshold"],
        contamination = config["checkm"]["contamination_threshold"],
    conda: "../../envs/binning.yaml"
    log: f"{OUT}/logs/binning/{{sample}}_filter.log"
    script:
        "../scripts/filter_bins.py"
