# =============================================================================
# Rule module: Report Generation
# Tool: MultiQC  – aggregate QC metrics from all steps into one HTML report
# =============================================================================

rule multiqc_report:
    """Aggregate all QC metrics into a single MultiQC HTML report."""
    input:
        # FastQC (raw)
        expand(
            f"{OUT}/qc/fastqc/{{sample}}_{{read}}_fastqc.zip",
            sample=SAMPLES, read=["R1", "R2"]
        ),
        # FastQC (trimmed)
        expand(
            f"{OUT}/qc/fastqc/{{sample}}_{{read}}_trimmed_fastqc.zip",
            sample=SAMPLES, read=["R1", "R2"]
        ),
        # fastp JSON reports
        expand(
            f"{OUT}/qc/fastp/{{sample}}/{{sample}}_fastp.json",
            sample=SAMPLES
        ),
        # Bowtie2 alignment stats (host removal)
        expand(
            f"{OUT}/host_removed/{{sample}}/{{sample}}_bowtie2.log",
            sample=SAMPLES
        ),
        # Kraken2 reports
        expand(
            f"{OUT}/taxonomy/{{sample}}/{{sample}}.kraken2.report",
            sample=SAMPLES
        ),
    output:
        report = f"{OUT}/report/multiqc_report.html",
        data   = directory(f"{OUT}/report/multiqc_data"),
    params:
        search_dirs = [
            f"{OUT}/qc",
            f"{OUT}/host_removed",
            f"{OUT}/taxonomy",
        ],
        outdir = f"{OUT}/report",
        extra  = config["multiqc"]["extra_args"],
    conda: "../../envs/qc.yaml"
    log: f"{OUT}/logs/multiqc.log"
    shell:
        """
        multiqc \
            {params.extra} \
            --outdir {params.outdir} \
            {params.search_dirs} \
            > {log} 2>&1
        """
