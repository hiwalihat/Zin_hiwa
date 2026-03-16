# =============================================================================
# Metagenomics Analysis Pipeline
# =============================================================================
# Full shotgun metagenomics workflow:
#   1. Quality control          (FastQC + fastp)
#   2. Host read removal        (Bowtie2)
#   3. De novo assembly         (MEGAHIT)
#   4. Taxonomic classification (Kraken2 + Bracken)
#   5. Functional profiling     (HUMAnN3)
#   6. Metagenomic binning      (MetaBAT2 + CheckM)
#   7. Report generation        (MultiQC)
# =============================================================================

configfile: "config/config.yaml"

# Make config values available under short names
SAMPLES    = config["samples"]
OUT        = config["output_dir"]

# ---------------------------------------------------------------------------
# Include rule modules
# ---------------------------------------------------------------------------
include: "workflow/rules/qc.smk"
include: "workflow/rules/host_removal.smk"
include: "workflow/rules/assembly.smk"
include: "workflow/rules/taxonomy.smk"
include: "workflow/rules/functional.smk"
include: "workflow/rules/binning.smk"
include: "workflow/rules/report.smk"

# ---------------------------------------------------------------------------
# Target rule – request all final outputs
# ---------------------------------------------------------------------------
rule all:
    input:
        # --- QC ---
        expand(
            "{out}/qc/fastqc/{sample}_{read}_fastqc.html",
            out=OUT, sample=SAMPLES, read=["R1", "R2"]
        ),
        expand(
            "{out}/qc/fastp/{sample}/{sample}_R1_trimmed.fastq.gz",
            out=OUT, sample=SAMPLES
        ),
        # --- Host removal ---
        expand(
            "{out}/host_removed/{sample}/{sample}_R1_clean.fastq.gz",
            out=OUT, sample=SAMPLES
        ),
        # --- Assembly ---
        expand(
            "{out}/assembly/{sample}/final.contigs.fa",
            out=OUT, sample=SAMPLES
        ),
        # --- Taxonomy ---
        expand(
            "{out}/taxonomy/{sample}/{sample}.kraken2.report",
            out=OUT, sample=SAMPLES
        ),
        expand(
            "{out}/taxonomy/{sample}/{sample}.bracken",
            out=OUT, sample=SAMPLES
        ),
        # --- Functional profiling ---
        expand(
            "{out}/functional/{sample}/{sample}_pathabundance.tsv",
            out=OUT, sample=SAMPLES
        ),
        # --- Binning ---
        expand(
            "{out}/binning/{sample}/checkm/checkm_summary.tsv",
            out=OUT, sample=SAMPLES
        ),
        # --- Report ---
        f"{OUT}/report/multiqc_report.html",
