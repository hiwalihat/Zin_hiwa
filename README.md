# Zin_hiwa – Full Metagenomics Analysis Pipeline

A production-ready, reproducible shotgun metagenomics workflow built with
[Snakemake](https://snakemake.readthedocs.io/) and managed with
[conda](https://docs.conda.io/).

---

## Pipeline Overview

```
Raw paired-end FASTQ reads
        │
        ▼
┌───────────────────────────────┐
│  1. Quality Control           │  FastQC (pre & post) + fastp (trimming)
└───────────────┬───────────────┘
                │
                ▼
┌───────────────────────────────┐
│  2. Host Read Removal         │  Bowtie2 → discard human/host reads
└───────────────┬───────────────┘
                │
        ┌───────┴────────┐
        │                │
        ▼                ▼
┌──────────────┐  ┌──────────────────────────────┐
│ 3. Assembly  │  │  4. Taxonomic Classification  │
│  (MEGAHIT)   │  │  Kraken2 → Bracken            │
└──────┬───────┘  └──────────────────────────────┘
       │
       ├──────────────────────────────┐
       │                              │
       ▼                              ▼
┌──────────────────────┐   ┌──────────────────────────┐
│ 5. Functional        │   │  6. Metagenomic Binning   │
│    Profiling         │   │  MetaBAT2 → CheckM QC     │
│    (HUMAnN3)         │   └──────────────────────────┘
└──────────────────────┘
        │
        ▼
┌───────────────────────────────┐
│  7. Aggregated Report         │  MultiQC HTML report
└───────────────────────────────┘
```

### Tools Used

| Step | Tool | Version |
|------|------|---------|
| Quality control | [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | 0.12.1 |
| Adapter trimming | [fastp](https://github.com/OpenGene/fastp) | 0.23.4 |
| Host removal | [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/) | 2.5.3 |
| Read alignment | [SAMtools](http://www.htslib.org/) | 1.19 |
| De novo assembly | [MEGAHIT](https://github.com/voutcn/megahit) | 1.2.9 |
| Taxonomic classification | [Kraken2](https://ccb.jhu.edu/software/kraken2/) | 2.1.3 |
| Abundance re-estimation | [Bracken](https://ccb.jhu.edu/software/bracken/) | 2.9 |
| Functional profiling | [HUMAnN3](https://huttenhower.sph.harvard.edu/humann/) | 3.8 |
| Metagenomic binning | [MetaBAT2](https://bitbucket.org/berkeleylab/metabat) | 2.15 |
| Bin quality assessment | [CheckM](https://ecogenomics.github.io/CheckM/) | 1.2.2 |
| Report aggregation | [MultiQC](https://multiqc.info/) | 1.21 |

---

## Repository Structure

```
Zin_hiwa/
├── Snakefile                     # Main workflow entry point
├── config/
│   └── config.yaml               # All user-configurable parameters
├── workflow/
│   ├── rules/
│   │   ├── qc.smk                # Rule: FastQC + fastp
│   │   ├── host_removal.smk      # Rule: Bowtie2 host removal
│   │   ├── assembly.smk          # Rule: MEGAHIT assembly + indexing
│   │   ├── taxonomy.smk          # Rule: Kraken2 + Bracken
│   │   ├── functional.smk        # Rule: HUMAnN3
│   │   ├── binning.smk           # Rule: MetaBAT2 + CheckM
│   │   └── report.smk            # Rule: MultiQC report
│   └── scripts/
│       ├── filter_bins.py        # Filter bins by quality thresholds
│       └── aggregate_stats.py    # Collect per-sample QC statistics
├── envs/
│   ├── qc.yaml                   # Conda env: QC tools
│   ├── host_removal.yaml         # Conda env: Bowtie2 + SAMtools
│   ├── assembly.yaml             # Conda env: MEGAHIT
│   ├── taxonomy.yaml             # Conda env: Kraken2 + Bracken
│   ├── functional.yaml           # Conda env: HUMAnN3
│   └── binning.yaml              # Conda env: MetaBAT2 + CheckM
└── tests/
    ├── test_pipeline.py          # Unit tests for helper scripts
    └── data/
        └── README.md             # Instructions for generating test reads
```

---

## Quick Start

### Prerequisites

- [Miniconda / Mambaforge](https://conda-forge.org/miniforge/) (recommended)
- [Snakemake ≥ 8.0](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

```bash
# Install Snakemake (if not already installed)
conda create -n snakemake -c conda-forge -c bioconda snakemake
conda activate snakemake
```

### 1. Clone and configure

```bash
git clone https://github.com/hiwalihat/Zin_hiwa.git
cd Zin_hiwa
```

Edit `config/config.yaml` to set your sample names, path to raw reads, and
paths to reference databases.

### 2. Prepare input data

Raw reads are expected as gzipped FASTQ pairs:

```
data/reads/
├── sample1_R1.fastq.gz
├── sample1_R2.fastq.gz
├── sample2_R1.fastq.gz
└── sample2_R2.fastq.gz
```

Update `config/config.yaml`:

```yaml
samples:
  - sample1
  - sample2
reads_dir: "data/reads"
```

### 3. Download reference databases

| Database | Required by | Download instructions |
|----------|------------|----------------------|
| Host genome Bowtie2 index (e.g. hg38) | Host removal | [UCSC / Ensembl](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/) then `bowtie2-build hg38.fa databases/host/hg38` |
| Kraken2 standard database | Kraken2, Bracken | `kraken2-build --standard --db databases/kraken2/standard` |
| ChocoPhlAn + UniRef (HUMAnN3) | HUMAnN3 | `humann_databases --download chocophlan full databases/humann3` |
| CheckM data | CheckM | `checkm data setRoot databases/checkm` |

Update the `databases:` section of `config/config.yaml` with the correct
paths.

### 4. Dry-run to validate

```bash
snakemake --use-conda --cores all -n
```

### 5. Run the full pipeline

```bash
snakemake \
    --use-conda \
    --cores all \
    --printshellcmds \
    --rerun-incomplete
```

### 6. Run on an HPC cluster (SLURM example)

```bash
snakemake \
    --use-conda \
    --jobs 50 \
    --executor slurm \
    --default-resources slurm_partition=standard mem_mb=16000 \
    --rerun-incomplete
```

---

## Configuration Reference (`config/config.yaml`)

| Parameter | Description | Default |
|-----------|-------------|---------|
| `samples` | List of sample IDs | — |
| `reads_dir` | Directory containing raw FASTQ files | `data/reads` |
| `databases.host_genome` | Bowtie2 index prefix for host genome | — |
| `databases.kraken2_db` | Kraken2 database directory | — |
| `databases.bracken_db` | Bracken database (same as Kraken2 DB) | — |
| `databases.humann3_nucleotide_db` | HUMAnN3 ChocoPhlAn DB | — |
| `databases.humann3_protein_db` | HUMAnN3 UniRef DB | — |
| `databases.checkm_data` | CheckM data root directory | — |
| `fastp.min_length` | Minimum read length after trimming | `50` |
| `fastp.qualified_quality_phred` | Minimum base quality | `20` |
| `megahit.min_contig_len` | Minimum contig length to report | `500` |
| `kraken2.confidence` | Kraken2 confidence threshold | `0.1` |
| `bracken.taxonomic_level` | Level for Bracken re-estimation | `S` (species) |
| `checkm.completeness_threshold` | Minimum bin completeness (%) | `50.0` |
| `checkm.contamination_threshold` | Maximum bin contamination (%) | `10.0` |

---

## Output Structure

```
results/
├── qc/
│   ├── fastqc/                   # FastQC reports (raw + trimmed)
│   └── fastp/{sample}/           # fastp HTML/JSON reports + trimmed reads
├── host_removed/{sample}/        # Non-host reads (clean FASTQ)
├── assembly/{sample}/
│   ├── final.contigs.fa          # Assembled contigs
│   ├── bowtie2_index/            # Bowtie2 index of contigs
│   └── reads_vs_contigs.sorted.bam
├── taxonomy/{sample}/
│   ├── {sample}.kraken2.report   # Kraken2 classification report
│   ├── {sample}.bracken          # Bracken abundance table
│   └── {sample}.mpa_style.txt    # MetaPhlAn-style report
├── functional/{sample}/
│   ├── {sample}_genefamilies.tsv
│   ├── {sample}_pathabundance.tsv
│   ├── {sample}_pathcoverage.tsv
│   ├── {sample}_genefamilies_cpm.tsv  # Normalised (CPM)
│   └── {sample}_pathabundance_cpm.tsv
├── binning/{sample}/
│   ├── depth.txt                 # Per-contig coverage depth
│   ├── bins/                     # MetaBAT2 bin FASTA files
│   ├── checkm/
│   │   └── checkm_summary.tsv    # CheckM QA report
│   └── high_quality_bins/        # Bins passing QC thresholds
└── report/
    └── multiqc_report.html       # Aggregated MultiQC report
```

---

## Running Tests

```bash
pip install pytest pandas
pytest tests/test_pipeline.py -v
```

---

## Aggregate Statistics Script

After a completed run, generate a per-sample statistics table:

```bash
python workflow/scripts/aggregate_stats.py \
    --results_dir results \
    --samples sample1 sample2 \
    --output results/summary_stats.tsv
```

The table includes:
- Raw and trimmed read counts
- Q30 base-quality rate
- Duplication rate
- Host read percentage
- Kraken2 classified read percentage
- Number of high-quality MAGs

---

## Citation

If you use this pipeline, please cite the individual tools it wraps
(FastQC, fastp, Bowtie2, MEGAHIT, Kraken2, Bracken, HUMAnN3, MetaBAT2,
CheckM, MultiQC).

---

## License

MIT
