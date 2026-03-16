# Zin_hiwa — Illumina WGS Metagenomics Pipeline

A fully automated, **database-free-setup** pipeline for processing Illumina WGS FASTQ files for metagenomics on Ubuntu Linux.  
All required databases are downloaded **automatically** on first run — no manual download steps are needed.

---

## What the Pipeline Does

| Step | Tool | Purpose |
|------|------|---------|
| 1 | **FastQC** | Quality assessment of raw reads |
| 2 | **fastp** | Adapter trimming & quality filtering |
| 3 | **FastQC** | Quality assessment of trimmed reads |
| 4 | **MetaPhlAn 4** | Taxonomic profiling (marker-gene DB auto-downloaded ~3 GB) |
| 5 | **MultiQC** | Aggregated HTML QC report |

---

## Requirements

- Ubuntu Linux (20.04 or later recommended)
- [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda

---

## Installation

### 1 — Install Miniconda (skip if already installed)

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p "$HOME/miniconda3"
source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda init bash && source ~/.bashrc
```

### 2 — Create the conda environment

```bash
conda env create -f environment.yml
conda activate metagenomics
```

> **Note:** The MetaPhlAn marker-gene database (~3 GB) is downloaded automatically
> to `~/metaphlan_databases/` the first time Step 4 runs.
> Subsequent runs reuse the cached database — no manual download required.

---

## Usage

```
bash metagenomics_pipeline.sh -1 READ1 [-2 READ2] -s SAMPLE_NAME \
     [-o OUTPUT_DIR] [-t THREADS] [-d DB_DIR]
```

| Flag | Description | Default |
|------|-------------|---------|
| `-1 READ1` | Path to R1 FASTQ file (or single-end) | **required** |
| `-2 READ2` | Path to R2 FASTQ file (paired-end) | *(single-end)* |
| `-s SAMPLE` | Short sample identifier (no spaces) | **required** |
| `-o OUTPUT_DIR` | Output directory | `./metagenomics_output` |
| `-t THREADS` | CPU threads to use | `4` |
| `-d DB_DIR` | MetaPhlAn database directory | `~/metaphlan_databases` |

### Paired-end example

```bash
conda activate metagenomics

bash metagenomics_pipeline.sh \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    -s MySample \
    -o results \
    -t 8
```

### Single-end example

```bash
conda activate metagenomics

bash metagenomics_pipeline.sh \
    -1 sample.fastq.gz \
    -s MySample \
    -o results \
    -t 4
```

---

## Output Structure

```
results/
├── 01_fastqc_raw/          # FastQC reports on raw reads
├── 02_trimmed/             # Adapter-trimmed FASTQ files
├── 03_fastqc_trimmed/      # FastQC reports on trimmed reads
├── 04_metaphlan/
│   ├── MySample_metaphlan_profile.txt  # Taxonomic abundance table
│   └── MySample_bowtie2.bam            # Bowtie2 alignment
├── 05_multiqc/
│   └── MySample_multiqc_report.html    # All-in-one QC report
└── logs/
    ├── MySample_pipeline.log           # Full pipeline log
    ├── MySample_fastp.json
    └── MySample_fastp.html
```

---

## How "No Database Download" Works

Traditional metagenomics workflows require manually downloading and unpacking
large reference databases (often tens of GB).  
This pipeline eliminates that friction:

- **fastp** detects Illumina adapters automatically — no adapter FASTA required.
- **MetaPhlAn 4** uses a compact (~3 GB) marker-gene database that is fetched
  automatically via `metaphlan --install` on the first run and cached locally.
  All subsequent runs use the cached copy with no extra action from the user.

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `command not found` for any tool | Run `conda activate metagenomics` first |
| MetaPhlAn database download fails | Check internet connection; re-run the script — download resumes |
| Low memory during MetaPhlAn step | Reduce `-t` threads or use a machine with ≥8 GB RAM |
| Reads fail length filter | Lower `-length_required` in the fastp call inside the script |