# Tools and Packages Reference

## Comprehensive Tool List for Stool Microbiome WGS Analysis Pipeline

### Linux Bioinformatics Tools

| Pipeline Stage | Tool/Package | Version | Purpose | Platform | Installation/Source |
|----------------|-------------|---------|---------|----------|-------------------|
| **Quality Control** | FastQC | v0.11.9+ | Quality assessment of raw sequencing reads, generates per-base quality scores, GC content, adapter content | Linux | [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) |
| **Quality Control** | MultiQC | v1.14+ | Aggregates QC reports from multiple samples into single interactive report | Linux/Python | [https://multiqc.info/](https://multiqc.info/) `pip install multiqc` |
| **Quality Control** | fastp | v0.23.2+ | Ultra-fast FASTQ preprocessing: adapter trimming, quality filtering, length filtering | Linux | [https://github.com/OpenGene/fastp](https://github.com/OpenGene/fastp) `conda install -c bioconda fastp` |
| **Quality Control** | cutadapt | v4.1+ | Alternative trimming tool for adapter removal and quality filtering | Linux/Python | [https://cutadapt.readthedocs.io/](https://cutadapt.readthedocs.io/) `pip install cutadapt` |
| **Host Decontamination** | bowtie2 | v2.4.5+ | Fast aligner for mapping reads to human reference genome for removal | Linux | [http://bowtie-bio.sourceforge.net/bowtie2/](http://bowtie-bio.sourceforge.net/bowtie2/) `conda install -c bioconda bowtie2` |
| **Host Decontamination** | kneaddata | v0.12+ | Integrated tool for quality control and human read removal (uses bowtie2 internally) | Linux/Python | [https://huttenhower.sph.harvard.edu/kneaddata/](https://huttenhower.sph.harvard.edu/kneaddata/) `pip install kneaddata` |
| **Taxonomic Classification** | Kraken2 | v2.1.2+ | Fast k-mer based taxonomic classification system (bacteria-focused for this pipeline) | Linux | [https://ccb.jhu.edu/software/kraken2/](https://ccb.jhu.edu/software/kraken2/) `conda install -c bioconda kraken2` |
| **Abundance Estimation** | Bracken | v2.7+ | Bayesian reestimation of abundance after classification (phylum to species level) | Linux | [https://ccb.jhu.edu/software/bracken/](https://ccb.jhu.edu/software/bracken/) `conda install -c bioconda bracken` |

### R Statistical Packages

| Pipeline Stage | Tool/Package | Version | Purpose | Platform | Installation/Source |
|----------------|-------------|---------|---------|----------|-------------------|
| **Data Management** | phyloseq | v1.42+ | Core package for microbiome data import, manipulation, and analysis | R/Bioconductor | [https://joey711.github.io/phyloseq/](https://joey711.github.io/phyloseq/) `BiocManager::install("phyloseq")` |
| **Data Wrangling** | tidyverse | v2.0+ | Collection of data manipulation packages (dplyr, tidyr, readr, etc.) | R/CRAN | [https://www.tidyverse.org/](https://www.tidyverse.org/) `install.packages("tidyverse")` |
| **Diversity Analysis** | vegan | v2.6+ | Community ecology package for diversity indices and ordination (PCoA, NMDS, PERMANOVA) | R/CRAN | [https://github.com/vegandevs/vegan](https://github.com/vegandevs/vegan) `install.packages("vegan")` |
| **Differential Abundance** | ANCOMBC | v2.0+ | Analysis of Compositions of Microbiomes with Bias Correction for differential abundance testing | R/Bioconductor | [https://bioconductor.org/packages/ANCOMBC/](https://bioconductor.org/packages/ANCOMBC/) `BiocManager::install("ANCOMBC")` |
| **Visualization** | ggplot2 | v3.4+ | Grammar of graphics plotting system (included in tidyverse) | R/CRAN | [https://ggplot2.tidyverse.org/](https://ggplot2.tidyverse.org/) `install.packages("ggplot2")` |
| **Microbiome Utilities** | microbiome | v1.20+ | Tools for microbiome analysis: diversity, composition, core microbiome analysis | R/Bioconductor | [https://microbiome.github.io/](https://microbiome.github.io/) `BiocManager::install("microbiome")` |
| **Compositional Analysis** | compositions | v2.0+ | Compositional data analysis, centered log-ratio (CLR) transformation | R/CRAN | [https://cran.r-project.org/package=compositions](https://cran.r-project.org/package=compositions) `install.packages("compositions")` |
| **Reporting** | Quarto | v1.3+ | Next-generation scientific publishing system for reproducible reports | R/Standalone | [https://quarto.org/](https://quarto.org/) Download from website |
| **Reporting** | RMarkdown | v2.20+ | Dynamic document generation system (alternative to Quarto) | R/CRAN | [https://rmarkdown.rstudio.com/](https://rmarkdown.rstudio.com/) `install.packages("rmarkdown")` |

## Additional Recommended Packages

### Optional R Packages for Enhanced Analysis

| Package | Purpose | Installation |
|---------|---------|-------------|
| **DESeq2** | Alternative differential abundance method (originally for RNA-seq) | `BiocManager::install("DESeq2")` |
| **ComplexHeatmap** | Advanced heatmap visualization with hierarchical clustering | `BiocManager::install("ComplexHeatmap")` |
| **patchwork** | Compose multiple ggplot2 plots into publication-ready figures | `install.packages("patchwork")` |
| **ggtree** | Phylogenetic tree visualization | `BiocManager::install("ggtree")` |
| **scales** | Scale functions for visualization (log scales, percentages) | `install.packages("scales")` |
| **RColorBrewer** | Color palettes for publication graphics | `install.packages("RColorBrewer")` |

## Database Requirements

| Resource | Purpose | Source |
|----------|---------|--------|
| **Kraken2 Standard Database** | Bacterial, archaeal, viral reference genomes for classification | [https://benlangmead.github.io/aws-indexes/k2](https://benlangmead.github.io/aws-indexes/k2) |
| **Human Reference Genome (GRCh38)** | For host read decontamination with bowtie2/kneaddata | [https://www.ncbi.nlm.nih.gov/genome/guide/human/](https://www.ncbi.nlm.nih.gov/genome/guide/human/) |
| **Bracken Database** | Pre-built from same database as Kraken2 for abundance estimation | Built using `bracken-build` after Kraken2 database |

## Software Versions Note

Version numbers provided are minimum recommended versions as of 2024-2025. Always check for the latest stable releases and compatibility with your system. Use package managers (conda, pip, CRAN, Bioconductor) to manage dependencies and ensure reproducibility.

## Platform-Specific Notes

- **Linux tools**: Best installed via conda/bioconda for dependency management
- **R packages**: Use `install.packages()` for CRAN, `BiocManager::install()` for Bioconductor
- **Python tools**: Use pip in isolated environments or conda for better reproducibility
