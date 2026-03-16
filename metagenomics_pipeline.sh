#!/usr/bin/env bash
# =============================================================================
# Metagenomics Pipeline for Illumina WGS FASTQ Files
# =============================================================================
# This pipeline performs:
#   1. Quality control with FastQC
#   2. Adapter/quality trimming with fastp
#   3. Post-trim quality control with FastQC
#   4. Taxonomic profiling with MetaPhlAn 4 (database auto-downloaded on first run)
#   5. Consolidated HTML report with MultiQC
#
# Usage:
#   bash metagenomics_pipeline.sh -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \
#        -s SAMPLE_NAME -o output_directory [-t THREADS] [-d DB_DIR]
#
# Single-end mode:
#   bash metagenomics_pipeline.sh -1 sample.fastq.gz \
#        -s SAMPLE_NAME -o output_directory [-t THREADS] [-d DB_DIR]
# =============================================================================

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
THREADS=4
DB_DIR="${HOME}/metaphlan_databases"
OUTPUT_DIR="./metagenomics_output"
SAMPLE_NAME=""
READ1=""
READ2=""
PAIRED=false

# ── Usage ─────────────────────────────────────────────────────────────────────
usage() {
    cat <<EOF
Usage: $(basename "$0") -1 READ1 [-2 READ2] -s SAMPLE_NAME [-o OUTPUT_DIR] [-t THREADS] [-d DB_DIR]

Required:
  -1 READ1        Path to FASTQ file (R1 for paired-end, or single-end file).
  -s SAMPLE_NAME  A short identifier for this sample (no spaces).

Optional:
  -2 READ2        Path to R2 FASTQ file (paired-end mode).
  -o OUTPUT_DIR   Directory to write all results (default: ./metagenomics_output).
  -t THREADS      Number of CPU threads (default: 4).
  -d DB_DIR       Directory for MetaPhlAn database (default: ~/metaphlan_databases).
                  The database (~3 GB) is downloaded automatically if not present.

Examples:
  # Paired-end
  bash metagenomics_pipeline.sh -1 sample_R1.fastq.gz -2 sample_R2.fastq.gz \\
       -s MySample -o results -t 8

  # Single-end
  bash metagenomics_pipeline.sh -1 sample.fastq.gz -s MySample -o results -t 8
EOF
    exit 1
}

# ── Argument parsing ──────────────────────────────────────────────────────────
while getopts ":1:2:s:o:t:d:h" opt; do
    case $opt in
        1) READ1="$OPTARG" ;;
        2) READ2="$OPTARG" ;;
        s) SAMPLE_NAME="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        d) DB_DIR="$OPTARG" ;;
        h) usage ;;
        :) echo "[ERROR] Option -$OPTARG requires an argument." >&2; usage ;;
        \?) echo "[ERROR] Unknown option: -$OPTARG" >&2; usage ;;
    esac
done

# ── Validate required arguments ───────────────────────────────────────────────
[[ -z "$READ1" ]]        && { echo "[ERROR] -1 READ1 is required."; usage; }
[[ -z "$SAMPLE_NAME" ]]  && { echo "[ERROR] -s SAMPLE_NAME is required."; usage; }
[[ ! -f "$READ1" ]]      && { echo "[ERROR] READ1 file not found: $READ1"; exit 1; }

if [[ -n "$READ2" ]]; then
    [[ ! -f "$READ2" ]] && { echo "[ERROR] READ2 file not found: $READ2"; exit 1; }
    PAIRED=true
fi

# ── Tool availability checks ──────────────────────────────────────────────────
check_tool() {
    if ! command -v "$1" &>/dev/null; then
        echo "[ERROR] '$1' not found. Please activate the conda environment:"
        echo "         conda activate metagenomics"
        exit 1
    fi
}

for tool in fastqc fastp metaphlan multiqc; do
    check_tool "$tool"
done

# ── Directory structure ───────────────────────────────────────────────────────
FASTQC_RAW_DIR="${OUTPUT_DIR}/01_fastqc_raw"
TRIMMED_DIR="${OUTPUT_DIR}/02_trimmed"
FASTQC_TRIM_DIR="${OUTPUT_DIR}/03_fastqc_trimmed"
METAPHLAN_DIR="${OUTPUT_DIR}/04_metaphlan"
MULTIQC_DIR="${OUTPUT_DIR}/05_multiqc"
LOGS_DIR="${OUTPUT_DIR}/logs"

mkdir -p "$FASTQC_RAW_DIR" "$TRIMMED_DIR" "$FASTQC_TRIM_DIR" \
         "$METAPHLAN_DIR" "$MULTIQC_DIR" "$LOGS_DIR" "$DB_DIR"

LOG_FILE="${LOGS_DIR}/${SAMPLE_NAME}_pipeline.log"

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "$LOG_FILE"; }

log "========================================================"
log " Metagenomics Pipeline — sample: ${SAMPLE_NAME}"
log "========================================================"
log "READ1       : $READ1"
log "READ2       : ${READ2:-<single-end>}"
log "OUTPUT_DIR  : $OUTPUT_DIR"
log "DB_DIR      : $DB_DIR"
log "THREADS     : $THREADS"
log "========================================================\n"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 1 — FastQC on raw reads
# ─────────────────────────────────────────────────────────────────────────────
log "STEP 1/5 — FastQC (raw reads)"

RAW_FILES=("$READ1")
$PAIRED && RAW_FILES+=("$READ2")

fastqc --threads "$THREADS" --outdir "$FASTQC_RAW_DIR" \
       "${RAW_FILES[@]}" \
       2>>"$LOG_FILE"

log "  FastQC raw — done. Reports in: $FASTQC_RAW_DIR"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 2 — Adapter trimming + quality filtering with fastp
# ─────────────────────────────────────────────────────────────────────────────
log "STEP 2/5 — Adapter trimming with fastp"

TRIM_JSON="${LOGS_DIR}/${SAMPLE_NAME}_fastp.json"
TRIM_HTML="${LOGS_DIR}/${SAMPLE_NAME}_fastp.html"

if $PAIRED; then
    TRIM_R1="${TRIMMED_DIR}/${SAMPLE_NAME}_R1_trimmed.fastq.gz"
    TRIM_R2="${TRIMMED_DIR}/${SAMPLE_NAME}_R2_trimmed.fastq.gz"
    fastp \
        --in1  "$READ1"   --in2  "$READ2" \
        --out1 "$TRIM_R1" --out2 "$TRIM_R2" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --thread "$THREADS" \
        --json "$TRIM_JSON" \
        --html "$TRIM_HTML" \
        2>>"$LOG_FILE"
    TRIMMED_INPUT="${TRIM_R1},${TRIM_R2}"
else
    TRIM_R1="${TRIMMED_DIR}/${SAMPLE_NAME}_trimmed.fastq.gz"
    fastp \
        --in1  "$READ1" \
        --out1 "$TRIM_R1" \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --thread "$THREADS" \
        --json "$TRIM_JSON" \
        --html "$TRIM_HTML" \
        2>>"$LOG_FILE"
    TRIMMED_INPUT="$TRIM_R1"
fi

log "  fastp trimming — done. Trimmed reads in: $TRIMMED_DIR"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 3 — FastQC on trimmed reads
# ─────────────────────────────────────────────────────────────────────────────
log "STEP 3/5 — FastQC (trimmed reads)"

TRIM_QC_FILES=("$TRIM_R1")
$PAIRED && TRIM_QC_FILES+=("$TRIM_R2")

fastqc --threads "$THREADS" --outdir "$FASTQC_TRIM_DIR" \
       "${TRIM_QC_FILES[@]}" \
       2>>"$LOG_FILE"

log "  FastQC trimmed — done. Reports in: $FASTQC_TRIM_DIR"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 4 — Taxonomic profiling with MetaPhlAn 4
#
# The MetaPhlAn marker-gene database (~3 GB) is downloaded automatically to
# DB_DIR the very first time this step runs. Subsequent runs reuse the cached
# database, so no manual download is ever required.
# ─────────────────────────────────────────────────────────────────────────────
log "STEP 4/5 — MetaPhlAn 4 taxonomic profiling"
log "  Database directory: $DB_DIR"
log "  (Database will be downloaded automatically if not already present)"

METAPHLAN_PROFILE="${METAPHLAN_DIR}/${SAMPLE_NAME}_metaphlan_profile.txt"
METAPHLAN_BOWTIE2="${METAPHLAN_DIR}/${SAMPLE_NAME}_bowtie2.bam"
METAPHLAN_SAM="${METAPHLAN_DIR}/${SAMPLE_NAME}.sam.bz2"

metaphlan \
    "$TRIMMED_INPUT" \
    --input_type  fastq \
    --bowtie2db   "$DB_DIR" \
    --bowtie2out  "$METAPHLAN_BOWTIE2" \
    --samout      "$METAPHLAN_SAM" \
    --nproc       "$THREADS" \
    --output_file "$METAPHLAN_PROFILE" \
    2>>"$LOG_FILE"

log "  MetaPhlAn 4 — done. Profile: $METAPHLAN_PROFILE"

# ─────────────────────────────────────────────────────────────────────────────
# STEP 5 — Aggregate QC reports with MultiQC
# ─────────────────────────────────────────────────────────────────────────────
log "STEP 5/5 — MultiQC report"

multiqc \
    "$FASTQC_RAW_DIR" \
    "$FASTQC_TRIM_DIR" \
    "$LOGS_DIR" \
    --outdir "$MULTIQC_DIR" \
    --filename "${SAMPLE_NAME}_multiqc_report" \
    --force \
    2>>"$LOG_FILE"

log "  MultiQC — done. Report: ${MULTIQC_DIR}/${SAMPLE_NAME}_multiqc_report.html"

# ─────────────────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────────────────
log ""
log "========================================================"
log " Pipeline complete!"
log "========================================================"
log "  Outputs:"
log "    Raw FastQC        : $FASTQC_RAW_DIR"
log "    Trimmed reads     : $TRIMMED_DIR"
log "    Trimmed FastQC    : $FASTQC_TRIM_DIR"
log "    MetaPhlAn profile : $METAPHLAN_PROFILE"
log "    MultiQC report    : ${MULTIQC_DIR}/${SAMPLE_NAME}_multiqc_report.html"
log "    Pipeline log      : $LOG_FILE"
log "========================================================"
