#!/usr/bin/env python3
"""
aggregate_stats.py
==================
Collect per-sample statistics from multiple pipeline steps and write a
consolidated TSV summary table.

Usage (standalone):
    python aggregate_stats.py \
        --results_dir results \
        --samples sample1 sample2 \
        --output results/summary_stats.tsv

The script reads:
  - fastp JSON reports      → total reads, Q30 rate, duplication rate
  - Bowtie2 host-removal logs → host read percentage
  - Kraken2 reports           → classified read percentage
  - CheckM summaries          → number of HQ bins per sample
"""

import argparse
import json
import os
import re
import sys
from pathlib import Path

import pandas as pd


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

def parse_fastp_json(path: str) -> dict:
    """Extract key metrics from a fastp JSON report."""
    with open(path) as fh:
        data = json.load(fh)

    summary = data.get("summary", {})
    before  = summary.get("before_filtering", {})
    after   = summary.get("after_filtering",  {})

    return {
        "raw_reads":        before.get("total_reads", None),
        "trimmed_reads":    after.get("total_reads",  None),
        "q30_rate":         after.get("q30_rate",     None),
        "duplication_rate": data.get("duplication", {}).get("rate", None),
    }


def parse_bowtie2_log(path: str) -> dict:
    """Parse the percentage of reads that aligned to the host genome."""
    host_pct = None
    try:
        with open(path) as fh:
            for line in fh:
                m = re.search(r"([\d.]+)% overall alignment rate", line)
                if m:
                    host_pct = float(m.group(1))
    except FileNotFoundError:
        pass
    return {"host_read_pct": host_pct}


def parse_kraken2_report(path: str) -> dict:
    """Return the percentage of reads classified by Kraken2."""
    classified_pct = None
    # Kraken2 report columns (0-based):
    #   0: pct_reads  1: covered_reads  2: direct_reads
    #   3: rank_code  4: taxid          5: name
    KRAKEN2_PCT_COL  = 0
    KRAKEN2_NAME_COL = 5
    KRAKEN2_MIN_COLS = 6
    try:
        with open(path) as fh:
            for line in fh:
                parts = line.strip().split("\t")
                if (len(parts) >= KRAKEN2_MIN_COLS and
                        parts[KRAKEN2_NAME_COL].strip() == "root"):
                    classified_pct = float(parts[KRAKEN2_PCT_COL])
                    break
    except FileNotFoundError:
        pass
    return {"kraken2_classified_pct": classified_pct}


def parse_checkm_summary(path: str) -> dict:
    """Count the number of high-quality bins listed in the CheckM summary."""
    try:
        df = pd.read_csv(path, sep="\t", comment="#")
        return {"hq_bins": len(df)}
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return {"hq_bins": 0}


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def build_table(results_dir: str, samples: list) -> pd.DataFrame:
    rows = []
    for sample in samples:
        row = {"sample": sample}

        fastp_json = os.path.join(
            results_dir, "qc", "fastp", sample, f"{sample}_fastp.json"
        )
        bowtie2_log = os.path.join(
            results_dir, "host_removed", sample, f"{sample}_bowtie2.log"
        )
        kraken2_report = os.path.join(
            results_dir, "taxonomy", sample, f"{sample}.kraken2.report"
        )
        checkm_summary = os.path.join(
            results_dir, "binning", sample, "checkm", "checkm_summary.tsv"
        )

        if os.path.exists(fastp_json):
            row.update(parse_fastp_json(fastp_json))
        if os.path.exists(bowtie2_log):
            row.update(parse_bowtie2_log(bowtie2_log))
        if os.path.exists(kraken2_report):
            row.update(parse_kraken2_report(kraken2_report))
        if os.path.exists(checkm_summary):
            row.update(parse_checkm_summary(checkm_summary))

        rows.append(row)

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results_dir", required=True,
                        help="Pipeline results directory")
    parser.add_argument("--samples", nargs="+", required=True,
                        help="List of sample IDs")
    parser.add_argument("--output", required=True,
                        help="Output TSV file path")
    args = parser.parse_args()

    df = build_table(args.results_dir, args.samples)
    df.to_csv(args.output, sep="\t", index=False)
    print(f"Summary written to {args.output}")


if __name__ == "__main__":
    main()
