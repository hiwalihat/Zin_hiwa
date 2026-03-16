#!/usr/bin/env python3
"""
filter_bins.py
==============
Snakemake script called from workflow/rules/binning.smk (rule: filter_bins).

Reads the CheckM tab-delimited summary file and copies bins that pass the
user-defined completeness / contamination thresholds into a separate
high-quality directory.

Snakemake objects available
---------------------------
snakemake.input.summary  – path to CheckM QA tab-separated output
snakemake.input.bin_dir  – directory containing MetaBAT2 bins (*.fa)
snakemake.output.hq_dir  – output directory for high-quality bins
snakemake.output.hq_list – text file listing high-quality bin names
snakemake.params.completeness   – minimum completeness threshold (float)
snakemake.params.contamination  – maximum contamination threshold (float)
snakemake.log[0]         – log file path
"""

import os
import shutil
import sys
import pandas as pd

log_path = str(snakemake.log[0])

def log(msg: str) -> None:
    with open(log_path, "a") as fh:
        fh.write(msg + "\n")


def main() -> None:
    summary_path   = str(snakemake.input.summary)
    bin_dir        = str(snakemake.input.bin_dir)
    hq_dir         = str(snakemake.output.hq_dir)
    hq_list_path   = str(snakemake.output.hq_list)
    min_complete   = float(snakemake.params.completeness)
    max_contam     = float(snakemake.params.contamination)

    os.makedirs(hq_dir, exist_ok=True)

    log(f"Reading CheckM summary: {summary_path}")
    df = pd.read_csv(summary_path, sep="\t", comment="#")

    # CheckM QA -o 2 column names (may vary slightly by version)
    # Typical columns: Bin Id, Completeness, Contamination, ...
    completeness_col  = "Completeness"
    contamination_col = "Contamination"
    bin_id_col        = "Bin Id"

    # Validate expected columns exist
    for col in [bin_id_col, completeness_col, contamination_col]:
        if col not in df.columns:
            log(f"ERROR: expected column '{col}' not found. "
                f"Available columns: {list(df.columns)}")
            sys.exit(1)

    hq_df = df[
        (df[completeness_col]  >= min_complete) &
        (df[contamination_col] <= max_contam)
    ]

    log(f"Total bins: {len(df)} | High-quality bins: {len(hq_df)} "
        f"(completeness >= {min_complete}%, contamination <= {max_contam}%)")

    hq_names = []
    for bin_id in hq_df[bin_id_col]:
        src = os.path.join(bin_dir, f"{bin_id}.fa")
        dst = os.path.join(hq_dir,  f"{bin_id}.fa")
        if os.path.exists(src):
            shutil.copy2(src, dst)
            hq_names.append(bin_id)
            log(f"Copied: {bin_id}.fa")
        else:
            log(f"WARNING: bin file not found: {src}")

    with open(hq_list_path, "w") as fh:
        fh.write("\n".join(hq_names) + "\n")

    log(f"Done. {len(hq_names)} high-quality bins written to {hq_dir}")


if __name__ == "__main__":
    main()
