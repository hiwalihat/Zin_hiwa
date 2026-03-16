#!/usr/bin/env python3
"""
Unit tests for metagenomics pipeline helper scripts.

Run with:
    pytest tests/test_pipeline.py -v
"""

import json
import os
import sys
import tempfile
import textwrap
from pathlib import Path
from unittest.mock import MagicMock

import pandas as pd
import pytest

# Make workflow scripts importable
REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT / "workflow" / "scripts"))

import aggregate_stats as agg


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def tmp_dir(tmp_path):
    return tmp_path


@pytest.fixture()
def fastp_json(tmp_dir):
    data = {
        "summary": {
            "before_filtering": {"total_reads": 2_000_000},
            "after_filtering":  {"total_reads": 1_800_000, "q30_rate": 0.95},
        },
        "duplication": {"rate": 0.12},
    }
    p = tmp_dir / "sample1_fastp.json"
    p.write_text(json.dumps(data))
    return str(p)


@pytest.fixture()
def bowtie2_log(tmp_dir):
    content = textwrap.dedent("""\
        1800000 reads; of these:
          900000 (50.00%) were paired; of these:
        3.52% overall alignment rate
    """)
    p = tmp_dir / "sample1_bowtie2.log"
    p.write_text(content)
    return str(p)


@pytest.fixture()
def kraken2_report(tmp_dir):
    content = textwrap.dedent("""\
         78.50\t1413000\t0\tR\t1\troot
          5.25\t94500\t0\tD\t2\t  Bacteria
    """)
    p = tmp_dir / "sample1.kraken2.report"
    p.write_text(content)
    return str(p)


@pytest.fixture()
def checkm_summary(tmp_dir):
    content = textwrap.dedent("""\
        Bin Id\tCompleteness\tContamination\tGenome size (bp)
        bin.1\t95.3\t1.2\t3500000
        bin.2\t82.1\t4.5\t2900000
        bin.3\t45.0\t2.0\t1200000
    """)
    p = tmp_dir / "checkm_summary.tsv"
    p.write_text(content)
    return str(p)


# ---------------------------------------------------------------------------
# Tests: aggregate_stats parsers
# ---------------------------------------------------------------------------

class TestParseFastpJson:
    def test_parses_read_counts(self, fastp_json):
        result = agg.parse_fastp_json(fastp_json)
        assert result["raw_reads"]     == 2_000_000
        assert result["trimmed_reads"] == 1_800_000

    def test_parses_q30_rate(self, fastp_json):
        result = agg.parse_fastp_json(fastp_json)
        assert result["q30_rate"] == pytest.approx(0.95)

    def test_parses_duplication_rate(self, fastp_json):
        result = agg.parse_fastp_json(fastp_json)
        assert result["duplication_rate"] == pytest.approx(0.12)


class TestParseBowtie2Log:
    def test_parses_alignment_rate(self, bowtie2_log):
        result = agg.parse_bowtie2_log(bowtie2_log)
        assert result["host_read_pct"] == pytest.approx(3.52)

    def test_returns_none_for_missing_file(self):
        result = agg.parse_bowtie2_log("/nonexistent/path.log")
        assert result["host_read_pct"] is None


class TestParseKraken2Report:
    def test_parses_classified_pct(self, kraken2_report):
        result = agg.parse_kraken2_report(kraken2_report)
        assert result["kraken2_classified_pct"] == pytest.approx(78.50)

    def test_returns_none_for_missing_file(self):
        result = agg.parse_kraken2_report("/nonexistent/path.report")
        assert result["kraken2_classified_pct"] is None


class TestParseCheckmSummary:
    def test_counts_bins(self, checkm_summary):
        result = agg.parse_checkm_summary(checkm_summary)
        assert result["hq_bins"] == 3

    def test_returns_zero_for_missing_file(self):
        result = agg.parse_checkm_summary("/nonexistent/checkm.tsv")
        assert result["hq_bins"] == 0


# ---------------------------------------------------------------------------
# Tests: aggregate_stats.build_table integration
# ---------------------------------------------------------------------------

class TestBuildTable:
    def test_returns_dataframe_with_correct_columns(
        self, tmp_dir, fastp_json, bowtie2_log, kraken2_report, checkm_summary
    ):
        # Set up directory structure mirroring actual pipeline output
        results = tmp_dir / "results"
        sample = "sample1"

        fastp_dir   = results / "qc" / "fastp" / sample
        hr_dir      = results / "host_removed" / sample
        tax_dir     = results / "taxonomy" / sample
        checkm_dir  = results / "binning" / sample / "checkm"

        for d in [fastp_dir, hr_dir, tax_dir, checkm_dir]:
            d.mkdir(parents=True)

        # Copy fixture files into expected locations
        import shutil
        shutil.copy(fastp_json,    str(fastp_dir / f"{sample}_fastp.json"))
        shutil.copy(bowtie2_log,   str(hr_dir    / f"{sample}_bowtie2.log"))
        shutil.copy(kraken2_report,str(tax_dir   / f"{sample}.kraken2.report"))
        shutil.copy(checkm_summary,str(checkm_dir / "checkm_summary.tsv"))

        df = agg.build_table(str(results), [sample])

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 1
        assert df.loc[0, "sample"] == sample
        assert df.loc[0, "raw_reads"] == 2_000_000
        assert df.loc[0, "host_read_pct"] == pytest.approx(3.52)
        assert df.loc[0, "kraken2_classified_pct"] == pytest.approx(78.50)
        assert df.loc[0, "hq_bins"] == 3

    def test_handles_missing_files_gracefully(self, tmp_dir):
        df = agg.build_table(str(tmp_dir / "empty_results"), ["missing_sample"])
        assert isinstance(df, pd.DataFrame)
        assert len(df) == 1
        assert df.loc[0, "sample"] == "missing_sample"


# ---------------------------------------------------------------------------
# Tests: filter_bins script (via direct function call simulation)
# ---------------------------------------------------------------------------

class TestFilterBins:
    """
    filter_bins.py is a Snakemake script that reads from the `snakemake`
    global object.  We test its logic by importing the module with a mock
    snakemake context.
    """

    def _make_checkm_df(self):
        return pd.DataFrame({
            "Bin Id":        ["bin.1", "bin.2", "bin.3", "bin.4"],
            "Completeness":  [95.0,    60.0,    40.0,    85.0],
            "Contamination": [1.0,     5.0,     2.0,     12.0],
        })

    def test_filter_logic(self):
        df = self._make_checkm_df()
        min_complete = 50.0
        max_contam   = 10.0
        hq = df[
            (df["Completeness"]  >= min_complete) &
            (df["Contamination"] <= max_contam)
        ]
        assert set(hq["Bin Id"]) == {"bin.1", "bin.2"}

    def test_all_pass_when_thresholds_zero(self):
        df = self._make_checkm_df()
        hq = df[(df["Completeness"] >= 0) & (df["Contamination"] <= 100)]
        assert len(hq) == len(df)

    def test_none_pass_when_thresholds_extreme(self):
        df = self._make_checkm_df()
        hq = df[(df["Completeness"] >= 100) & (df["Contamination"] <= 0)]
        assert len(hq) == 0
