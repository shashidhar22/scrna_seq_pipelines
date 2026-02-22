"""
Tests for pipeline/00_cellranger_multi.py: FASTQ resolution, sample table, config generation.
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))

from conftest import PIPELINE_DIR  # noqa: F401  (fixture)

# Import the module under test
import importlib.util

spec = importlib.util.spec_from_file_location(
    "cellranger_multi",
    str(Path(__file__).resolve().parents[1] / "pipeline" / "00_cellranger_multi.py"),
)
cr_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cr_mod)

resolve_fastq_paths = cr_mod.resolve_fastq_paths
build_sample_table = cr_mod.build_sample_table
write_multi_config = cr_mod.write_multi_config
verify_cellranger_outputs = cr_mod.verify_cellranger_outputs


# ===== resolve_fastq_paths =====

class TestResolveFastqPaths:
    def test_finds_standard_path(self, tmp_path):
        """Standard mkfastq output: {base}/{batch}/outs/fastq_path/{batch}/{repoID}/"""
        fastq_dir = tmp_path / "BATCH1" / "outs" / "fastq_path" / "BATCH1" / "REPO1"
        fastq_dir.mkdir(parents=True)
        result = resolve_fastq_paths("REPO1", "BATCH1", str(tmp_path))
        assert len(result) == 1
        assert "REPO1" in result[0]

    def test_finds_flat_path(self, tmp_path):
        """Flat structure: {base}/{batch}/{repoID}/"""
        fastq_dir = tmp_path / "BATCH1" / "REPO1"
        fastq_dir.mkdir(parents=True)
        result = resolve_fastq_paths("REPO1", "BATCH1", str(tmp_path))
        assert len(result) >= 1

    def test_no_match_returns_empty(self, tmp_path):
        result = resolve_fastq_paths("REPO_MISSING", "BATCH1", str(tmp_path))
        assert result == []

    def test_multiple_matches(self, tmp_path):
        """Both standard and flat paths exist."""
        (tmp_path / "B" / "outs" / "fastq_path" / "B" / "R").mkdir(parents=True)
        (tmp_path / "B" / "R").mkdir(parents=True)
        result = resolve_fastq_paths("R", "B", str(tmp_path))
        assert len(result) == 2


# ===== build_sample_table =====

class TestBuildSampleTable:
    def test_pairs_gex_and_vdj(self, sample_metadata_df, minimal_config, tmp_path):
        """GEX + VDJ rows for the same sampleName should be paired."""
        # Create fake FASTQ dirs so paths resolve
        for row in sample_metadata_df.itertuples():
            if row.Batch != "AAAMMK2M5":
                d = tmp_path / "fastq" / row.Batch / row.repoID
                d.mkdir(parents=True, exist_ok=True)
        minimal_config["paths"]["fastq_base"] = str(tmp_path / "fastq")

        table = build_sample_table(sample_metadata_df, minimal_config)
        assert isinstance(table, pd.DataFrame)
        assert "sample_name" in table.columns
        assert "has_gex" in table.columns
        assert "has_vdj" in table.columns

        # 008_216_V01 has both GEX and VDJ
        row_216 = table[table["sample_name"] == "008_216_V01"]
        assert len(row_216) == 1
        assert row_216.iloc[0]["has_gex"]
        assert row_216.iloc[0]["has_vdj"]

    def test_skips_missing_batch(self, sample_metadata_df, minimal_config, tmp_path):
        """Rows in AAAMMK2M5 should be skipped."""
        minimal_config["paths"]["fastq_base"] = str(tmp_path / "fastq")
        table = build_sample_table(sample_metadata_df, minimal_config)
        # 008_216_V09 only has a row in AAAMMK2M5, so it should not have GEX fastqs
        row_v09 = table[table["sample_name"] == "008_216_V09"]
        if len(row_v09) > 0:
            assert not row_v09.iloc[0]["has_gex"]

    def test_gex_only_sample(self, sample_metadata_df, minimal_config, tmp_path):
        """008_220_V01 has GEX but no VDJ in the test fixture."""
        for row in sample_metadata_df.itertuples():
            if row.Batch != "AAAMMK2M5":
                d = tmp_path / "fastq" / row.Batch / row.repoID
                d.mkdir(parents=True, exist_ok=True)
        minimal_config["paths"]["fastq_base"] = str(tmp_path / "fastq")

        table = build_sample_table(sample_metadata_df, minimal_config)
        row_220 = table[table["sample_name"] == "008_220_V01"]
        assert len(row_220) == 1
        assert row_220.iloc[0]["has_gex"]
        assert not row_220.iloc[0]["has_vdj"]


# ===== write_multi_config =====

class TestWriteMultiConfig:
    def test_writes_valid_config(self, minimal_config, tmp_path):
        sample_row = pd.Series({
            "sample_name": "SAMPLE_A",
            "patient_id": "008_100",
            "expected_cells": 3000,
            "gex_fastqs": ["/path/to/gex"],
            "gex_fastq_ids": ["GEX_REPO"],
            "vdj_fastqs": ["/path/to/vdj"],
            "vdj_fastq_ids": ["VDJ_REPO"],
            "has_gex": True,
            "has_vdj": True,
        })
        output_dir = str(tmp_path / "configs")
        path = write_multi_config(sample_row, minimal_config, output_dir)

        content = Path(path).read_text()
        assert "[gene-expression]" in content
        assert "[vdj]" in content
        assert "[libraries]" in content
        assert "Gene Expression" in content
        assert "VDJ-T" in content
        assert "GEX_REPO" in content
        assert "3000" in content

    def test_gex_only_config(self, minimal_config, tmp_path):
        sample_row = pd.Series({
            "sample_name": "SAMPLE_B",
            "patient_id": "008_200",
            "expected_cells": 3000,
            "gex_fastqs": ["/gex"],
            "gex_fastq_ids": ["GEX_ID"],
            "vdj_fastqs": [],
            "vdj_fastq_ids": [],
            "has_gex": True,
            "has_vdj": False,
        })
        path = write_multi_config(sample_row, minimal_config, str(tmp_path))
        content = Path(path).read_text()
        assert "[gene-expression]" in content
        assert "[vdj]" not in content


# ===== verify_cellranger_outputs =====

class TestVerifyCellrangerOutputs:
    def test_missing_outputs_flagged(self, minimal_config, sample_metadata_df):
        """All samples should be MISSING when no Cell Ranger output exists."""
        result = verify_cellranger_outputs(minimal_config, sample_metadata_df)
        assert isinstance(result, pd.DataFrame)
        assert (result["status"] == "MISSING").all()

    def test_ok_sample(self, minimal_config, sample_metadata_df, tmp_path):
        """A sample with good metrics should be OK."""
        sample = "008_216_V01"
        metrics_dir = Path(minimal_config["paths"]["cellranger_output"]) / sample / "outs" / "multi" / "count"
        metrics_dir.mkdir(parents=True)
        metrics = pd.DataFrame({
            "Estimated Number of Cells": [3000],
            "Median Genes per Cell": [1500],
            "Sequencing Saturation": ["75%"],
        })
        metrics.to_csv(metrics_dir / "summary.csv", index=False)

        result = verify_cellranger_outputs(minimal_config, sample_metadata_df)
        row = result[result["sample"] == sample]
        assert len(row) == 1
        assert row.iloc[0]["status"] == "OK"

    def test_low_cells_flagged(self, minimal_config, sample_metadata_df, tmp_path):
        """Sample with <1000 cells should be flagged."""
        sample = "008_216_V01"
        metrics_dir = Path(minimal_config["paths"]["cellranger_output"]) / sample / "outs" / "multi" / "count"
        metrics_dir.mkdir(parents=True)
        metrics = pd.DataFrame({
            "Estimated Number of Cells": [500],
            "Median Genes per Cell": [1500],
            "Sequencing Saturation": ["75%"],
        })
        metrics.to_csv(metrics_dir / "summary.csv", index=False)

        result = verify_cellranger_outputs(minimal_config, sample_metadata_df)
        row = result[result["sample"] == sample]
        assert row.iloc[0]["status"] == "FLAG"
        assert "LOW_CELLS" in row.iloc[0].get("flags", "")
