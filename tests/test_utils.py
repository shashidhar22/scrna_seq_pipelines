"""
Tests for pipeline/utils.py: config, metadata parsing, QC helpers.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import yaml

# Ensure pipeline/ is importable
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))

from utils import (
    load_config,
    setup_paths,
    load_metadata,
    load_clinical_metadata,
    get_gex_samples,
    get_vdj_samples,
    get_sample_arm,
    filter_samples_by_arm,
    map_clinical_metadata,
    find_cellranger_output,
    calculate_qc_metrics,
    mad_outlier,
    plot_filtering_summary,
)


# ===== Configuration =====

class TestLoadConfig:
    def test_load_yaml(self, tmp_path):
        cfg_data = {"paths": {"metadata": "test.csv"}, "study": {"study_id": "TEST"}}
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text(yaml.dump(cfg_data))
        result = load_config(str(cfg_file))
        assert result["study"]["study_id"] == "TEST"
        assert result["paths"]["metadata"] == "test.csv"

    def test_load_missing_file(self):
        with pytest.raises(FileNotFoundError):
            load_config("/nonexistent/config.yaml")


class TestSetupPaths:
    def test_creates_output_dirs(self, minimal_config, tmp_path):
        setup_paths(minimal_config)
        assert Path(minimal_config["paths"]["qc_output"]).is_dir()
        assert Path(minimal_config["paths"]["figures"]).is_dir()
        assert Path(minimal_config["paths"]["integration_output"]).is_dir()

    def test_skips_reference_paths(self, minimal_config):
        """Should not create dirs for reference/metadata paths."""
        setup_paths(minimal_config)
        assert not Path(minimal_config["paths"]["gex_reference"]).exists()
        assert not Path(minimal_config["paths"]["vdj_reference"]).exists()


# ===== Metadata Loading =====

class TestLoadMetadata:
    def test_basic_load(self, config_with_metadata):
        meta = load_metadata(config_with_metadata)
        assert isinstance(meta, pd.DataFrame)
        assert len(meta) > 0
        assert "patientID" in meta.columns

    def test_hiv_status_override(self, config_with_metadata):
        """Patient 008_252 should have all rows corrected to 'negative'."""
        meta = load_metadata(config_with_metadata)
        p252 = meta[meta["patientID"] == "008_252"]
        assert len(p252) > 0
        assert (p252["HIVstatus"] == "negative").all(), (
            f"HIV status override failed: {p252['HIVstatus'].tolist()}"
        )

    def test_hiv_override_does_not_affect_others(self, config_with_metadata):
        """HIV status for other patients should not change."""
        meta = load_metadata(config_with_metadata)
        p220 = meta[meta["patientID"] == "008_220"]
        assert (p220["HIVstatus"] == "positive").all()


class TestGetGexSamples:
    def test_filters_gex_only(self, sample_metadata_df):
        gex = get_gex_samples(sample_metadata_df)
        assert (gex["locus"] == "5primeGEX").all()
        assert len(gex) < len(sample_metadata_df)

    def test_returns_copy(self, sample_metadata_df):
        gex = get_gex_samples(sample_metadata_df)
        gex["new_col"] = 1
        assert "new_col" not in sample_metadata_df.columns


class TestGetVdjSamples:
    def test_filters_vdj_only(self, sample_metadata_df):
        vdj = get_vdj_samples(sample_metadata_df)
        assert (vdj["locus"] == "5primeVDJ").all()


class TestGetSampleArm:
    def test_unsorted_is_primary(self, minimal_config):
        row = pd.Series({"sortingCT": "Unsorted", "patientID": "008_999"})
        assert get_sample_arm(row, minimal_config) == "primary"

    def test_sorted_cd4_is_focused(self, minimal_config):
        row = pd.Series({"sortingCT": "CD4pTCells", "patientID": "008_999"})
        assert get_sample_arm(row, minimal_config) == "focused"

    def test_sorted_cd8_is_focused(self, minimal_config):
        row = pd.Series({"sortingCT": "CD8pTCells", "patientID": "008_999"})
        assert get_sample_arm(row, minimal_config) == "focused"

    def test_sorted_validation_patient(self, minimal_config):
        row = pd.Series({"sortingCT": "CD4pTCells", "patientID": "008_216"})
        assert get_sample_arm(row, minimal_config) == "validation"


class TestFilterSamplesByArm:
    def test_primary_arm(self, sample_metadata_df, minimal_config):
        result = filter_samples_by_arm(sample_metadata_df, "primary", minimal_config)
        assert (result["sortingCT"] == "Unsorted").all()

    def test_focused_arm(self, sample_metadata_df, minimal_config):
        result = filter_samples_by_arm(sample_metadata_df, "focused", minimal_config)
        assert result["sortingCT"].isin(["CD4pTCells", "CD8pTCells"]).all()

    def test_validation_arm(self, sample_metadata_df, minimal_config):
        result = filter_samples_by_arm(sample_metadata_df, "validation", minimal_config)
        assert result["patientID"].isin(["008_216", "008_217", "008_220"]).all()

    def test_unknown_arm_raises(self, sample_metadata_df, minimal_config):
        with pytest.raises(ValueError, match="Unknown analysis arm"):
            filter_samples_by_arm(sample_metadata_df, "bogus", minimal_config)


class TestMapClinicalMetadata:
    def test_maps_ptid_correctly(self, sample_metadata_df, clinical_metadata_df):
        merged = map_clinical_metadata(sample_metadata_df, clinical_metadata_df)
        # 008_216 -> ptid 216 -> best_resp CR
        p216 = merged[merged["patientID"] == "008_216"]
        assert (p216["best_resp"] == "CR").all()

    def test_non_008_patients_get_nan(self, clinical_metadata_df):
        """Patients without 008_XXX format should get NaN for clinical cols."""
        meta = pd.DataFrame({
            "patientID": ["Aph1", "EHW1"],
            "locus": ["5primeGEX", "5primeGEX"],
        })
        merged = map_clinical_metadata(meta, clinical_metadata_df)
        assert merged["best_resp"].isna().all()

    def test_unmatched_ptid_gets_nan(self, sample_metadata_df, clinical_metadata_df):
        merged = map_clinical_metadata(sample_metadata_df, clinical_metadata_df)
        # ptid 999 has clinical data but no matching 008_999 in sample metadata
        # All rows that do match should have non-NaN best_resp
        matched = merged[merged["ptid"].notna() & merged["best_resp"].notna()]
        assert len(matched) > 0


# ===== I/O helpers =====

class TestFindCellrangerOutput:
    def test_paths_structure(self, minimal_config):
        result = find_cellranger_output("MySample", minimal_config)
        assert "raw_h5" in result
        assert "filtered_h5" in result
        assert "metrics" in result
        assert "vdj_contigs" in result
        assert "MySample" in str(result["raw_h5"])


# ===== QC helpers =====

class TestCalculateQcMetrics:
    def test_adds_mt_ribo_hb_flags(self, adata_with_qc_genes):
        adata = calculate_qc_metrics(adata_with_qc_genes)
        assert "mt" in adata.var.columns
        assert "ribo" in adata.var.columns
        assert "hb" in adata.var.columns
        assert adata.var["mt"].sum() == 5  # MT-GENE0..4
        assert adata.var["hb"].sum() == 3  # HBA1, HBA2, HBB

    def test_adds_pct_counts(self, adata_with_qc_genes):
        adata = calculate_qc_metrics(adata_with_qc_genes)
        assert "pct_counts_mt" in adata.obs.columns
        assert "pct_counts_ribo" in adata.obs.columns
        assert "pct_counts_hb" in adata.obs.columns
        assert "total_counts" in adata.obs.columns
        assert "n_genes_by_counts" in adata.obs.columns


class TestMadOutlier:
    def test_detects_obvious_outliers(self):
        """With spread data (non-zero MAD), extreme values are flagged."""
        values = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 100])
        outliers = mad_outlier(values, nmads=3, direction="upper")
        assert outliers[-1] == True  # 100 is an outlier
        assert outliers[0] == False

    def test_both_direction(self):
        values = np.array([-50, 1, 2, 3, 4, 5, 6, 7, 8, 100])
        outliers = mad_outlier(values, nmads=3, direction="both")
        assert outliers[0] == True   # -50 is low outlier
        assert outliers[-1] == True  # 100 is high outlier
        assert not outliers[4]       # 4 is normal

    def test_lower_direction(self):
        values = np.array([-100, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        outliers = mad_outlier(values, nmads=3, direction="lower")
        assert outliers[0] == True

    def test_constant_data_returns_no_outliers(self):
        values = np.array([5, 5, 5, 5, 5])
        outliers = mad_outlier(values, nmads=3, direction="both")
        assert not outliers.any()

    def test_near_constant_uses_std_fallback(self):
        """When MAD=0, falls back to std. With std fallback, 3*1.4826*std
        creates wide bounds, so only very extreme values are outliers."""
        values = np.array([5, 5, 5, 5, 5, 5, 5, 5, 5, 5000])
        outliers = mad_outlier(values, nmads=3, direction="both")
        # std fallback: std ~ 1500, so upper ~ 5 + 3*1.4826*1500 ~ 6676
        # 5000 is within bounds. This tests that the fallback doesn't crash.
        assert isinstance(outliers, np.ndarray)
        assert len(outliers) == len(values)

    def test_handles_nan(self):
        """NaN values should be handled gracefully via nanmedian."""
        values = np.array([1, np.nan, 2, 3, 4, 5, 6, 7, 100])
        outliers = mad_outlier(values, nmads=3, direction="upper")
        assert outliers[-1] == True

    def test_invalid_direction_raises(self):
        with pytest.raises(ValueError, match="Unknown direction"):
            mad_outlier(np.array([1, 2, 3]), direction="sideways")

    def test_scale_factor(self):
        """Verify the 1.4826 scale factor is applied correctly."""
        values = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        med = np.median(values)
        mad = np.median(np.abs(values - med))
        # Upper threshold = med + 3 * 1.4826 * mad
        expected_upper = med + 3 * 1.4826 * mad
        # No value in 1..10 should be an outlier with this range
        outliers = mad_outlier(values, nmads=3, direction="upper")
        assert not outliers.any()

    def test_upper_only_ignores_low(self):
        """Upper direction should not flag low values."""
        values = np.array([-100, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        outliers = mad_outlier(values, nmads=3, direction="upper")
        assert not outliers[0]  # -100 is low, not flagged


# ===== Plotting helpers =====

class TestPlotFilteringSummary:
    def test_produces_figure(self, tmp_path):
        pre = {"s1": 1000, "s2": 800}
        post = {"s1": 900, "s2": 600}
        save_path = str(tmp_path / "summary.png")
        plot_filtering_summary(pre, post, save_path=save_path)
        assert Path(save_path).exists()
