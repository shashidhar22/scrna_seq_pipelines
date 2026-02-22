"""
Tests against the real KSTME metadata files.
These tests verify actual data properties and known issues.
Skipped if data files are not present.
"""

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))

REPO_ROOT = Path(__file__).resolve().parents[1]
META_PATH = REPO_ROOT / "data" / "kstme_sc_metadata.csv"
CLINICAL_PATH = REPO_ROOT / "data" / "KSTME_updated_clinical_metadata.csv"
CONFIG_PATH = REPO_ROOT / "pipeline" / "config.yaml"

skip_if_no_data = pytest.mark.skipif(
    not META_PATH.exists(), reason="Real metadata not available"
)
skip_if_no_clinical = pytest.mark.skipif(
    not CLINICAL_PATH.exists(), reason="Real clinical metadata not available"
)
skip_if_no_config = pytest.mark.skipif(
    not CONFIG_PATH.exists(), reason="Pipeline config not available"
)


@skip_if_no_data
class TestRealMetadata:
    @pytest.fixture(autouse=True)
    def _load(self):
        self.meta = pd.read_csv(META_PATH)

    def test_expected_columns(self):
        expected = [
            "repoName", "sampleName", "patientID", "visitCode",
            "expected_cells", "locus", "sortingCT", "tissueType",
            "Batch", "HIVstatus", "repoID",
        ]
        for col in expected:
            assert col in self.meta.columns, f"Missing column: {col}"

    def test_row_count(self):
        # 100 rows: 50 GEX + 50 VDJ
        assert len(self.meta) == 100

    def test_locus_values(self):
        assert set(self.meta["locus"].unique()) == {"5primeGEX", "5primeVDJ"}

    def test_gex_count(self):
        gex = self.meta[self.meta["locus"] == "5primeGEX"]
        assert len(gex) == 50

    def test_vdj_count(self):
        vdj = self.meta[self.meta["locus"] == "5primeVDJ"]
        assert len(vdj) == 50

    def test_known_batches(self):
        known = {"H2JJHDRX2", "HTVWKDRXX", "HGMLLDRXX", "HVFNLDRXX",
                 "AAAMMK2M5", "AH7572DRX2", "HJTLFDRX2", "HMVGWDRX2"}
        actual = set(self.meta["Batch"].unique())
        # All actual batches should be in known set
        assert actual.issubset(known), f"Unknown batches: {actual - known}"

    def test_missing_batch_samples(self):
        """Samples in AAAMMK2M5 should exist in metadata."""
        missing = self.meta[self.meta["Batch"] == "AAAMMK2M5"]
        assert len(missing) > 0, "Expected samples in AAAMMK2M5 batch"

    def test_hiv_status_values(self):
        assert set(self.meta["HIVstatus"].unique()).issubset({"positive", "negative"})

    def test_patient_252_hiv_inconsistency(self):
        """Patient 008_252 has one mislabeled VDJ row (positive instead of negative)."""
        p252 = self.meta[self.meta["patientID"] == "008_252"]
        hiv_vals = set(p252["HIVstatus"].unique())
        # Before correction, should have both positive and negative
        assert hiv_vals == {"positive", "negative"}, (
            f"Expected both labels for 008_252 before correction, got {hiv_vals}"
        )
        # Exactly one row should be positive (the mislabeled VDJ row)
        assert (p252["HIVstatus"] == "positive").sum() == 1

    def test_tissue_types(self):
        expected = {"PBMC", "Single Cell Suspension", "CD4 sorted PBMC", "CD8 sorted PBMC"}
        actual = set(self.meta["tissueType"].unique())
        assert actual.issubset(expected), f"Unexpected tissue types: {actual - expected}"

    def test_sorting_types(self):
        expected = {"Unsorted", "CD4pTCells", "CD8pTCells"}
        actual = set(self.meta["sortingCT"].unique())
        assert actual.issubset(expected), f"Unexpected sorting: {actual - expected}"

    def test_sorted_sample_patients(self):
        """Sorted samples should come from known patients."""
        sorted_samples = self.meta[self.meta["sortingCT"].isin(["CD4pTCells", "CD8pTCells"])]
        patients = sorted_samples["patientID"].unique()
        # Main sorted patients, plus CVEU_0985_CD8_SPLIT2 which is also CD8-sorted
        expected = {"008_216", "008_217", "008_220", "CVEU_0985_CD8_SPLIT2"}
        assert set(patients).issubset(expected)


@skip_if_no_clinical
class TestRealClinicalMetadata:
    @pytest.fixture(autouse=True)
    def _load(self):
        self.clinical = pd.read_csv(CLINICAL_PATH)

    def test_expected_columns(self):
        for col in ["ptid", "best_resp"]:
            assert col in self.clinical.columns

    def test_ptid_range(self):
        """ptid should be numeric."""
        assert pd.api.types.is_numeric_dtype(self.clinical["ptid"])


@skip_if_no_data
@skip_if_no_config
class TestRealMetadataWithConfig:
    def test_load_metadata_applies_override(self):
        from utils import load_config, load_metadata
        cfg = load_config(str(CONFIG_PATH))
        meta = load_metadata(cfg)
        p252 = meta[meta["patientID"] == "008_252"]
        # After override, all rows should be "negative"
        assert (p252["HIVstatus"] == "negative").all()

    def test_gex_vdj_sample_filtering(self):
        from utils import load_config, load_metadata, get_gex_samples, get_vdj_samples
        cfg = load_config(str(CONFIG_PATH))
        meta = load_metadata(cfg)
        gex = get_gex_samples(meta)
        vdj = get_vdj_samples(meta)
        assert len(gex) + len(vdj) == len(meta)

    def test_arm_filtering(self):
        from utils import load_config, load_metadata, get_gex_samples, filter_samples_by_arm
        cfg = load_config(str(CONFIG_PATH))
        meta = load_metadata(cfg)
        gex = get_gex_samples(meta)
        primary = filter_samples_by_arm(gex, "primary", cfg)
        focused = filter_samples_by_arm(gex, "focused", cfg)
        assert (primary["sortingCT"] == "Unsorted").all()
        assert focused["sortingCT"].isin(["CD4pTCells", "CD8pTCells"]).all()
        assert len(primary) + len(focused) <= len(gex)

    def test_clinical_mapping(self):
        from utils import (
            load_config, load_metadata, load_clinical_metadata,
            get_gex_samples, map_clinical_metadata,
        )
        cfg = load_config(str(CONFIG_PATH))
        meta = load_metadata(cfg)
        clinical = load_clinical_metadata(cfg)
        gex = get_gex_samples(meta)
        merged = map_clinical_metadata(gex, clinical)
        assert "best_resp" in merged.columns
        # Patients with 008_XXX format should have some clinical matches
        has_clinical = merged["best_resp"].notna()
        assert has_clinical.sum() > 0
