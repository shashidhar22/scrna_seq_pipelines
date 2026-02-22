"""
Tests for pipeline/07_differential_expression.py: pseudobulk, cell proportions.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import anndata as ad
from scipy.sparse import csr_matrix

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))

import importlib.util

spec = importlib.util.spec_from_file_location(
    "differential_expression",
    str(Path(__file__).resolve().parents[1] / "pipeline" / "07_differential_expression.py"),
)
de_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(de_mod)

compute_pseudobulk = de_mod.compute_pseudobulk
analyze_cell_proportions = de_mod.analyze_cell_proportions


# ===== compute_pseudobulk =====

class TestComputePseudobulk:
    def _make_adata_for_pb(self, n_patients=4, n_cells_per=50, seed=42):
        rng = np.random.default_rng(seed)
        patients = [f"P{i}" for i in range(n_patients)]
        cell_types = ["CD4 T", "CD8 T", "B cells"]
        n_obs = n_patients * n_cells_per
        n_vars = 100

        counts = rng.poisson(10, (n_obs, n_vars)).astype(np.float32)
        obs_data = []
        for p in patients:
            for _ in range(n_cells_per):
                obs_data.append({
                    "patientID": p,
                    "cell_type": rng.choice(cell_types),
                    "HIVstatus": "positive" if int(p[1]) % 2 == 0 else "negative",
                    "Batch": "B1",
                    "tissueType": "PBMC",
                })
        obs = pd.DataFrame(obs_data, index=[f"c{i}" for i in range(n_obs)])

        adata = ad.AnnData(
            X=csr_matrix(counts),
            obs=obs,
            var=pd.DataFrame(index=[f"Gene{i}" for i in range(n_vars)]),
        )
        adata.layers["counts"] = adata.X.copy()
        return adata

    def test_returns_counts_and_meta(self):
        adata = self._make_adata_for_pb()
        result = compute_pseudobulk(
            adata, sample_col="patientID", groups_col="cell_type",
            min_cells=5, min_counts=50,
        )
        # Manual aggregation path returns (counts_df, meta_df)
        assert isinstance(result, tuple)
        counts_df, meta_df = result
        assert isinstance(counts_df, pd.DataFrame)
        assert isinstance(meta_df, pd.DataFrame)
        assert counts_df.shape[1] == adata.n_vars

    def test_filters_low_cell_count(self):
        """Groups below min_cells threshold should be excluded."""
        adata = self._make_adata_for_pb(n_patients=4, n_cells_per=100)
        # With 100 cells per patient across 3 cell types, ~33 cells per group
        counts_df, meta_df = compute_pseudobulk(
            adata, sample_col="patientID", groups_col="cell_type",
            min_cells=20, min_counts=10,
        )
        # All retained groups should have >= 20 cells
        assert len(meta_df) > 0
        assert all(meta_df["n_cells"] >= 20)

    def test_pseudobulk_sums_correctly(self):
        """Pseudobulk counts should equal sum of single-cell counts."""
        rng = np.random.default_rng(42)
        n_obs, n_vars = 20, 10
        counts = rng.poisson(5, (n_obs, n_vars)).astype(np.float32)
        adata = ad.AnnData(
            X=csr_matrix(counts),
            obs=pd.DataFrame({
                "patientID": ["P1"] * 20,
                "cell_type": ["T"] * 20,
                "HIVstatus": ["pos"] * 20,
                "Batch": ["B1"] * 20,
                "tissueType": ["PBMC"] * 20,
            }, index=[f"c{i}" for i in range(n_obs)]),
            var=pd.DataFrame(index=[f"G{i}" for i in range(n_vars)]),
        )
        adata.layers["counts"] = adata.X.copy()
        counts_df, meta_df = compute_pseudobulk(
            adata, sample_col="patientID", groups_col="cell_type",
            min_cells=1, min_counts=1,
        )
        # Should have one row: P1__T
        assert len(counts_df) == 1
        # Sum should match
        expected = np.asarray(adata.X.sum(axis=0)).flatten()
        actual = counts_df.iloc[0].values
        np.testing.assert_array_almost_equal(actual, expected, decimal=0)


# ===== analyze_cell_proportions =====

class TestAnalyzeCellProportions:
    def _make_adata_for_props(self, seed=42):
        rng = np.random.default_rng(seed)
        n_obs = 400
        patients = rng.choice(["P1", "P2", "P3", "P4"], n_obs)
        cell_types = rng.choice(["CD4 T", "CD8 T", "B cells"], n_obs)
        hiv = pd.Series(patients).map({
            "P1": "positive", "P2": "positive",
            "P3": "negative", "P4": "negative",
        }).values

        adata = ad.AnnData(
            X=csr_matrix(rng.poisson(3, (n_obs, 20)).astype(np.float32)),
            obs=pd.DataFrame({
                "patientID": patients,
                "cell_type": cell_types,
                "HIVstatus": hiv,
            }, index=[f"c{i}" for i in range(n_obs)]),
            var=pd.DataFrame(index=[f"g{i}" for i in range(20)]),
        )
        return adata

    def test_returns_dataframe(self):
        adata = self._make_adata_for_props()
        result = analyze_cell_proportions(adata)
        assert isinstance(result, pd.DataFrame)
        assert "cell_type" in result.columns
        assert "pval" in result.columns
        assert "mean_positive" in result.columns
        assert "mean_negative" in result.columns

    def test_all_cell_types_present(self):
        adata = self._make_adata_for_props()
        result = analyze_cell_proportions(adata)
        assert set(result["cell_type"]) == {"CD4 T", "CD8 T", "B cells"}

    def test_proportions_sum_near_one(self):
        """Mean proportions for each group should approximately sum to 1."""
        adata = self._make_adata_for_props()
        result = analyze_cell_proportions(adata)
        total_pos = result["mean_positive"].sum()
        total_neg = result["mean_negative"].sum()
        assert abs(total_pos - 1.0) < 0.1
        assert abs(total_neg - 1.0) < 0.1

    def test_fdr_correction_applied(self):
        adata = self._make_adata_for_props()
        result = analyze_cell_proportions(adata)
        if "padj" in result.columns:
            # padj should be >= pval
            valid = result.dropna(subset=["pval", "padj"])
            assert (valid["padj"] >= valid["pval"]).all()
