"""
Tests for pipeline/08_trajectory.py: pseudotime analysis, clonotype trajectories.
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
    "trajectory",
    str(Path(__file__).resolve().parents[1] / "pipeline" / "08_trajectory.py"),
)
traj_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(traj_mod)

analyze_expansion_along_pseudotime = traj_mod.analyze_expansion_along_pseudotime
track_clonotype_trajectories = traj_mod.track_clonotype_trajectories


# ===== analyze_expansion_along_pseudotime =====

class TestAnalyzeExpansionAlongPseudotime:
    def test_returns_none_without_pseudotime(self):
        rng = np.random.default_rng(42)
        n = 100
        adata = ad.AnnData(
            X=csr_matrix(rng.poisson(3, (n, 20)).astype(np.float32)),
            obs=pd.DataFrame({
                "clonal_expansion": rng.choice(["1", "2", "3"], n),
            }, index=[f"c{i}" for i in range(n)]),
            var=pd.DataFrame(index=[f"g{i}" for i in range(20)]),
        )
        assert analyze_expansion_along_pseudotime(adata) is None

    def test_returns_none_without_expansion(self):
        rng = np.random.default_rng(42)
        n = 100
        adata = ad.AnnData(
            X=csr_matrix(rng.poisson(3, (n, 20)).astype(np.float32)),
            obs=pd.DataFrame({
                "dpt_pseudotime": rng.uniform(0, 1, n),
            }, index=[f"c{i}" for i in range(n)]),
            var=pd.DataFrame(index=[f"g{i}" for i in range(20)]),
        )
        assert analyze_expansion_along_pseudotime(adata) is None

    def test_returns_binned_results(self, adata_tcell):
        result = analyze_expansion_along_pseudotime(adata_tcell)
        assert result is not None
        assert isinstance(result, pd.DataFrame)
        assert "pseudotime_bin" in result.columns
        assert "frac_expanded" in result.columns
        assert "n_cells" in result.columns
        # Should have up to 20 bins
        assert len(result) <= 20

    def test_expansion_fractions_valid(self, adata_tcell):
        result = analyze_expansion_along_pseudotime(adata_tcell)
        assert (result["frac_expanded"] >= 0).all()
        assert (result["frac_expanded"] <= 1).all()

    def test_mean_pseudotime_increases(self, adata_tcell):
        result = analyze_expansion_along_pseudotime(adata_tcell)
        if len(result) > 1:
            # Mean pseudotime should generally increase with bin number
            sorted_result = result.sort_values("pseudotime_bin")
            assert sorted_result["mean_pseudotime"].iloc[-1] > sorted_result["mean_pseudotime"].iloc[0]


# ===== track_clonotype_trajectories =====

class TestTrackClonotypeTrajectories:
    def test_returns_none_without_clone_id(self):
        rng = np.random.default_rng(42)
        n = 100
        adata = ad.AnnData(
            X=csr_matrix(rng.poisson(3, (n, 20)).astype(np.float32)),
            obs=pd.DataFrame({
                "dpt_pseudotime": rng.uniform(0, 1, n),
            }, index=[f"c{i}" for i in range(n)]),
            var=pd.DataFrame(index=[f"g{i}" for i in range(20)]),
        )
        assert track_clonotype_trajectories(adata) is None

    def test_returns_none_without_pseudotime(self):
        rng = np.random.default_rng(42)
        n = 100
        adata = ad.AnnData(
            X=csr_matrix(rng.poisson(3, (n, 20)).astype(np.float32)),
            obs=pd.DataFrame({
                "clone_id": rng.choice(["c1", "c2", "c3"], n),
            }, index=[f"c{i}" for i in range(n)]),
            var=pd.DataFrame(index=[f"g{i}" for i in range(20)]),
        )
        assert track_clonotype_trajectories(adata) is None

    def test_tracks_expanded_clones(self, adata_tcell):
        result = track_clonotype_trajectories(adata_tcell, top_n=10)
        assert result is not None
        assert isinstance(result, pd.DataFrame)
        assert "clone_id" in result.columns
        assert "n_cells" in result.columns
        assert "mean_pseudotime" in result.columns
        assert "pseudotime_span" in result.columns
        # All tracked clones should have n_cells > 1
        assert (result["n_cells"] > 1).all()

    def test_top_n_limits_output(self, adata_tcell):
        result5 = track_clonotype_trajectories(adata_tcell, top_n=5)
        result20 = track_clonotype_trajectories(adata_tcell, top_n=20)
        assert len(result5) <= 5
        assert len(result20) <= 20

    def test_pseudotime_span_nonnegative(self, adata_tcell):
        result = track_clonotype_trajectories(adata_tcell)
        if result is not None:
            assert (result["pseudotime_span"] >= 0).all()

    def test_dominant_subtype_recorded(self, adata_tcell):
        result = track_clonotype_trajectories(adata_tcell)
        if result is not None:
            assert "dominant_subtype" in result.columns
            assert "n_subtypes" in result.columns
            assert (result["n_subtypes"] >= 1).all()
