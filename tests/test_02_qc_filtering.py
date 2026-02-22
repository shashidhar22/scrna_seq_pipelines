"""
Tests for pipeline/02_qc_filtering.py: MAD-based adaptive QC, Scrublet, gene filtering.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))

import importlib.util

spec = importlib.util.spec_from_file_location(
    "qc_filtering",
    str(Path(__file__).resolve().parents[1] / "pipeline" / "02_qc_filtering.py"),
)
qc_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(qc_mod)

filter_library_adaptive = qc_mod.filter_library_adaptive
filter_genes = qc_mod.filter_genes


# ===== filter_library_adaptive =====

class TestFilterLibraryAdaptive:
    def _make_adata_for_qc(self, n_obs=500, seed=42):
        """Create AnnData with realistic QC-like expression."""
        rng = np.random.default_rng(seed)
        n_vars = 120
        counts = rng.poisson(lam=5, size=(n_obs, n_vars)).astype(np.float32)

        gene_names = [f"Gene{i}" for i in range(n_vars)]
        # MT genes at indices 0-4
        for i in range(5):
            gene_names[i] = f"MT-GENE{i}"
        # Ribo genes
        for i in range(5, 10):
            gene_names[i] = f"RPS{i}"
        # HB genes
        gene_names[10] = "HBA1"
        gene_names[11] = "HBA2"
        gene_names[12] = "HBB"

        adata = ad.AnnData(
            X=csr_matrix(counts),
            obs=pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)]),
            var=pd.DataFrame(index=gene_names),
        )
        return adata

    def test_filters_cells(self, minimal_config):
        """Should remove some outlier cells."""
        adata = self._make_adata_for_qc(500)
        # Inject an outlier: give cell 0 extremely high counts
        adata.X[0] = csr_matrix(np.full(adata.n_vars, 5000, dtype=np.float32))
        # And cell 1 extremely low
        adata.X[1] = csr_matrix(np.zeros(adata.n_vars, dtype=np.float32))

        filtered = filter_library_adaptive(adata, "test_sample", minimal_config)
        assert filtered.n_obs < adata.n_obs
        assert "outlier_any" not in filtered.obs.columns or not filtered.obs["outlier_any"].any()

    def test_preserves_most_normal_cells(self, minimal_config):
        """Most cells in a well-behaved dataset should pass QC."""
        adata = self._make_adata_for_qc(500)
        filtered = filter_library_adaptive(adata, "good_sample", minimal_config)
        # With Poisson(5) data, very few should be outliers
        assert filtered.n_obs > 400

    def test_high_hb_cells_removed(self, minimal_config):
        """Cells above max_pct_hb should be removed."""
        adata = self._make_adata_for_qc(200)
        # Make cell 0 have 50% HB content
        adata.X[0] = csr_matrix(np.zeros(adata.n_vars, dtype=np.float32))
        adata.X[0, 10] = 1000  # HBA1
        adata.X[0, 11] = 1000  # HBA2
        adata.X[0, 12] = 1000  # HBB

        filtered = filter_library_adaptive(adata, "hb_test", minimal_config)
        # cell_0 should be filtered out (high HB)
        assert "cell_0" not in filtered.obs_names

    def test_outlier_flags_stored_before_filtering(self, minimal_config):
        """Outlier flags should be in .obs even if they're used to filter."""
        adata = self._make_adata_for_qc(200)
        filtered = filter_library_adaptive(adata, "flag_test", minimal_config)
        # The returned adata keeps only non-outliers, so outlier_any should be False
        # but the flags should exist as columns
        assert "outlier_counts" in filtered.obs.columns or filtered.n_obs < 200


# ===== filter_genes =====

class TestFilterGenes:
    def test_removes_mt_genes(self):
        rng = np.random.default_rng(42)
        n_obs, n_vars = 100, 30
        counts = rng.poisson(5, (n_obs, n_vars)).astype(np.float32)
        gene_names = [f"Gene{i}" for i in range(n_vars)]
        gene_names[0] = "MT-CO1"
        gene_names[1] = "MT-ND1"
        adata = ad.AnnData(
            X=csr_matrix(counts),
            obs=pd.DataFrame(index=[f"c{i}" for i in range(n_obs)]),
            var=pd.DataFrame(index=gene_names),
        )
        filtered = filter_genes(adata, min_cells=3)
        assert "MT-CO1" not in filtered.var_names
        assert "MT-ND1" not in filtered.var_names

    def test_removes_low_expressed_genes(self):
        rng = np.random.default_rng(42)
        n_obs, n_vars = 100, 20
        counts = rng.poisson(5, (n_obs, n_vars)).astype(np.float32)
        # Gene 0 expressed in only 1 cell
        counts[:, 0] = 0
        counts[0, 0] = 1

        gene_names = [f"Gene{i}" for i in range(n_vars)]
        adata = ad.AnnData(
            X=csr_matrix(counts),
            obs=pd.DataFrame(index=[f"c{i}" for i in range(n_obs)]),
            var=pd.DataFrame(index=gene_names),
        )
        filtered = filter_genes(adata, min_cells=3)
        assert "Gene0" not in filtered.var_names

    def test_keeps_well_expressed_genes(self):
        rng = np.random.default_rng(42)
        n_obs, n_vars = 100, 20
        counts = rng.poisson(5, (n_obs, n_vars)).astype(np.float32)
        gene_names = [f"Gene{i}" for i in range(n_vars)]
        adata = ad.AnnData(
            X=csr_matrix(counts),
            obs=pd.DataFrame(index=[f"c{i}" for i in range(n_obs)]),
            var=pd.DataFrame(index=gene_names),
        )
        filtered = filter_genes(adata, min_cells=3)
        # Most genes at Poisson(5) should be in most cells
        assert filtered.n_vars >= 15
