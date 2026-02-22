"""
Tests for pipeline/03_merge_integrate.py: merge, normalization, HVG, clustering.
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
    "merge_integrate",
    str(Path(__file__).resolve().parents[1] / "pipeline" / "03_merge_integrate.py"),
)
merge_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(merge_mod)

merge_libraries = merge_mod.merge_libraries
normalize_and_hvg = merge_mod.normalize_and_hvg
cluster_and_umap = merge_mod.cluster_and_umap


# ===== Helpers =====

def _make_library(name, n_obs=100, n_vars=80, batch="B1", patient="008_216",
                   hiv="negative", seed=42):
    rng = np.random.default_rng(seed)
    counts = rng.poisson(5, (n_obs, n_vars)).astype(np.float32)
    gene_names = [f"Gene{i}" for i in range(n_vars)]
    adata = ad.AnnData(
        X=csr_matrix(counts),
        obs=pd.DataFrame({"sample_id": name}, index=[f"{name}_{i}" for i in range(n_obs)]),
        var=pd.DataFrame(index=gene_names),
    )
    return adata


def _sample_meta():
    """Metadata matching the libraries."""
    return pd.DataFrame([
        {"sampleName": "S1", "patientID": "008_216", "Batch": "B1",
         "sortingCT": "Unsorted", "tissueType": "PBMC", "HIVstatus": "negative",
         "locus": "5primeGEX"},
        {"sampleName": "S2", "patientID": "008_220", "Batch": "B2",
         "sortingCT": "Unsorted", "tissueType": "PBMC", "HIVstatus": "positive",
         "locus": "5primeGEX"},
    ])


# ===== merge_libraries =====

class TestMergeLibraries:
    def test_inner_join_on_genes(self):
        """Merge should keep only genes present in all libraries."""
        adata1 = _make_library("S1", n_obs=50, n_vars=80, seed=1)
        adata2 = _make_library("S2", n_obs=60, n_vars=80, seed=2)
        # Remove some genes from adata2
        adata2 = adata2[:, :70].copy()

        meta = _sample_meta()
        merged = merge_libraries([adata1, adata2], meta)
        assert merged.n_obs == 110
        assert merged.n_vars == 70  # inner join

    def test_adds_metadata_columns(self):
        adata1 = _make_library("S1", seed=1)
        adata2 = _make_library("S2", seed=2)
        meta = _sample_meta()
        merged = merge_libraries([adata1, adata2], meta)

        for col in ["patientID", "Batch", "sortingCT", "tissueType", "HIVstatus"]:
            assert col in merged.obs.columns

    def test_stores_counts_layer(self):
        adata1 = _make_library("S1", seed=1)
        meta = _sample_meta()
        merged = merge_libraries([adata1], meta)
        assert "counts" in merged.layers

    def test_clinical_metadata_merge(self):
        adata1 = _make_library("S1", seed=1)
        meta = _sample_meta()
        clinical = pd.DataFrame({
            "ptid": [216, 220],
            "best_resp": ["CR", "PR"],
            "bl_plasma_grp": ["low", "high"],
            "best_resp_new": ["CR", "PR"],
        })
        merged = merge_libraries([adata1], meta, clinical)
        assert "best_resp" in merged.obs.columns


# ===== normalize_and_hvg =====

class TestNormalizeAndHvg:
    def test_selects_hvgs(self, minimal_config):
        """Should select highly variable genes."""
        rng = np.random.default_rng(42)
        n_obs, n_vars = 300, 500
        # Generate count data with varying dispersion across genes
        counts = rng.negative_binomial(5, 0.3, (n_obs, n_vars)).astype(np.float32)
        # Make some genes very noisy (highly variable)
        counts[:, :50] = rng.negative_binomial(2, 0.1, (n_obs, 50)).astype(np.float32)

        gene_names = [f"Gene{i}" for i in range(n_vars)]
        gene_names[0] = "MT-CO1"
        gene_names[1] = "RPS1"
        gene_names[2] = "HBA1"

        adata = ad.AnnData(
            X=csr_matrix(counts),
            obs=pd.DataFrame({
                "sample_id": rng.choice(["S1", "S2"], n_obs),
            }, index=[f"c{i}" for i in range(n_obs)]),
            var=pd.DataFrame(index=gene_names),
        )
        adata.layers["counts"] = adata.X.copy()

        # Mark QC gene types
        adata.var["mt"] = adata.var_names.str.startswith("MT-")
        adata.var["ribo"] = adata.var_names.str.match(r"^RP[SL]\d")
        adata.var["hb"] = adata.var_names.str.match(r"^HB[^P]")

        cfg = minimal_config.copy()
        cfg["integration"]["n_top_genes"] = 100

        result = normalize_and_hvg(adata, cfg)
        assert "highly_variable" in result.var.columns
        assert result.var["highly_variable"].sum() > 0
        # MT/ribo/HB genes should not be HVGs
        assert not result.var.loc["MT-CO1", "highly_variable"]
        assert not result.var.loc["HBA1", "highly_variable"]

    def test_stores_raw(self, minimal_config):
        rng = np.random.default_rng(42)
        n_obs, n_vars = 200, 300
        counts = rng.negative_binomial(5, 0.3, (n_obs, n_vars)).astype(np.float32)
        adata = ad.AnnData(
            X=csr_matrix(counts),
            obs=pd.DataFrame({"sample_id": rng.choice(["S1", "S2"], n_obs)},
                             index=[f"c{i}" for i in range(n_obs)]),
            var=pd.DataFrame(index=[f"G{i}" for i in range(n_vars)]),
        )
        adata.layers["counts"] = adata.X.copy()

        cfg = minimal_config.copy()
        cfg["integration"]["n_top_genes"] = 50

        result = normalize_and_hvg(adata, cfg)
        assert result.raw is not None


# ===== cluster_and_umap =====

class TestClusterAndUmap:
    def test_multi_resolution_leiden(self, minimal_config):
        """Should create Leiden clusters at multiple resolutions."""
        rng = np.random.default_rng(42)
        n_obs = 300
        n_latent = 10
        # Create two distinct clusters in latent space
        latent = np.vstack([
            rng.normal(loc=-3, size=(150, n_latent)),
            rng.normal(loc=3, size=(150, n_latent)),
        ]).astype(np.float32)

        adata = ad.AnnData(
            X=csr_matrix(rng.poisson(3, (n_obs, 50)).astype(np.float32)),
            obs=pd.DataFrame(index=[f"c{i}" for i in range(n_obs)]),
            var=pd.DataFrame(index=[f"g{i}" for i in range(50)]),
        )
        adata.obsm["X_scVI"] = latent

        result = cluster_and_umap(adata, minimal_config, use_rep="X_scVI")
        for res in minimal_config["integration"]["leiden_resolutions"]:
            key = f"leiden_{res}"
            assert key in result.obs.columns
        # At least the highest resolution should produce >1 cluster
        assert result.obs["leiden_1.5"].nunique() > 1

    def test_umap_computed(self, minimal_config):
        rng = np.random.default_rng(42)
        n_obs = 200
        adata = ad.AnnData(
            X=csr_matrix(rng.poisson(3, (n_obs, 50)).astype(np.float32)),
            obs=pd.DataFrame(index=[f"c{i}" for i in range(n_obs)]),
            var=pd.DataFrame(index=[f"g{i}" for i in range(50)]),
        )
        adata.obsm["X_scVI"] = rng.normal(size=(n_obs, 10)).astype(np.float32)

        result = cluster_and_umap(adata, minimal_config, use_rep="X_scVI")
        assert "X_umap" in result.obsm
        assert result.obsm["X_umap"].shape == (n_obs, 2)
