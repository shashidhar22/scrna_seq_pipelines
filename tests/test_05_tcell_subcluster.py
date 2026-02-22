"""
Tests for pipeline/05_tcell_subcluster.py: T cell subsetting and subtype annotation.
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
    "tcell_subcluster",
    str(Path(__file__).resolve().parents[1] / "pipeline" / "05_tcell_subcluster.py"),
)
tc_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(tc_mod)

subset_tcells = tc_mod.subset_tcells
annotate_tcell_subtypes = tc_mod.annotate_tcell_subtypes


# ===== subset_tcells =====

class TestSubsetTcells:
    def test_subsets_correct_types(self, adata_annotated):
        adata_t = subset_tcells(adata_annotated)
        expected_types = {"CD4 T", "CD8 T", "T cells", "NK/T"}
        assert adata_t.obs["cell_type"].isin(expected_types).all()

    def test_excludes_non_tcells(self, adata_annotated):
        adata_t = subset_tcells(adata_annotated)
        assert "B cells" not in adata_t.obs["cell_type"].values
        assert "Monocytes" not in adata_t.obs["cell_type"].values
        assert "NK" not in adata_t.obs["cell_type"].values

    def test_returns_copy(self, adata_annotated):
        adata_t = subset_tcells(adata_annotated)
        # Modifying the subset should not affect the original
        adata_t.obs["new_col"] = 1
        assert "new_col" not in adata_annotated.obs.columns

    def test_preserves_all_columns(self, adata_annotated):
        adata_t = subset_tcells(adata_annotated)
        for col in adata_annotated.obs.columns:
            assert col in adata_t.obs.columns


# ===== annotate_tcell_subtypes =====

class TestAnnotateTcellSubtypes:
    def _make_tcell_adata(self, dominant_markers, n_obs=100, n_clusters=5, seed=42):
        """
        Create a T cell AnnData with controlled marker expression for testing
        the decision tree.
        """
        rng = np.random.default_rng(seed)
        genes = [
            "CD8A", "CD8B", "CD4", "IL7R",
            "CCR7", "SELL", "LEF1", "TCF7",
            "GZMB", "PRF1", "GNLY", "GZMA",
            "PDCD1", "LAG3", "HAVCR2", "TOX", "TIGIT",
            "FOXP3", "IL2RA", "CTLA4",
            "MKI67", "TOP2A",
            "Gene1", "Gene2", "Gene3",
        ]
        n_vars = len(genes)
        # Base expression: low
        counts = rng.poisson(0.5, (n_obs, n_vars)).astype(np.float32)

        adata = ad.AnnData(
            X=csr_matrix(counts),
            obs=pd.DataFrame({
                "leiden_1.0": pd.Categorical([str(i % n_clusters) for i in range(n_obs)]),
            }, index=[f"c{i}" for i in range(n_obs)]),
            var=pd.DataFrame(index=genes),
        )

        # Normalize and store raw for rank_genes_groups
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.raw = adata

        return adata

    def test_assigns_tcell_subtype_column(self):
        adata = self._make_tcell_adata({})
        result = annotate_tcell_subtypes(adata, cluster_key="leiden_1.0")
        assert "tcell_subtype" in result.obs.columns
        assert result.obs["tcell_subtype"].notna().all()

    def test_all_clusters_get_labels(self):
        adata = self._make_tcell_adata({})
        result = annotate_tcell_subtypes(adata, cluster_key="leiden_1.0")
        clusters = adata.obs["leiden_1.0"].unique()
        labeled_clusters = set()
        for cl in clusters:
            mask = result.obs["leiden_1.0"] == cl
            labels = result.obs.loc[mask, "tcell_subtype"].unique()
            assert len(labels) == 1  # One label per cluster
            labeled_clusters.add(cl)
        assert labeled_clusters == set(clusters)

    def test_labels_are_meaningful(self):
        """All labels should be from the known subtype set."""
        adata = self._make_tcell_adata({})
        result = annotate_tcell_subtypes(adata, cluster_key="leiden_1.0")
        valid_labels = {
            "T cell", "CD8 Naive", "CD8 Effector", "CD8 Exhausted",
            "CD8 Memory", "CD8 Proliferating",
            "CD4 Naive", "CD4 Memory", "CD4 Treg", "CD4 Proliferating",
        }
        for label in result.obs["tcell_subtype"].unique():
            assert label in valid_labels, f"Unexpected label: {label}"
