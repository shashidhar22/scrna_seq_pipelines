"""
Tests for pipeline/04_cluster_annotate.py: coarse label assignment, marker genes.
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
    "cluster_annotate",
    str(Path(__file__).resolve().parents[1] / "pipeline" / "04_cluster_annotate.py"),
)
ann_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(ann_mod)

CANONICAL_MARKERS = ann_mod.CANONICAL_MARKERS
assign_coarse_labels = ann_mod.assign_coarse_labels
find_marker_genes = ann_mod.find_marker_genes


# ===== CANONICAL_MARKERS =====

class TestCanonicalMarkers:
    def test_contains_expected_cell_types(self):
        expected = ["T cells", "CD4 T", "CD8 T", "NK", "B cells", "Monocytes CD14+"]
        for ct in expected:
            assert ct in CANONICAL_MARKERS, f"Missing cell type: {ct}"

    def test_markers_are_lists_of_strings(self):
        for ct, genes in CANONICAL_MARKERS.items():
            assert isinstance(genes, list), f"{ct}: markers should be a list"
            for g in genes:
                assert isinstance(g, str), f"{ct}: marker {g} should be a string"


# ===== assign_coarse_labels =====

class TestAssignCoarseLabels:
    def _make_annotated_adata(self, celltypist_labels):
        """Create AnnData with CellTypist predictions."""
        n = len(celltypist_labels)
        rng = np.random.default_rng(42)
        adata = ad.AnnData(
            X=csr_matrix(rng.poisson(3, (n, 30)).astype(np.float32)),
            obs=pd.DataFrame({
                "celltypist_Immune_All_Low_majority": pd.Categorical(celltypist_labels),
                "leiden_0.8": pd.Categorical([str(i % 5) for i in range(n)]),
            }, index=[f"c{i}" for i in range(n)]),
            var=pd.DataFrame(index=[f"g{i}" for i in range(30)]),
        )
        return adata

    def test_t_cell_mapping(self):
        labels = ["CD4+ Tcm", "CD8+ Tem", "Treg", "T helper cells"]
        adata = self._make_annotated_adata(labels)
        result = assign_coarse_labels(adata)
        assert result.obs.loc["c0", "cell_type"] == "CD4 T"
        assert result.obs.loc["c1", "cell_type"] == "CD8 T"
        assert result.obs.loc["c2", "cell_type"] == "CD4 T"  # Treg -> CD4 T
        assert result.obs.loc["c3", "cell_type"] == "T cells"

    def test_nk_mapping(self):
        labels = ["NK cells", "NKT cells"]
        adata = self._make_annotated_adata(labels)
        result = assign_coarse_labels(adata)
        assert result.obs.loc["c0", "cell_type"] == "NK"
        assert result.obs.loc["c1", "cell_type"] == "NK/T"

    def test_b_cell_mapping(self):
        labels = ["B cells", "Plasma cells"]
        adata = self._make_annotated_adata(labels)
        result = assign_coarse_labels(adata)
        assert result.obs.loc["c0", "cell_type"] == "B cells"
        assert result.obs.loc["c1", "cell_type"] == "Plasma"

    def test_myeloid_mapping(self):
        labels = ["Classical Monocytes", "DC2", "pDC precursors"]
        adata = self._make_annotated_adata(labels)
        result = assign_coarse_labels(adata)
        assert result.obs.loc["c0", "cell_type"] == "Monocytes"
        assert result.obs.loc["c1", "cell_type"] == "DCs"
        assert result.obs.loc["c2", "cell_type"] == "pDCs"

    def test_unknown_label(self):
        """A label with no matching key should map to Unknown."""
        labels = ["Granulocyes"]  # No substring match for any key in coarse_map
        adata = self._make_annotated_adata(labels)
        result = assign_coarse_labels(adata)
        assert result.obs.loc["c0", "cell_type"] == "Unknown"

    def test_cell_subtype_stores_fine_label(self):
        labels = ["CD4+ Tcm", "Classical Monocytes"]
        adata = self._make_annotated_adata(labels)
        result = assign_coarse_labels(adata)
        assert result.obs.loc["c0", "cell_subtype"] == "CD4+ Tcm"
        assert result.obs.loc["c1", "cell_subtype"] == "Classical Monocytes"

    def test_fallback_when_no_celltypist(self):
        """Without CellTypist columns, should use cluster IDs."""
        rng = np.random.default_rng(42)
        n = 50
        adata = ad.AnnData(
            X=csr_matrix(rng.poisson(3, (n, 20)).astype(np.float32)),
            obs=pd.DataFrame({
                "leiden_0.8": pd.Categorical([str(i % 5) for i in range(n)]),
            }, index=[f"c{i}" for i in range(n)]),
            var=pd.DataFrame(index=[f"g{i}" for i in range(20)]),
        )
        result = assign_coarse_labels(adata, cluster_key="leiden_0.8")
        # cell_type should be cluster IDs
        assert "cell_type" in result.obs.columns


# ===== find_marker_genes =====

class TestFindMarkerGenes:
    def test_returns_dataframe(self, adata_annotated):
        """Should return a DataFrame with marker genes per cluster."""
        adata = adata_annotated.copy()
        # Need normalized data in .raw for rank_genes_groups
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        adata.raw = adata

        result = find_marker_genes(adata, cluster_key="leiden_0.8", n_genes=5)
        assert isinstance(result, pd.DataFrame)
        assert "cluster" in result.columns
        assert "gene" in result.columns
        assert "logfoldchange" in result.columns
        assert "pval_adj" in result.columns
        # Should have 5 genes per cluster
        for cl in adata.obs["leiden_0.8"].unique():
            assert len(result[result["cluster"] == cl]) == 5
