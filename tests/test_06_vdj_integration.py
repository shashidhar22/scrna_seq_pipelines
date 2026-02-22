"""
Tests for pipeline/06_vdj_integration.py: multi-chain doublet flagging.
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
    "vdj_integration",
    str(Path(__file__).resolve().parents[1] / "pipeline" / "06_vdj_integration.py"),
)
vdj_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(vdj_mod)

flag_multi_chain_doublets = vdj_mod.flag_multi_chain_doublets


# ===== flag_multi_chain_doublets =====

class TestFlagMultiChainDoublets:
    def _make_vdj_adata(self, n_obs=100, dual_frac=0.1, seed=42):
        rng = np.random.default_rng(seed)
        adata = ad.AnnData(
            X=csr_matrix(rng.poisson(3, (n_obs, 20)).astype(np.float32)),
            obs=pd.DataFrame({
                "IR_VDJ_1_junction_aa": [
                    "CASSLAPGATNEKLFF" if rng.random() > 0.3 else None
                    for _ in range(n_obs)
                ],
                "IR_VDJ_2_junction_aa": [
                    "CASSLGGQETQYF" if rng.random() < dual_frac else None
                    for _ in range(n_obs)
                ],
                "doublet_score": rng.uniform(0, 0.5, n_obs),
            }, index=[f"c{i}" for i in range(n_obs)]),
            var=pd.DataFrame(index=[f"g{i}" for i in range(20)]),
        )
        return adata

    def test_adds_multi_trb_column(self):
        adata = self._make_vdj_adata()
        result = flag_multi_chain_doublets(adata)
        assert "multi_trb" in result.obs.columns
        assert result.obs["multi_trb"].dtype == bool

    def test_detects_dual_trb(self):
        """Cells with IR_VDJ_2_junction_aa should be flagged."""
        adata = self._make_vdj_adata(dual_frac=0.5)
        result = flag_multi_chain_doublets(adata)
        # At least some should be flagged
        assert result.obs["multi_trb"].sum() > 0

    def test_no_dual_trb(self):
        """When no cell has dual TRB, no flags."""
        adata = self._make_vdj_adata(dual_frac=0.0)
        result = flag_multi_chain_doublets(adata)
        assert result.obs["multi_trb"].sum() == 0

    def test_works_without_doublet_score(self):
        """Should still work if doublet_score column is absent."""
        adata = self._make_vdj_adata()
        del adata.obs["doublet_score"]
        result = flag_multi_chain_doublets(adata)
        assert "multi_trb" in result.obs.columns
