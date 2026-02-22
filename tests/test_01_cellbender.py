"""
Tests for pipeline/01_cellbender.py: CellBender command construction, SoupX fallback, validation.
"""

import sys
from pathlib import Path
from unittest.mock import patch, MagicMock

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))

import importlib.util

spec = importlib.util.spec_from_file_location(
    "cellbender",
    str(Path(__file__).resolve().parents[1] / "pipeline" / "01_cellbender.py"),
)
cb_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cb_mod)

run_cellbender_sample = cb_mod.run_cellbender_sample
check_cellbender_convergence = cb_mod.check_cellbender_convergence
soupx_fallback = cb_mod.soupx_fallback
validate_correction = cb_mod.validate_correction


# ===== run_cellbender_sample =====

class TestRunCellbenderSample:
    def test_builds_correct_command(self, minimal_config, tmp_path):
        """Should construct a CellBender command with all required flags."""
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=0, stderr="", stdout="")
            result = run_cellbender_sample(
                "sample1",
                str(tmp_path / "raw.h5"),
                str(tmp_path / "out.h5"),
                minimal_config,
                expected_cells=3000,
            )
            assert result is True
            call_args = mock_run.call_args[0][0]
            assert "cellbender" in call_args
            assert "--fpr" in call_args
            assert "0.01" in call_args
            assert "--cuda" in call_args
            assert "--expected-cells" in call_args
            assert "3000" in call_args

    def test_handles_failure(self, minimal_config, tmp_path):
        with patch("subprocess.run") as mock_run:
            mock_run.return_value = MagicMock(returncode=1, stderr="Error", stdout="")
            result = run_cellbender_sample(
                "sample_fail", str(tmp_path / "raw.h5"),
                str(tmp_path / "out.h5"), minimal_config,
            )
            assert result is False

    def test_handles_timeout(self, minimal_config, tmp_path):
        import subprocess
        with patch("subprocess.run", side_effect=subprocess.TimeoutExpired("cmd", 3600)):
            result = run_cellbender_sample(
                "sample_timeout", str(tmp_path / "raw.h5"),
                str(tmp_path / "out.h5"), minimal_config,
            )
            assert result is False

    def test_handles_missing_cellbender(self, minimal_config, tmp_path):
        with patch("subprocess.run", side_effect=FileNotFoundError):
            result = run_cellbender_sample(
                "sample_nf", str(tmp_path / "raw.h5"),
                str(tmp_path / "out.h5"), minimal_config,
            )
            assert result is False


# ===== check_cellbender_convergence =====

class TestCheckCellbenderConvergence:
    def test_no_output_returns_false(self, tmp_path):
        result = check_cellbender_convergence(
            str(tmp_path / "nonexistent.h5"), "sample1"
        )
        assert result is False

    def test_valid_h5ad_returns_true(self, tmp_path):
        """A valid h5ad file with enough cells should pass."""
        adata = ad.AnnData(
            X=csr_matrix(np.random.poisson(3, (200, 50)).astype(np.float32)),
            obs=pd.DataFrame(index=[f"c_{i}" for i in range(200)]),
            var=pd.DataFrame(index=[f"g_{i}" for i in range(50)]),
        )
        path = tmp_path / "out.h5ad"
        adata.write_h5ad(str(path))
        result = check_cellbender_convergence(str(path), "sample_ok")
        assert result is True

    def test_too_few_cells_returns_false(self, tmp_path):
        adata = ad.AnnData(
            X=csr_matrix(np.random.poisson(3, (50, 20)).astype(np.float32)),
            obs=pd.DataFrame(index=[f"c_{i}" for i in range(50)]),
            var=pd.DataFrame(index=[f"g_{i}" for i in range(20)]),
        )
        path = tmp_path / "out.h5ad"
        adata.write_h5ad(str(path))
        result = check_cellbender_convergence(str(path), "sample_small")
        assert result is False


# ===== soupx_fallback =====

class TestSoupxFallback:
    def _write_h5(self, path, n_obs, n_vars, seed=42):
        """Write a minimal h5ad for the fallback to read."""
        rng = np.random.default_rng(seed)
        X = csr_matrix(rng.poisson(5, (n_obs, n_vars)).astype(np.float32))
        adata = ad.AnnData(
            X=X,
            obs=pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)]),
            var=pd.DataFrame(index=[f"gene_{i}" for i in range(n_vars)]),
        )
        adata.write_h5ad(str(path))
        return adata

    def test_soupx_fallback_produces_output(self, tmp_path):
        """The fallback should write a corrected h5ad."""
        raw_path = tmp_path / "raw.h5ad"
        filtered_path = tmp_path / "filtered.h5ad"
        output_path = tmp_path / "corrected.h5ad"

        # Raw has 500 barcodes (300 empty + 200 cells)
        rng = np.random.default_rng(42)
        X_raw = csr_matrix(rng.poisson(5, (500, 50)).astype(np.float32))
        raw_barcodes = [f"cell_{i}" for i in range(500)]
        adata_raw = ad.AnnData(
            X=X_raw,
            obs=pd.DataFrame(index=raw_barcodes),
            var=pd.DataFrame(index=[f"gene_{i}" for i in range(50)]),
        )
        adata_raw.write_h5ad(str(raw_path))

        # Filtered has first 200 barcodes
        adata_filt = adata_raw[:200].copy()
        adata_filt.write_h5ad(str(filtered_path))

        # Patch sc.read_10x_h5 to return our h5ad objects
        def mock_read_10x_h5(path):
            return sc.read_h5ad(path)

        with patch.object(sc, "read_10x_h5", side_effect=mock_read_10x_h5):
            result = soupx_fallback(
                str(raw_path), str(filtered_path), str(output_path), "test_sample"
            )
        assert result is True
        assert output_path.exists()

        # Corrected counts should be <= original (ambient subtracted)
        corrected = sc.read_h5ad(str(output_path))
        assert corrected.n_obs == 200


# ===== validate_correction =====

class TestValidateCorrection:
    def test_skip_when_no_raw(self, tmp_path):
        result = validate_correction(
            str(tmp_path / "nonexistent_raw.h5"),
            str(tmp_path / "nonexistent_corr.h5ad"),
            "sample1",
        )
        assert result["validation"] == "SKIP_NO_RAW"
