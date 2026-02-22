#!/usr/bin/env python3
"""
Phase 1: CellBender ambient RNA removal (GPU required).

Runs CellBender remove-background per-library on Cell Ranger raw matrices.
Falls back to a simple SoupX-like approach if CellBender fails to converge.

Usage:
    python pipeline/01_cellbender.py [--config pipeline/config.yaml]
    python pipeline/01_cellbender.py --sample 008_216_V01_L1  # Single sample
"""

import argparse
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

from utils import (
    load_config, load_metadata, get_gex_samples, setup_logging,
    setup_paths, find_cellranger_output, save_adata,
)

logger = setup_logging("01_cellbender")


def run_cellbender_sample(
    sample_name: str,
    raw_h5_path: str,
    output_h5_path: str,
    cfg: dict,
    expected_cells: int | None = None,
) -> bool:
    """
    Run CellBender remove-background on a single sample.

    Returns True if successful, False if failed.
    """
    cb_cfg = cfg["cellbender"]
    fpr = cb_cfg["fpr"]
    total_droplets = cb_cfg["total_droplets_included"]
    epochs = cb_cfg["epochs"]
    lr = cb_cfg["learning_rate"]

    cmd = [
        "cellbender", "remove-background",
        "--input", str(raw_h5_path),
        "--output", str(output_h5_path),
        "--fpr", str(fpr),
        "--total-droplets-included", str(total_droplets),
        "--epochs", str(epochs),
        "--learning-rate", str(lr),
        "--cuda",
    ]

    if expected_cells is not None:
        cmd.extend(["--expected-cells", str(expected_cells)])

    logger.info(f"Running CellBender on {sample_name}")
    logger.info(f"  Input: {raw_h5_path}")
    logger.info(f"  Output: {output_h5_path}")
    logger.info(f"  Command: {' '.join(cmd)}")

    try:
        result = subprocess.run(
            cmd, capture_output=True, text=True, timeout=3600  # 1hr timeout
        )
        if result.returncode != 0:
            logger.error(f"CellBender failed for {sample_name}:")
            logger.error(result.stderr[-2000:] if len(result.stderr) > 2000 else result.stderr)
            return False
        logger.info(f"CellBender completed for {sample_name}")
        return True
    except subprocess.TimeoutExpired:
        logger.error(f"CellBender timed out for {sample_name}")
        return False
    except FileNotFoundError:
        logger.error("CellBender not found. Install with: pip install cellbender")
        return False


def check_cellbender_convergence(output_h5_path: str, sample_name: str) -> bool:
    """
    Check if CellBender training converged by examining the output.
    Returns True if the output file exists and is valid.
    """
    output_path = Path(output_h5_path)
    # CellBender appends _filtered.h5 to the output
    filtered_path = output_path.with_name(
        output_path.stem + "_filtered.h5"
    )

    if not filtered_path.exists() and not output_path.exists():
        logger.warning(f"No CellBender output found for {sample_name}")
        return False

    # Try to read and do a basic sanity check
    try:
        check_path = filtered_path if filtered_path.exists() else output_path
        adata = sc.read_h5ad(str(check_path)) if str(check_path).endswith(".h5ad") else sc.read_10x_h5(str(check_path))
        if adata.n_obs < 100:
            logger.warning(f"CellBender output has only {adata.n_obs} cells for {sample_name}")
            return False
        return True
    except Exception as e:
        logger.warning(f"Cannot read CellBender output for {sample_name}: {e}")
        return False


def soupx_fallback(
    raw_h5_path: str,
    filtered_h5_path: str,
    output_path: str,
    sample_name: str,
) -> bool:
    """
    Simple ambient RNA removal fallback when CellBender fails.

    Uses the ratio of empty droplet profile to subtract ambient contamination.
    This is a simplified version of the SoupX approach implemented in pure Python.
    """
    logger.info(f"Running SoupX-style fallback for {sample_name}")

    try:
        raw = sc.read_10x_h5(raw_h5_path)
        filtered = sc.read_10x_h5(filtered_h5_path)

        # Identify empty droplets (in raw but not in filtered)
        filtered_barcodes = set(filtered.obs_names)
        empty_mask = ~np.isin(raw.obs_names, list(filtered_barcodes))
        empty = raw[empty_mask]

        if empty.n_obs < 50:
            logger.warning(f"Too few empty droplets for {sample_name}, skipping correction")
            save_adata(filtered, output_path)
            return True

        # Estimate ambient profile from empty droplets
        from scipy.sparse import issparse
        if issparse(empty.X):
            ambient_profile = np.asarray(empty.X.sum(axis=0)).flatten()
        else:
            ambient_profile = empty.X.sum(axis=0).flatten()
        ambient_profile = ambient_profile / ambient_profile.sum()

        # Conservative contamination fraction estimate
        rho = 0.05  # 5% contamination - conservative default

        # Subtract estimated ambient counts
        if issparse(filtered.X):
            corrected = filtered.X.toarray().astype(np.float64)
        else:
            corrected = filtered.X.astype(np.float64).copy()

        total_counts = corrected.sum(axis=1, keepdims=True)
        ambient_counts = total_counts * rho * ambient_profile[np.newaxis, :]
        corrected = corrected - ambient_counts
        corrected = np.maximum(corrected, 0)  # Floor at zero
        corrected = np.round(corrected).astype(np.float32)

        from scipy.sparse import csr_matrix
        filtered.X = csr_matrix(corrected)

        save_adata(filtered, output_path)
        logger.info(f"SoupX fallback completed for {sample_name}")
        return True

    except Exception as e:
        logger.error(f"SoupX fallback failed for {sample_name}: {e}")
        return False


def validate_correction(
    raw_h5_path: str,
    corrected_path: str,
    sample_name: str,
) -> dict:
    """
    Validate CellBender correction by checking marker gene expression.
    Returns a dict of validation metrics.
    """
    markers = {
        "T_cell": ["CD3E", "CD3D"],
        "B_cell": ["CD19", "MS4A1"],
        "Monocyte": ["CD14", "LYZ"],
    }

    try:
        filtered_path = Path(raw_h5_path).parent.parent / "filtered_feature_bc_matrix.h5"
        if filtered_path.exists():
            raw = sc.read_10x_h5(str(filtered_path))
        else:
            return {"sample": sample_name, "validation": "SKIP_NO_RAW"}

        corrected = sc.read_h5ad(corrected_path) if corrected_path.endswith(".h5ad") else sc.read_10x_h5(corrected_path)

        results = {"sample": sample_name}
        for cell_type, genes in markers.items():
            for gene in genes:
                if gene in raw.var_names and gene in corrected.var_names:
                    from scipy.sparse import issparse
                    raw_expr = raw[:, gene].X
                    corr_expr = corrected[:, gene].X
                    if issparse(raw_expr):
                        raw_expr = raw_expr.toarray()
                    if issparse(corr_expr):
                        corr_expr = corr_expr.toarray()
                    raw_mean = float(np.mean(raw_expr))
                    corr_mean = float(np.mean(corr_expr))
                    ratio = corr_mean / raw_mean if raw_mean > 0 else 0
                    results[f"{gene}_ratio"] = round(ratio, 3)
                    if ratio < 0.3:
                        logger.warning(
                            f"{sample_name}: {gene} expression ratio = {ratio:.3f} "
                            "(possible over-correction)"
                        )

        return results
    except Exception as e:
        return {"sample": sample_name, "validation": f"ERROR: {e}"}


def main():
    parser = argparse.ArgumentParser(description="CellBender ambient RNA removal")
    parser.add_argument("--config", default="pipeline/config.yaml")
    parser.add_argument("--sample", default=None, help="Run single sample")
    parser.add_argument("--no-gpu", action="store_true", help="Run without CUDA")
    parser.add_argument("--validate-only", action="store_true",
                        help="Only validate existing CellBender outputs")
    args = parser.parse_args()

    cfg = load_config(args.config)
    setup_paths(cfg)
    meta = load_metadata(cfg)
    gex = get_gex_samples(meta)

    # Deduplicate to unique sample names (excluding missing batch duplicates)
    missing_batch = cfg["study"]["missing_batch"]
    gex = gex[gex["Batch"] != missing_batch]
    sample_names = gex["sampleName"].unique()

    if args.sample:
        sample_names = [args.sample]

    cb_output_dir = Path(cfg["paths"]["cellbender_output"])
    cb_output_dir.mkdir(parents=True, exist_ok=True)

    # Track results
    results = []
    failed_samples = []

    for sample_name in sample_names:
        logger.info(f"\n{'='*60}")
        logger.info(f"Processing: {sample_name}")
        logger.info(f"{'='*60}")

        cr_paths = find_cellranger_output(sample_name, cfg)
        raw_h5 = cr_paths["raw_h5"]
        filtered_h5 = cr_paths["filtered_h5"]

        if not Path(raw_h5).exists():
            logger.error(f"Raw matrix not found: {raw_h5}")
            results.append({"sample": sample_name, "status": "NO_INPUT"})
            continue

        output_h5 = cb_output_dir / f"{sample_name}_cellbender.h5"
        output_h5ad = cb_output_dir / f"{sample_name}_cellbender.h5ad"

        if args.validate_only:
            val = validate_correction(str(raw_h5), str(output_h5ad), sample_name)
            results.append(val)
            continue

        # Skip if already done
        if output_h5ad.exists():
            logger.info(f"Output exists, skipping: {output_h5ad}")
            results.append({"sample": sample_name, "status": "EXISTS"})
            continue

        # Get expected cells from Cell Ranger output if available
        expected_cells = None
        metrics_path = cr_paths["metrics"]
        if Path(metrics_path).exists():
            try:
                metrics_df = pd.read_csv(metrics_path)
                for col in metrics_df.columns:
                    if "estimated" in col.lower() and "cell" in col.lower():
                        val = str(metrics_df[col].iloc[0]).replace(",", "")
                        expected_cells = int(float(val))
                        break
            except Exception:
                pass

        # Run CellBender
        success = run_cellbender_sample(
            sample_name, str(raw_h5), str(output_h5), cfg, expected_cells
        )

        if success:
            # Check convergence
            converged = check_cellbender_convergence(str(output_h5), sample_name)
            if converged:
                # Convert to h5ad for downstream
                try:
                    # CellBender outputs _filtered.h5
                    cb_filtered = output_h5.with_name(output_h5.stem + "_filtered.h5")
                    if cb_filtered.exists():
                        adata = sc.read_10x_h5(str(cb_filtered))
                    else:
                        adata = sc.read_10x_h5(str(output_h5))
                    adata.obs["sample_id"] = sample_name
                    save_adata(adata, str(output_h5ad))
                    results.append({"sample": sample_name, "status": "CELLBENDER_OK"})
                except Exception as e:
                    logger.error(f"Failed to convert CellBender output: {e}")
                    results.append({"sample": sample_name, "status": f"CONVERT_ERROR: {e}"})
            else:
                logger.warning(f"CellBender did not converge for {sample_name}, trying fallback")
                failed_samples.append(sample_name)
        else:
            failed_samples.append(sample_name)

    # Run SoupX fallback for failed samples
    for sample_name in failed_samples:
        cr_paths = find_cellranger_output(sample_name, cfg)
        output_h5ad = cb_output_dir / f"{sample_name}_cellbender.h5ad"
        success = soupx_fallback(
            str(cr_paths["raw_h5"]),
            str(cr_paths["filtered_h5"]),
            str(output_h5ad),
            sample_name,
        )
        results.append({
            "sample": sample_name,
            "status": "SOUPX_FALLBACK_OK" if success else "FAILED",
        })

    # Save summary
    results_df = pd.DataFrame(results)
    summary_path = cb_output_dir / "cellbender_summary.csv"
    results_df.to_csv(summary_path, index=False)
    logger.info(f"\nSummary saved to {summary_path}")
    logger.info(f"Total: {len(results_df)}, "
                f"CellBender OK: {(results_df['status'] == 'CELLBENDER_OK').sum()}, "
                f"Fallback: {(results_df['status'] == 'SOUPX_FALLBACK_OK').sum()}, "
                f"Failed: {(results_df['status'] == 'FAILED').sum()}")


if __name__ == "__main__":
    main()
