#!/usr/bin/env python3
"""
Phase 2: Per-library QC filtering + Scrublet doublet detection.

Applies MAD-based adaptive thresholds per-library for:
  - log1p(total_counts): both directions, 3 MAD
  - log1p(n_genes_by_counts): both directions, 3 MAD
  - pct_counts_mt: upper only, 3 MAD
  - pct_counts_hb: hard ceiling at 5%

Also runs Scrublet for doublet detection and applies gene-level filters.

Usage:
    python pipeline/02_qc_filtering.py [--config pipeline/config.yaml]
    python pipeline/02_qc_filtering.py --sample 008_216_V01_L1
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
from scipy.sparse import issparse

from utils import (
    load_config, load_metadata, get_gex_samples, setup_logging,
    setup_paths, calculate_qc_metrics, mad_outlier, save_adata,
    plot_qc_violins, plot_filtering_summary,
)

logger = setup_logging("02_qc_filtering")


def filter_library_adaptive(
    adata: sc.AnnData,
    sample_name: str,
    cfg: dict,
) -> sc.AnnData:
    """
    Apply MAD-based adaptive QC filtering to a single library.

    Returns filtered AnnData with QC annotations in .obs.
    """
    nmads = cfg["qc"]["mad_threshold"]
    max_hb = cfg["qc"]["max_pct_hb"]

    n_before = adata.n_obs
    logger.info(f"  {sample_name}: {n_before} cells before filtering")

    # Calculate QC metrics
    adata = calculate_qc_metrics(adata)

    # --- Outlier detection (all on log1p-transformed where appropriate) ---
    # Total counts
    outlier_counts = mad_outlier(
        adata.obs["log1p_total_counts"].values, nmads=nmads, direction="both"
    )
    # Gene counts
    outlier_genes = mad_outlier(
        adata.obs["log1p_n_genes_by_counts"].values, nmads=nmads, direction="both"
    )
    # Mito % (upper only - low mito is fine)
    outlier_mt = mad_outlier(
        adata.obs["pct_counts_mt"].values, nmads=nmads, direction="upper"
    )
    # Hemoglobin - hard ceiling
    outlier_hb = adata.obs["pct_counts_hb"].values > max_hb

    # Combined outlier flag
    outlier = outlier_counts | outlier_genes | outlier_mt | outlier_hb

    # Store flags before filtering
    adata.obs["outlier_counts"] = outlier_counts
    adata.obs["outlier_genes"] = outlier_genes
    adata.obs["outlier_mt"] = outlier_mt
    adata.obs["outlier_hb"] = outlier_hb
    adata.obs["outlier_any"] = outlier

    # Log threshold values for reference
    log_counts = adata.obs["log1p_total_counts"].values
    log_genes = adata.obs["log1p_n_genes_by_counts"].values
    mt_vals = adata.obs["pct_counts_mt"].values

    def _mad_bounds(values, nmads):
        med = np.nanmedian(values)
        mad = np.nanmedian(np.abs(values - med))
        if mad == 0:
            mad = np.nanstd(values)
        return med - nmads * 1.4826 * mad, med + nmads * 1.4826 * mad

    lo_c, hi_c = _mad_bounds(log_counts, nmads)
    lo_g, hi_g = _mad_bounds(log_genes, nmads)
    _, hi_mt = _mad_bounds(mt_vals, nmads)

    logger.info(f"    log1p_total_counts: [{lo_c:.2f}, {hi_c:.2f}]")
    logger.info(f"    log1p_n_genes:      [{lo_g:.2f}, {hi_g:.2f}]")
    logger.info(f"    pct_counts_mt:      <= {hi_mt:.2f}")
    logger.info(f"    pct_counts_hb:      <= {max_hb}")

    # Filter
    adata = adata[~outlier].copy()
    n_after = adata.n_obs
    logger.info(f"    Removed {n_before - n_after} cells "
                f"({(n_before - n_after)/n_before*100:.1f}%), "
                f"retained {n_after}")

    return adata


def run_scrublet(
    adata: sc.AnnData,
    sample_name: str,
    expected_doublet_rate: float = 0.046,
    random_state: int = 42,
) -> sc.AnnData:
    """
    Run Scrublet doublet detection on a single library.
    Stores doublet scores and predictions in .obs.
    """
    logger.info(f"  Running Scrublet on {sample_name} ({adata.n_obs} cells)")

    counts = adata.X
    if issparse(counts):
        counts = counts.toarray()

    try:
        scrub = scr.Scrublet(
            counts,
            expected_doublet_rate=expected_doublet_rate,
            random_state=random_state,
        )
        doublet_scores, predicted_doublets = scrub.scrub_doublets(
            min_counts=2, min_cells=3, min_gene_variability_pctl=85,
            n_prin_comps=30, verbose=False,
        )

        adata.obs["doublet_score"] = doublet_scores
        adata.obs["predicted_doublet"] = predicted_doublets

        n_doublets = predicted_doublets.sum()
        logger.info(f"    Scrublet: {n_doublets} doublets detected "
                    f"({n_doublets/adata.n_obs*100:.1f}%), "
                    f"threshold={scrub.threshold_:.3f}")

    except Exception as e:
        logger.warning(f"    Scrublet failed for {sample_name}: {e}")
        logger.warning("    Setting doublet_score=0 and predicted_doublet=False")
        adata.obs["doublet_score"] = 0.0
        adata.obs["predicted_doublet"] = False

    return adata


def filter_genes(adata: sc.AnnData, min_cells: int = 3) -> sc.AnnData:
    """
    Gene-level filtering:
    - Remove genes detected in fewer than min_cells
    - Remove mitochondrial genes from feature set (keep metrics in .obs)
    - Flag ribosomal genes (keep in data)
    """
    n_before = adata.n_vars

    # Filter by min cells
    sc.pp.filter_genes(adata, min_cells=min_cells)

    # Remove mito genes from the matrix (metrics already stored in .obs)
    mt_mask = adata.var_names.str.startswith("MT-")
    n_mt = mt_mask.sum()
    if n_mt > 0:
        adata = adata[:, ~mt_mask].copy()

    n_after = adata.n_vars
    logger.info(f"    Gene filtering: {n_before} -> {n_after} genes "
                f"(removed {n_before - n_after}, incl. {n_mt} MT genes)")

    return adata


def process_single_sample(
    sample_name: str,
    cfg: dict,
) -> dict:
    """Process a single sample through QC + doublet detection."""
    cb_output = Path(cfg["paths"]["cellbender_output"])
    qc_output = Path(cfg["paths"]["qc_output"])

    input_path = cb_output / f"{sample_name}_cellbender.h5ad"
    output_path = qc_output / f"{sample_name}_qc.h5ad"

    if not input_path.exists():
        logger.error(f"Input not found: {input_path}")
        return {"sample": sample_name, "status": "NO_INPUT"}

    if output_path.exists():
        logger.info(f"Output exists, skipping: {output_path}")
        adata = sc.read_h5ad(output_path)
        return {
            "sample": sample_name,
            "status": "EXISTS",
            "cells_after": adata.n_obs,
        }

    # Load CellBender output
    adata = sc.read_h5ad(str(input_path))
    adata.obs["sample_id"] = sample_name
    adata.var_names_make_unique()
    n_cells_raw = adata.n_obs

    # QC filtering (MAD-based, adaptive per-library)
    adata = filter_library_adaptive(adata, sample_name, cfg)

    # Gene-level filtering
    adata = filter_genes(adata, min_cells=cfg["qc"]["min_cells_per_gene"])

    # Doublet detection
    adata = run_scrublet(
        adata, sample_name,
        expected_doublet_rate=cfg["qc"]["expected_doublet_rate"],
        random_state=cfg["compute"]["random_seed"],
    )

    # Remove predicted doublets
    n_before_doublet = adata.n_obs
    adata = adata[~adata.obs["predicted_doublet"]].copy()
    n_doublets_removed = n_before_doublet - adata.n_obs
    logger.info(f"    Removed {n_doublets_removed} predicted doublets")

    # Check minimum cells threshold
    min_cells = cfg["qc"]["min_cells_per_library"]
    if adata.n_obs < min_cells:
        logger.warning(
            f"    {sample_name}: only {adata.n_obs} cells after QC "
            f"(below threshold of {min_cells})"
        )

    # Save
    save_adata(adata, str(output_path))

    return {
        "sample": sample_name,
        "status": "OK",
        "cells_raw": n_cells_raw,
        "cells_after_qc": n_before_doublet,
        "doublets_removed": n_doublets_removed,
        "cells_final": adata.n_obs,
        "pct_retained": round(adata.n_obs / n_cells_raw * 100, 1) if n_cells_raw > 0 else 0,
        "n_genes": adata.n_vars,
    }


def main():
    parser = argparse.ArgumentParser(description="Per-library QC + doublet detection")
    parser.add_argument("--config", default="pipeline/config.yaml")
    parser.add_argument("--sample", default=None, help="Process single sample")
    args = parser.parse_args()

    cfg = load_config(args.config)
    setup_paths(cfg)
    meta = load_metadata(cfg)
    gex = get_gex_samples(meta)

    # Deduplicate sample names (exclude missing batch)
    missing_batch = cfg["study"]["missing_batch"]
    gex = gex[gex["Batch"] != missing_batch]
    sample_names = gex["sampleName"].unique()

    if args.sample:
        sample_names = [args.sample]

    logger.info(f"Processing {len(sample_names)} samples for QC")

    # Process each sample
    all_results = []
    pre_counts = {}
    post_counts = {}

    for sample_name in sample_names:
        logger.info(f"\n{'='*60}")
        logger.info(f"Sample: {sample_name}")
        logger.info(f"{'='*60}")

        result = process_single_sample(sample_name, cfg)
        all_results.append(result)

        if "cells_raw" in result:
            pre_counts[sample_name] = result["cells_raw"]
            post_counts[sample_name] = result["cells_final"]

    # Save summary
    qc_output = Path(cfg["paths"]["qc_output"])
    results_df = pd.DataFrame(all_results)
    summary_path = qc_output / "qc_summary.csv"
    results_df.to_csv(summary_path, index=False)
    logger.info(f"\nQC summary saved to {summary_path}")

    # Report
    ok_results = results_df[results_df["status"] == "OK"]
    if len(ok_results) > 0:
        logger.info(f"Processed: {len(ok_results)} samples")
        logger.info(f"Total cells retained: {ok_results['cells_final'].sum():,}")
        logger.info(f"Median retention: {ok_results['pct_retained'].median():.1f}%")

        # Flag samples that lost >50% cells
        low_retention = ok_results[ok_results["pct_retained"] < 50]
        if len(low_retention) > 0:
            logger.warning(f"\n{len(low_retention)} samples lost >50% cells:")
            for _, row in low_retention.iterrows():
                logger.warning(f"  {row['sample']}: {row['pct_retained']}% retained")

    # Generate QC summary plots
    if pre_counts and post_counts:
        fig_dir = Path(cfg["paths"]["figures"]) / "qc"
        fig_dir.mkdir(parents=True, exist_ok=True)
        plot_filtering_summary(
            pre_counts, post_counts,
            save_path=str(fig_dir / "filtering_summary.png"),
        )
        logger.info(f"Filtering summary plot saved to {fig_dir}/filtering_summary.png")


if __name__ == "__main__":
    main()
