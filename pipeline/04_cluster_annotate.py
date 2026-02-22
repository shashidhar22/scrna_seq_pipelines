#!/usr/bin/env python3
"""
Phase 4: Cell type annotation with CellTypist + manual marker validation.

Steps:
  1. Load integrated AnnData
  2. Run CellTypist (Immune_All_Low + Immune_All_High)
  3. Validate with canonical immune markers (dot plots)
  4. Refine clusters: merge/split based on marker profiles
  5. Store final annotations in cell_type (coarse) and cell_subtype (fine)

Usage:
    python pipeline/04_cluster_annotate.py [--config pipeline/config.yaml]
    python pipeline/04_cluster_annotate.py --arm primary
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

from utils import (
    load_config, setup_logging, setup_paths, save_adata,
    load_adata, set_plotting_defaults,
)

logger = setup_logging("04_cluster_annotate")

# Canonical immune marker sets for validation
CANONICAL_MARKERS = {
    "T cells": ["CD3D", "CD3E"],
    "CD4 T": ["CD4", "IL7R"],
    "CD8 T": ["CD8A", "CD8B"],
    "NK": ["NKG7", "GNLY", "KLRD1"],
    "B cells": ["CD19", "MS4A1", "CD79A"],
    "Monocytes CD14+": ["CD14", "LYZ"],
    "Monocytes CD16+": ["FCGR3A", "MS4A7"],
    "DCs": ["FCER1A", "CLEC10A"],
    "pDCs": ["IRF7", "LILRA4", "IL3RA"],
    "Platelets": ["PPBP", "PF4"],
    "Plasma": ["MZB1", "SDC1", "JCHAIN"],
    "Proliferating": ["MKI67", "TOP2A"],
}


def run_celltypist(adata: sc.AnnData, cfg: dict) -> sc.AnnData:
    """Run CellTypist with configured models."""
    import celltypist

    models = cfg["annotation"]["celltypist_models"]
    majority_voting = cfg["annotation"]["majority_voting"]

    for model_name in models:
        logger.info(f"Running CellTypist with model: {model_name}")

        # CellTypist expects normalized, log-transformed data
        # Use .raw if available, otherwise use current .X
        model = celltypist.models.Model.load(model=model_name)

        predictions = celltypist.annotate(
            adata, model=model, majority_voting=majority_voting
        )

        # Store predictions
        short_name = model_name.replace(".pkl", "")
        adata.obs[f"celltypist_{short_name}"] = (
            predictions.predicted_labels.predicted_labels.values
        )
        if majority_voting:
            adata.obs[f"celltypist_{short_name}_majority"] = (
                predictions.predicted_labels.majority_voting.values
            )

        # Log distribution
        col = f"celltypist_{short_name}_majority" if majority_voting else f"celltypist_{short_name}"
        value_counts = adata.obs[col].value_counts()
        logger.info(f"  {model_name}: {len(value_counts)} cell types")
        for ct, count in value_counts.head(10).items():
            logger.info(f"    {ct}: {count} ({count/adata.n_obs*100:.1f}%)")

    return adata


def validate_markers(adata: sc.AnnData, cfg: dict, cluster_key: str = "leiden_0.8",
                     save_dir: str | None = None) -> None:
    """Generate marker validation plots (dot plots, feature plots)."""
    # Collect all markers that exist in the data
    all_markers = []
    for cell_type, genes in CANONICAL_MARKERS.items():
        for gene in genes:
            if gene in adata.var_names:
                all_markers.append(gene)
            elif adata.raw is not None and gene in adata.raw.var_names:
                all_markers.append(gene)

    all_markers = list(dict.fromkeys(all_markers))  # Deduplicate preserving order

    if not all_markers:
        logger.warning("No canonical markers found in dataset")
        return

    if save_dir:
        save_path = Path(save_dir)
        save_path.mkdir(parents=True, exist_ok=True)

    # Dot plot by cluster
    logger.info(f"Generating marker dot plot for {cluster_key}...")
    fig, ax = plt.subplots(figsize=(max(14, len(all_markers) * 0.6),
                                     max(6, adata.obs[cluster_key].nunique() * 0.4)))
    sc.pl.dotplot(adata, var_names=all_markers, groupby=cluster_key,
                  standard_scale="var", show=False, ax=ax)
    plt.tight_layout()
    if save_dir:
        fig.savefig(save_path / f"dotplot_markers_{cluster_key}.png",
                    bbox_inches="tight", dpi=150)
        plt.close(fig)

    # Dot plot by CellTypist annotation
    ct_col = None
    for col in adata.obs.columns:
        if "celltypist" in col and "majority" in col:
            ct_col = col
            break

    if ct_col:
        fig, ax = plt.subplots(figsize=(max(14, len(all_markers) * 0.6),
                                         max(8, adata.obs[ct_col].nunique() * 0.3)))
        sc.pl.dotplot(adata, var_names=all_markers, groupby=ct_col,
                      standard_scale="var", show=False, ax=ax)
        plt.tight_layout()
        if save_dir:
            fig.savefig(save_path / f"dotplot_markers_{ct_col}.png",
                        bbox_inches="tight", dpi=150)
            plt.close(fig)

    # UMAP feature plots for key markers
    key_markers = ["CD3E", "CD4", "CD8A", "MS4A1", "CD14", "NKG7", "PPBP", "IRF7"]
    available_markers = [m for m in key_markers if m in adata.var_names or
                         (adata.raw is not None and m in adata.raw.var_names)]

    if available_markers:
        n_cols = 4
        n_rows = (len(available_markers) + n_cols - 1) // n_cols
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(4 * n_cols, 4 * n_rows))
        axes = np.atleast_2d(axes)
        for idx, gene in enumerate(available_markers):
            row, col = divmod(idx, n_cols)
            sc.pl.umap(adata, color=gene, ax=axes[row, col], show=False,
                       title=gene, frameon=False)
        # Hide empty axes
        for idx in range(len(available_markers), n_rows * n_cols):
            row, col = divmod(idx, n_cols)
            axes[row, col].set_visible(False)
        plt.tight_layout()
        if save_dir:
            fig.savefig(save_path / "umap_marker_genes.png",
                        bbox_inches="tight", dpi=150)
            plt.close(fig)


def assign_coarse_labels(adata: sc.AnnData, cluster_key: str = "leiden_0.8") -> sc.AnnData:
    """
    Assign coarse cell type labels based on CellTypist majority voting.
    Maps fine-grained CellTypist labels to broad categories.
    """
    # Coarse label mapping from CellTypist fine labels
    coarse_map = {
        # T cells
        "T": "T cells", "CD4": "CD4 T", "CD8": "CD8 T", "Treg": "CD4 T",
        "Th1": "CD4 T", "Th2": "CD4 T", "Th17": "CD4 T", "Tfh": "CD4 T",
        "Tcm": "T cells", "Tem": "T cells", "Tnaive": "T cells",
        "NKT": "NK/T", "MAIT": "T cells", "gdT": "T cells",
        # NK
        "NK": "NK",
        # B cells
        "B": "B cells", "Plasma": "Plasma",
        # Myeloid
        "Mono": "Monocytes", "DC": "DCs", "pDC": "pDCs",
        "Macro": "Macrophages",
        # Other
        "Platelet": "Platelets", "Mega": "Platelets",
        "Erythro": "Erythrocytes",
        "ILC": "ILCs",
    }

    # Find CellTypist majority column
    ct_col = None
    for col in adata.obs.columns:
        if "celltypist" in col and "majority" in col and "Low" in col:
            ct_col = col
            break
    if ct_col is None:
        for col in adata.obs.columns:
            if "celltypist" in col and "majority" in col:
                ct_col = col
                break

    if ct_col is None:
        logger.warning("No CellTypist predictions found. Using cluster IDs as labels.")
        adata.obs["cell_type"] = adata.obs[cluster_key].astype(str)
        adata.obs["cell_subtype"] = adata.obs[cluster_key].astype(str)
        return adata

    # Map fine labels to coarse
    fine_labels = adata.obs[ct_col].astype(str)

    def _map_coarse(label: str) -> str:
        # Check longer keys first so specific matches (e.g., "CD4") beat
        # generic ones (e.g., "T").
        for key in sorted(coarse_map, key=len, reverse=True):
            if key.lower() in label.lower():
                return coarse_map[key]
        return "Unknown"

    adata.obs["cell_type"] = fine_labels.map(_map_coarse)
    adata.obs["cell_subtype"] = fine_labels

    # Log distribution
    logger.info("Coarse cell type distribution:")
    for ct, count in adata.obs["cell_type"].value_counts().items():
        logger.info(f"  {ct}: {count} ({count/adata.n_obs*100:.1f}%)")

    return adata


def find_marker_genes(adata: sc.AnnData, cluster_key: str = "leiden_0.8",
                      n_genes: int = 20) -> pd.DataFrame:
    """Run rank_genes_groups to find cluster markers."""
    logger.info(f"Finding marker genes per cluster ({cluster_key})...")
    sc.tl.rank_genes_groups(adata, groupby=cluster_key, method="wilcoxon",
                            use_raw=True, pts=True)

    # Extract results into a DataFrame
    result = adata.uns["rank_genes_groups"]
    groups = result["names"].dtype.names
    rows = []
    for group in groups:
        for i in range(n_genes):
            rows.append({
                "cluster": group,
                "gene": result["names"][group][i],
                "logfoldchange": result["logfoldchanges"][group][i],
                "pval_adj": result["pvals_adj"][group][i],
                "pct_in": result["pts"][group][i] if "pts" in result else None,
                "pct_out": result["pts_rest"][group][i] if "pts_rest" in result else None,
            })

    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description="Cell type annotation")
    parser.add_argument("--config", default="pipeline/config.yaml")
    parser.add_argument("--arm", default="primary")
    parser.add_argument("--cluster-key", default="leiden_0.8")
    args = parser.parse_args()

    cfg = load_config(args.config)
    setup_paths(cfg)
    set_plotting_defaults()

    # Load integrated data
    in_path = Path(cfg["paths"]["integration_output"]) / f"integrated_{args.arm}.h5ad"
    logger.info(f"Loading integrated data from {in_path}")
    adata = load_adata(str(in_path))

    # Step 1: CellTypist annotation
    adata = run_celltypist(adata, cfg)

    # Step 2: Assign coarse labels
    adata = assign_coarse_labels(adata, cluster_key=args.cluster_key)

    # Step 3: Find cluster marker genes
    markers_df = find_marker_genes(adata, cluster_key=args.cluster_key)
    markers_path = Path(cfg["paths"]["annotation_output"]) / f"cluster_markers_{args.arm}.csv"
    markers_path.parent.mkdir(parents=True, exist_ok=True)
    markers_df.to_csv(markers_path, index=False)
    logger.info(f"Cluster markers saved to {markers_path}")

    # Step 4: Validation plots
    fig_dir = Path(cfg["paths"]["figures"]) / "annotation"
    validate_markers(adata, cfg, cluster_key=args.cluster_key, save_dir=str(fig_dir))

    # UMAP by annotation
    fig, axes = plt.subplots(1, 3, figsize=(24, 7))
    sc.pl.umap(adata, color="cell_type", ax=axes[0], show=False,
               title="Cell type (coarse)")
    sc.pl.umap(adata, color=args.cluster_key, ax=axes[1], show=False,
               title=f"Leiden ({args.cluster_key})")
    ct_col = [c for c in adata.obs.columns if "celltypist" in c and "majority" in c]
    if ct_col:
        sc.pl.umap(adata, color=ct_col[0], ax=axes[2], show=False,
                   title="CellTypist (fine)")
    plt.tight_layout()
    fig.savefig(fig_dir / f"umap_annotation_{args.arm}.png",
                bbox_inches="tight", dpi=150)
    plt.close(fig)

    # Save annotated data
    out_path = Path(cfg["paths"]["annotation_output"]) / f"annotated_{args.arm}.h5ad"
    save_adata(adata, str(out_path))

    logger.info(f"\nDone. Annotated data: {out_path}")
    logger.info(f"Cell types found: {adata.obs['cell_type'].nunique()}")


if __name__ == "__main__":
    main()
