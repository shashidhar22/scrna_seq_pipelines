#!/usr/bin/env python3
"""
Phase 8: Trajectory analysis - PAGA + diffusion pseudotime.

Steps:
  1. PAGA connectivity on T cell subtypes
  2. Diffusion pseudotime (DPT) with naive T cell root
  3. Overlay clonal expansion on pseudotime
  4. Track clonotype lineages through differentiation states

Usage:
    python pipeline/08_trajectory.py [--config pipeline/config.yaml]
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

logger = setup_logging("08_trajectory")


def run_paga(adata_t: sc.AnnData, groups_key: str = "tcell_subtype") -> sc.AnnData:
    """Run PAGA to infer connectivity between T cell subtypes."""
    logger.info(f"Running PAGA on {groups_key}...")

    sc.tl.paga(adata_t, groups=groups_key)

    # Log connections
    if "paga" in adata_t.uns:
        logger.info("PAGA connectivity computed")
        connectivities = adata_t.uns["paga"]["connectivities"]
        cats = adata_t.obs[groups_key].cat.categories if hasattr(adata_t.obs[groups_key], "cat") else adata_t.obs[groups_key].unique()
        logger.info(f"  Groups: {list(cats)}")

    return adata_t


def run_diffusion_pseudotime(
    adata_t: sc.AnnData,
    cfg: dict,
    root_subtype: str | None = None,
) -> sc.AnnData:
    """
    Compute diffusion pseudotime.
    Sets root to the most naive-like cell.
    """
    traj_cfg = cfg["trajectory"]

    logger.info("Computing diffusion map...")
    sc.tl.diffmap(adata_t, n_comps=traj_cfg["n_dcs"])

    # Determine root cell
    root_type = root_subtype or traj_cfg.get("root_cell_type", "Naive")

    if "tcell_subtype" in adata_t.obs.columns:
        # Find naive T cell cluster
        naive_mask = adata_t.obs["tcell_subtype"].str.contains(root_type, case=False, na=False)
        if naive_mask.any():
            # Pick the cell with lowest diffusion component 1 among naive cells
            naive_idx = np.where(naive_mask)[0]
            dc1_naive = adata_t.obsm["X_diffmap"][naive_idx, 0]
            root_cell_local = naive_idx[np.argmin(dc1_naive)]
            adata_t.uns["iroot"] = root_cell_local
            logger.info(f"Root cell set to index {root_cell_local} "
                        f"(from '{root_type}' cells)")
        else:
            logger.warning(f"No '{root_type}' cells found. Using first cell as root.")
            adata_t.uns["iroot"] = 0
    else:
        adata_t.uns["iroot"] = 0

    logger.info("Computing diffusion pseudotime...")
    sc.tl.dpt(adata_t)

    logger.info(f"DPT range: {adata_t.obs['dpt_pseudotime'].min():.3f} - "
                f"{adata_t.obs['dpt_pseudotime'].max():.3f}")

    return adata_t


def analyze_expansion_along_pseudotime(adata_t: sc.AnnData) -> pd.DataFrame | None:
    """Analyze clonal expansion as a function of pseudotime."""
    if "dpt_pseudotime" not in adata_t.obs.columns:
        return None
    if "clonal_expansion" not in adata_t.obs.columns:
        return None

    # Bin pseudotime and compute expansion fraction per bin
    adata_t.obs["dpt_bin"] = pd.cut(adata_t.obs["dpt_pseudotime"], bins=20, labels=False)

    results = []
    for bin_val, group in adata_t.obs.groupby("dpt_bin"):
        total = len(group)
        expanded = (group["clonal_expansion"] != "1").sum() if "1" in group["clonal_expansion"].values else 0
        results.append({
            "pseudotime_bin": bin_val,
            "mean_pseudotime": group["dpt_pseudotime"].mean(),
            "n_cells": total,
            "n_expanded": expanded,
            "frac_expanded": expanded / total if total > 0 else 0,
        })

    return pd.DataFrame(results)


def track_clonotype_trajectories(adata_t: sc.AnnData, top_n: int = 20) -> pd.DataFrame | None:
    """
    Track top expanded clonotypes along pseudotime.
    Returns per-clonotype pseudotime statistics.
    """
    if "clone_id" not in adata_t.obs.columns or "dpt_pseudotime" not in adata_t.obs.columns:
        return None

    # Get top expanded clonotypes
    clone_counts = adata_t.obs["clone_id"].value_counts()
    top_clones = clone_counts[clone_counts > 1].head(top_n).index

    results = []
    for clone in top_clones:
        mask = adata_t.obs["clone_id"] == clone
        cells = adata_t.obs.loc[mask]
        dpt = cells["dpt_pseudotime"]

        subtypes = cells["tcell_subtype"].value_counts() if "tcell_subtype" in cells.columns else pd.Series(dtype=int)
        dominant_subtype = subtypes.index[0] if len(subtypes) > 0 else "unknown"

        results.append({
            "clone_id": clone,
            "n_cells": mask.sum(),
            "mean_pseudotime": dpt.mean(),
            "min_pseudotime": dpt.min(),
            "max_pseudotime": dpt.max(),
            "pseudotime_span": dpt.max() - dpt.min(),
            "dominant_subtype": dominant_subtype,
            "n_subtypes": cells["tcell_subtype"].nunique() if "tcell_subtype" in cells.columns else 0,
        })

    return pd.DataFrame(results)


def plot_trajectory_results(adata_t: sc.AnnData, cfg: dict, save_dir: str) -> None:
    """Generate trajectory analysis plots."""
    save_path = Path(save_dir)
    save_path.mkdir(parents=True, exist_ok=True)

    # PAGA graph
    if "paga" in adata_t.uns:
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        sc.pl.paga(adata_t, ax=axes[0], show=False, title="PAGA connectivity",
                   threshold=0.1)
        sc.pl.umap(adata_t, color="tcell_subtype", ax=axes[1], show=False,
                   title="T cell subtypes")
        plt.tight_layout()
        fig.savefig(save_path / "paga_connectivity.png", bbox_inches="tight", dpi=150)
        plt.close(fig)

    # Pseudotime on UMAP
    if "dpt_pseudotime" in adata_t.obs.columns:
        fig, axes = plt.subplots(1, 3, figsize=(21, 6))
        sc.pl.umap(adata_t, color="dpt_pseudotime", ax=axes[0], show=False,
                   title="Diffusion pseudotime", cmap="viridis")
        sc.pl.umap(adata_t, color="tcell_subtype", ax=axes[1], show=False,
                   title="T cell subtype")
        if "clonal_expansion" in adata_t.obs.columns:
            sc.pl.umap(adata_t, color="clonal_expansion", ax=axes[2], show=False,
                       title="Clonal expansion")
        else:
            axes[2].set_visible(False)
        plt.tight_layout()
        fig.savefig(save_path / "pseudotime_umap.png", bbox_inches="tight", dpi=150)
        plt.close(fig)

        # Pseudotime distribution per subtype
        if "tcell_subtype" in adata_t.obs.columns:
            fig, ax = plt.subplots(figsize=(10, 6))
            subtypes = adata_t.obs["tcell_subtype"].unique()
            for st in sorted(subtypes):
                mask = adata_t.obs["tcell_subtype"] == st
                vals = adata_t.obs.loc[mask, "dpt_pseudotime"].dropna()
                if len(vals) > 0:
                    ax.hist(vals, bins=30, alpha=0.4, label=st, density=True)
            ax.set_xlabel("Pseudotime")
            ax.set_ylabel("Density")
            ax.set_title("Pseudotime distribution per T cell subtype")
            ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
            plt.tight_layout()
            fig.savefig(save_path / "pseudotime_by_subtype.png",
                        bbox_inches="tight", dpi=150)
            plt.close(fig)

    # Expansion along pseudotime
    expansion_df = analyze_expansion_along_pseudotime(adata_t)
    if expansion_df is not None and len(expansion_df) > 0:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(expansion_df["mean_pseudotime"], expansion_df["frac_expanded"],
                "o-", color="firebrick")
        ax.set_xlabel("Pseudotime")
        ax.set_ylabel("Fraction expanded clones")
        ax.set_title("Clonal expansion along differentiation trajectory")
        plt.tight_layout()
        fig.savefig(save_path / "expansion_vs_pseudotime.png",
                    bbox_inches="tight", dpi=150)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Trajectory analysis")
    parser.add_argument("--config", default="pipeline/config.yaml")
    parser.add_argument("--arm", default="primary")
    args = parser.parse_args()

    cfg = load_config(args.config)
    setup_paths(cfg)
    set_plotting_defaults()

    # Load VDJ-integrated T cell data (or T cell subcluster data)
    vdj_path = Path(cfg["paths"]["vdj_output"]) / f"tcell_vdj_{args.arm}.h5ad"
    tcell_path = Path(cfg["paths"]["tcell_output"]) / f"tcell_{args.arm}.h5ad"

    if vdj_path.exists():
        logger.info(f"Loading VDJ-integrated T cell data from {vdj_path}")
        adata_t = load_adata(str(vdj_path))
    elif tcell_path.exists():
        logger.info(f"Loading T cell data from {tcell_path}")
        adata_t = load_adata(str(tcell_path))
    else:
        logger.error("No T cell data found. Run phase 5 or 6 first.")
        return

    # Step 1: PAGA
    groups_key = "tcell_subtype" if "tcell_subtype" in adata_t.obs.columns else "leiden_1.0"
    adata_t = run_paga(adata_t, groups_key=groups_key)

    # Step 2: Diffusion pseudotime
    adata_t = run_diffusion_pseudotime(adata_t, cfg)

    # Step 3: Clonotype trajectory tracking
    clone_traj = track_clonotype_trajectories(adata_t)
    if clone_traj is not None:
        traj_output = Path(cfg["paths"]["trajectory_output"])
        traj_output.mkdir(parents=True, exist_ok=True)
        clone_traj.to_csv(traj_output / f"clonotype_trajectories_{args.arm}.csv", index=False)
        logger.info(f"Clonotype trajectory data saved")

        # Log top clones spanning multiple subtypes
        multi_subtype = clone_traj[clone_traj["n_subtypes"] > 1]
        if len(multi_subtype) > 0:
            logger.info(f"\nClones spanning multiple subtypes ({len(multi_subtype)}):")
            for _, row in multi_subtype.head(10).iterrows():
                logger.info(f"  {row['clone_id']}: {row['n_cells']} cells, "
                            f"{row['n_subtypes']} subtypes, "
                            f"span={row['pseudotime_span']:.3f}")

    # Step 4: Plots
    fig_dir = Path(cfg["paths"]["figures"]) / "trajectory"
    plot_trajectory_results(adata_t, cfg, str(fig_dir))

    # Save
    out_dir = Path(cfg["paths"]["trajectory_output"])
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"trajectory_{args.arm}.h5ad"
    save_adata(adata_t, str(out_path))

    logger.info(f"\nDone. Trajectory data: {out_path}")


if __name__ == "__main__":
    main()
