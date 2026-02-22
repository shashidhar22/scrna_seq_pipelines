#!/usr/bin/env python3
"""
Phase 7: Pseudobulk differential expression with PyDESeq2.

Steps:
  1. Aggregate counts to pseudobulk (per patient per cell type)
  2. Run DESeq2 per cell type: HIV+ vs HIV-
  3. Optional: scVI-based DE as complementary approach
  4. Cell type proportion analysis (compositional)
  5. Volcano plots and results export

Usage:
    python pipeline/07_differential_expression.py [--config pipeline/config.yaml]
    python pipeline/07_differential_expression.py --contrast HIV_status
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.sparse import issparse

from utils import (
    load_config, setup_logging, setup_paths, save_adata,
    load_adata, set_plotting_defaults,
)

logger = setup_logging("07_differential_expression")


def compute_pseudobulk(
    adata: sc.AnnData,
    sample_col: str,
    groups_col: str,
    layer: str = "counts",
    min_cells: int = 10,
    min_counts: int = 1000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Aggregate single-cell counts to pseudobulk per sample per cell type.

    Returns:
        counts_df: Genes x samples pseudobulk count matrix
        meta_df: Sample-level metadata
    """
    logger.info(f"Computing pseudobulk (group by {sample_col} x {groups_col})...")

    # Use decoupler if available, otherwise manual aggregation
    try:
        import decoupler as dc
        pdata = dc.get_pseudobulk(
            adata, sample_col=sample_col, groups_col=groups_col,
            layer=layer, min_cells=min_cells, min_counts=min_counts,
        )
        logger.info(f"Pseudobulk: {pdata.n_obs} samples x {pdata.n_vars} genes")
        return pdata
    except ImportError:
        logger.info("decoupler not available, using manual aggregation")

    # Manual pseudobulk aggregation
    X = adata.layers[layer] if layer in adata.layers else adata.X

    groups = adata.obs.groupby([sample_col, groups_col])
    pb_counts = {}
    pb_meta = []

    for (sample, cell_type), obs_names in groups.groups.items():
        if len(obs_names) < min_cells:
            continue

        # Convert obs_names to integer positions for sparse matrix indexing
        idx = np.where(adata.obs_names.isin(obs_names))[0]
        subset_X = X[idx]
        if issparse(subset_X):
            total = np.asarray(subset_X.sum(axis=0)).flatten()
        else:
            total = subset_X.sum(axis=0).flatten()

        if total.sum() < min_counts:
            continue

        key = f"{sample}__{cell_type}"
        pb_counts[key] = total

        # Get metadata from first cell in group
        row_meta = adata.obs.loc[obs_names[0]]
        pb_meta.append({
            "pb_sample": key,
            sample_col: sample,
            groups_col: cell_type,
            "n_cells": len(idx),
            "HIVstatus": row_meta.get("HIVstatus", "unknown"),
            "Batch": row_meta.get("Batch", "unknown"),
            "tissueType": row_meta.get("tissueType", "unknown"),
        })

        # Add clinical metadata if present
        for col in ["best_resp", "bl_plasma_grp", "best_resp_new"]:
            if col in adata.obs.columns:
                pb_meta[-1][col] = row_meta.get(col, "unknown")

    counts_df = pd.DataFrame(pb_counts, index=adata.var_names).T
    meta_df = pd.DataFrame(pb_meta).set_index("pb_sample")

    logger.info(f"Pseudobulk: {len(meta_df)} samples x {counts_df.shape[1]} genes")
    logger.info(f"Cell types: {meta_df[groups_col].nunique()}")

    return counts_df, meta_df


def run_deseq2(
    counts_df: pd.DataFrame,
    meta_df: pd.DataFrame,
    design_col: str,
    reference: str,
    target: str,
    cell_type: str | None = None,
    groups_col: str = "cell_type",
) -> pd.DataFrame | None:
    """
    Run PyDESeq2 for a specific contrast.
    If cell_type is specified, subset to that cell type first.
    """
    from pydeseq2.dds import DeseqDataSet
    from pydeseq2.ds import DeseqStats

    # Subset to cell type if specified
    if cell_type:
        mask = meta_df[groups_col] == cell_type
        counts_sub = counts_df.loc[mask]
        meta_sub = meta_df.loc[mask]
    else:
        counts_sub = counts_df
        meta_sub = meta_df

    # Check we have samples in both conditions
    condition_counts = meta_sub[design_col].value_counts()
    if reference not in condition_counts or target not in condition_counts:
        logger.warning(f"Insufficient samples for contrast {target} vs {reference} "
                       f"in {cell_type}: {dict(condition_counts)}")
        return None

    if condition_counts[reference] < 2 or condition_counts[target] < 2:
        logger.warning(f"Need >=2 samples per group for {cell_type}: {dict(condition_counts)}")
        return None

    label = cell_type or "all"
    logger.info(f"Running DESeq2 for {label}: {target} vs {reference} "
                f"(n={condition_counts[target]} vs {condition_counts[reference]})")

    try:
        # Ensure counts are integers
        counts_int = counts_sub.astype(int)

        # Filter low-count genes
        gene_sums = counts_int.sum(axis=0)
        keep_genes = gene_sums >= 10
        counts_int = counts_int.loc[:, keep_genes]

        dds = DeseqDataSet(
            counts=counts_int,
            metadata=meta_sub,
            design_factors=design_col,
            refit_cooks=True,
        )
        dds.deseq2()

        stat_res = DeseqStats(dds, contrast=[design_col, target, reference])
        stat_res.summary()

        results = stat_res.results_df.copy()
        results["cell_type"] = label
        results["contrast"] = f"{target}_vs_{reference}"
        results["gene"] = results.index

        n_sig = (results["padj"] < 0.05).sum()
        n_up = ((results["padj"] < 0.05) & (results["log2FoldChange"] > 0)).sum()
        n_down = ((results["padj"] < 0.05) & (results["log2FoldChange"] < 0)).sum()
        logger.info(f"  {label}: {n_sig} DE genes (FDR<0.05): {n_up} up, {n_down} down")

        return results

    except Exception as e:
        logger.error(f"DESeq2 failed for {label}: {e}")
        return None


def run_scvi_de(adata: sc.AnnData, cfg: dict,
                groupby: str = "HIVstatus",
                group1: str = "positive",
                group2: str = "negative") -> pd.DataFrame | None:
    """
    Run scVI-based differential expression as a complementary approach.
    Requires a trained scVI model.
    """
    import scvi as scvi_module

    model_dir = Path(cfg["paths"]["models"]) / "scvi_model"
    if not model_dir.exists():
        logger.warning("scVI model not found, skipping scVI DE")
        return None

    try:
        model = scvi_module.model.SCVI.load(str(model_dir), adata=adata)
        de_results = model.differential_expression(
            groupby=groupby,
            group1=group1,
            group2=group2,
        )
        logger.info(f"scVI DE: {len(de_results)} genes tested")
        return de_results
    except Exception as e:
        logger.error(f"scVI DE failed: {e}")
        return None


def analyze_cell_proportions(adata: sc.AnnData, sample_col: str = "patientID",
                             cell_type_col: str = "cell_type",
                             condition_col: str = "HIVstatus") -> pd.DataFrame:
    """
    Compute cell type proportions per sample and test for differences.
    """
    from scipy.stats import mannwhitneyu

    # Compute proportions
    ct_counts = adata.obs.groupby([sample_col, cell_type_col]).size().unstack(fill_value=0)
    ct_props = ct_counts.div(ct_counts.sum(axis=1), axis=0)

    # Get condition per sample
    sample_condition = adata.obs.groupby(sample_col)[condition_col].first()

    results = []
    for ct in ct_props.columns:
        props = ct_props[ct]
        group1 = props[sample_condition == "positive"].values
        group2 = props[sample_condition == "negative"].values

        if len(group1) >= 2 and len(group2) >= 2:
            stat, pval = mannwhitneyu(group1, group2, alternative="two-sided")
        else:
            stat, pval = np.nan, np.nan

        results.append({
            "cell_type": ct,
            "mean_positive": np.mean(group1),
            "mean_negative": np.mean(group2),
            "log2fc": np.log2(np.mean(group1) / np.mean(group2))
            if np.mean(group2) > 0 else np.nan,
            "pval": pval,
            "n_positive": len(group1),
            "n_negative": len(group2),
        })

    results_df = pd.DataFrame(results)
    # Multiple testing correction
    from scipy.stats import false_discovery_control
    valid_pvals = results_df["pval"].dropna()
    if len(valid_pvals) > 0:
        results_df.loc[valid_pvals.index, "padj"] = false_discovery_control(valid_pvals.values)

    return results_df


def plot_volcano(results: pd.DataFrame, title: str, save_path: str | None = None) -> None:
    """Generate volcano plot from DE results."""
    fig, ax = plt.subplots(figsize=(8, 6))

    # Color by significance
    sig = (results["padj"] < 0.05) & (results["log2FoldChange"].abs() > 1)
    up = sig & (results["log2FoldChange"] > 0)
    down = sig & (results["log2FoldChange"] < 0)
    ns = ~sig

    neg_log10_p = -np.log10(results["padj"].clip(lower=1e-300))

    ax.scatter(results.loc[ns, "log2FoldChange"], neg_log10_p[ns],
               c="gray", alpha=0.3, s=5, label="NS")
    ax.scatter(results.loc[up, "log2FoldChange"], neg_log10_p[up],
               c="firebrick", alpha=0.6, s=8, label=f"Up ({up.sum()})")
    ax.scatter(results.loc[down, "log2FoldChange"], neg_log10_p[down],
               c="steelblue", alpha=0.6, s=8, label=f"Down ({down.sum()})")

    # Label top genes
    top_genes = results.nsmallest(10, "padj")
    for _, row in top_genes.iterrows():
        if row["padj"] < 0.05:
            ax.annotate(row["gene"], (row["log2FoldChange"], -np.log10(max(row["padj"], 1e-300))),
                        fontsize=6, alpha=0.8)

    ax.axhline(-np.log10(0.05), ls="--", color="gray", alpha=0.5)
    ax.axvline(-1, ls="--", color="gray", alpha=0.5)
    ax.axvline(1, ls="--", color="gray", alpha=0.5)
    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(adjusted p-value)")
    ax.set_title(title)
    ax.legend(loc="upper right")
    plt.tight_layout()

    if save_path:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, bbox_inches="tight", dpi=150)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Pseudobulk DE analysis")
    parser.add_argument("--config", default="pipeline/config.yaml")
    parser.add_argument("--arm", default="primary")
    parser.add_argument("--contrast", default="HIV_status",
                        help="Contrast name from config")
    parser.add_argument("--scvi-de", action="store_true",
                        help="Also run scVI-based DE")
    args = parser.parse_args()

    cfg = load_config(args.config)
    setup_paths(cfg)
    set_plotting_defaults()

    # Load annotated data
    ann_path = Path(cfg["paths"]["annotation_output"]) / f"annotated_{args.arm}.h5ad"
    logger.info(f"Loading annotated data from {ann_path}")
    adata = load_adata(str(ann_path))

    de_cfg = cfg["de"]
    de_output = Path(cfg["paths"]["de_output"])
    de_output.mkdir(parents=True, exist_ok=True)
    fig_dir = Path(cfg["paths"]["figures"]) / "de"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # Find contrast config
    contrast_cfg = None
    for c in de_cfg["contrasts"]:
        if c["name"] == args.contrast:
            contrast_cfg = c
            break
    if contrast_cfg is None:
        logger.error(f"Contrast '{args.contrast}' not found in config")
        return

    design_col = contrast_cfg["column"]
    reference = contrast_cfg["reference"]
    target = contrast_cfg["target"]

    # Step 1: Pseudobulk aggregation
    try:
        pb_result = compute_pseudobulk(
            adata,
            sample_col=de_cfg["sample_col"],
            groups_col="cell_type",
            min_cells=de_cfg["min_cells_per_pseudobulk"],
            min_counts=de_cfg["min_counts_per_pseudobulk"],
        )
        # Handle decoupler AnnData or manual (counts_df, meta_df) return
        if isinstance(pb_result, tuple):
            counts_df, meta_df = pb_result
        else:
            # decoupler returns AnnData
            if issparse(pb_result.X):
                counts_df = pd.DataFrame(
                    pb_result.X.toarray(), index=pb_result.obs_names,
                    columns=pb_result.var_names
                )
            else:
                counts_df = pd.DataFrame(
                    pb_result.X, index=pb_result.obs_names,
                    columns=pb_result.var_names
                )
            meta_df = pb_result.obs
    except Exception as e:
        logger.error(f"Pseudobulk failed: {e}")
        return

    # Step 2: DE per cell type
    all_de_results = []
    cell_types = meta_df["cell_type"].unique() if "cell_type" in meta_df.columns else ["all"]

    for ct in cell_types:
        results = run_deseq2(
            counts_df, meta_df,
            design_col=design_col,
            reference=reference,
            target=target,
            cell_type=ct,
        )
        if results is not None:
            all_de_results.append(results)

            # Volcano plot
            plot_volcano(
                results,
                title=f"{ct}: {target} vs {reference}",
                save_path=str(fig_dir / f"volcano_{ct.replace(' ', '_')}_{args.contrast}.png"),
            )

    # Combine and save all DE results
    if all_de_results:
        combined = pd.concat(all_de_results, ignore_index=True)
        combined_path = de_output / f"de_results_{args.contrast}_{args.arm}.csv"
        combined.to_csv(combined_path, index=False)
        logger.info(f"Combined DE results: {combined_path}")

        # Summary
        sig = combined[combined["padj"] < 0.05]
        logger.info(f"\nDE Summary ({args.contrast}):")
        for ct in sig["cell_type"].unique():
            ct_sig = sig[sig["cell_type"] == ct]
            logger.info(f"  {ct}: {len(ct_sig)} DE genes")

    # Step 3: scVI DE (optional)
    if args.scvi_de:
        scvi_results = run_scvi_de(adata, cfg, groupby=design_col,
                                    group1=target, group2=reference)
        if scvi_results is not None:
            scvi_path = de_output / f"scvi_de_{args.contrast}_{args.arm}.csv"
            scvi_results.to_csv(scvi_path)
            logger.info(f"scVI DE results: {scvi_path}")

    # Step 4: Cell type proportions
    logger.info("\nAnalyzing cell type proportions...")
    prop_results = analyze_cell_proportions(
        adata, sample_col=de_cfg["sample_col"],
        condition_col=design_col,
    )
    prop_path = de_output / f"cell_proportions_{args.contrast}_{args.arm}.csv"
    prop_results.to_csv(prop_path, index=False)
    logger.info(f"Cell proportion results: {prop_path}")

    sig_props = prop_results[prop_results.get("padj", pd.Series(dtype=float)) < 0.05]
    if len(sig_props) > 0:
        logger.info("Significant cell type proportion differences:")
        for _, row in sig_props.iterrows():
            logger.info(f"  {row['cell_type']}: log2FC={row['log2fc']:.2f}, padj={row['padj']:.4f}")

    logger.info("\nDone.")


if __name__ == "__main__":
    main()
