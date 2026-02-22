#!/usr/bin/env python3
"""
Phase 5: T cell subclustering with re-integration.

Steps:
  1. Subset T cells from annotated data (CD4 T, CD8 T, T cells)
  2. Re-select HVGs on T cell subset (captures T-specific variation)
  3. Re-integrate with scVI on T cell subset
  4. Leiden clustering at higher resolution
  5. Annotate T cell subtypes with canonical markers:
     - CD8: Naive, Effector, Memory, Exhausted, TPEX
     - CD4: Naive, Th1, Th2, Treg, Tfh

Usage:
    python pipeline/05_tcell_subcluster.py [--config pipeline/config.yaml]
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

logger = setup_logging("05_tcell_subcluster")


def subset_tcells(adata: sc.AnnData) -> sc.AnnData:
    """Extract T cells based on cell_type annotation."""
    t_cell_types = ["CD4 T", "CD8 T", "T cells", "NK/T"]
    mask = adata.obs["cell_type"].isin(t_cell_types)
    adata_t = adata[mask].copy()
    logger.info(f"Subset T cells: {adata_t.n_obs} cells "
                f"({adata_t.n_obs/adata.n_obs*100:.1f}% of total)")

    # Distribution
    for ct, count in adata_t.obs["cell_type"].value_counts().items():
        logger.info(f"  {ct}: {count}")

    return adata_t


def reprocess_tcells(adata_t: sc.AnnData, cfg: dict) -> sc.AnnData:
    """Re-run HVG selection and scVI integration on T cell subset."""
    tcell_cfg = cfg["tcell"]

    # HVG selection on T cells using raw counts
    logger.info(f"Selecting {tcell_cfg['n_top_genes']} HVGs for T cell subset...")
    sc.pp.highly_variable_genes(
        adata_t,
        n_top_genes=tcell_cfg["n_top_genes"],
        flavor="seurat_v3",
        batch_key="sample_id",
        layer="counts",
    )

    # Exclude MT/ribo/HB from HVGs
    for col in ["mt", "ribo", "hb"]:
        if col in adata_t.var.columns:
            adata_t.var.loc[adata_t.var[col], "highly_variable"] = False

    n_hvg = adata_t.var["highly_variable"].sum()
    logger.info(f"  {n_hvg} HVGs selected")

    return adata_t


def integrate_tcells_scvi(adata_t: sc.AnnData, cfg: dict) -> sc.AnnData:
    """Run scVI integration on T cell subset."""
    import scvi as scvi_module

    tcell_cfg = cfg["tcell"]["scvi"]

    scvi_module.model.SCVI.setup_anndata(
        adata_t,
        layer="counts",
        batch_key=tcell_cfg["batch_key"],
    )

    model = scvi_module.model.SCVI(
        adata_t,
        n_latent=tcell_cfg["n_latent"],
    )

    logger.info(f"Training T cell scVI (max_epochs={tcell_cfg['max_epochs']})...")
    model.train(
        max_epochs=tcell_cfg["max_epochs"],
        early_stopping=True,
    )

    adata_t.obsm["X_scVI"] = model.get_latent_representation()

    # Save model
    model_dir = Path(cfg["paths"]["models"]) / "scvi_tcell_model"
    model_dir.mkdir(parents=True, exist_ok=True)
    model.save(str(model_dir), overwrite=True)
    logger.info(f"Saved T cell scVI model to {model_dir}")

    return adata_t


def cluster_tcells(adata_t: sc.AnnData, cfg: dict) -> sc.AnnData:
    """Cluster T cells and compute UMAP."""
    tcell_cfg = cfg["tcell"]

    sc.pp.neighbors(adata_t, use_rep="X_scVI",
                    n_neighbors=tcell_cfg["n_neighbors"])

    # Multiple resolutions
    for res in [0.5, 0.8, 1.0, 1.5, 2.0]:
        key = f"leiden_{res}"
        sc.tl.leiden(adata_t, resolution=res, key_added=key)
        n_clusters = adata_t.obs[key].nunique()
        logger.info(f"  T cell Leiden res={res}: {n_clusters} clusters")

    sc.tl.umap(adata_t, min_dist=0.3)

    return adata_t


def score_tcell_subtypes(adata_t: sc.AnnData, cfg: dict) -> sc.AnnData:
    """
    Score T cell subtypes using marker gene sets from config.
    Uses scanpy's score_genes to assign continuous scores.
    """
    for lineage in ["cd8_markers", "cd4_markers"]:
        markers = cfg["tcell"].get(lineage, {})
        for subtype, genes in markers.items():
            available = [g for g in genes if g in adata_t.var_names or
                         (adata_t.raw is not None and g in adata_t.raw.var_names)]
            if not available:
                continue
            score_name = f"score_{lineage.split('_')[0]}_{subtype}"
            try:
                sc.tl.score_genes(adata_t, gene_list=available, score_name=score_name)
                logger.info(f"  Scored {score_name}: {len(available)} genes")
            except Exception as e:
                logger.warning(f"  Failed to score {score_name}: {e}")

    return adata_t


def annotate_tcell_subtypes(adata_t: sc.AnnData, cluster_key: str = "leiden_1.0") -> sc.AnnData:
    """
    Annotate T cell subtypes based on marker scores and cluster expression.
    This provides initial labels; manual refinement is expected.
    """
    # Use rank_genes_groups to find cluster markers
    sc.tl.rank_genes_groups(adata_t, groupby=cluster_key, method="wilcoxon",
                            use_raw=True, pts=True)

    # Key discriminating markers
    cd8_markers = {"CD8A", "CD8B"}
    cd4_markers = {"CD4", "IL7R"}
    naive_markers = {"CCR7", "SELL", "LEF1", "TCF7"}
    effector_markers = {"GZMB", "PRF1", "GNLY", "GZMA"}
    exhausted_markers = {"PDCD1", "LAG3", "HAVCR2", "TOX", "TIGIT"}
    treg_markers = {"FOXP3", "IL2RA", "CTLA4"}
    proliferating_markers = {"MKI67", "TOP2A"}

    # Score each cluster for these marker sets
    cluster_labels = {}
    clusters = adata_t.obs[cluster_key].unique()

    for cl in clusters:
        mask = adata_t.obs[cluster_key] == cl
        subset = adata_t[mask]

        def _mean_expr(genes):
            available = [g for g in genes if g in adata_t.var_names]
            if not available:
                return 0
            from scipy.sparse import issparse
            expr = subset[:, available].X
            if issparse(expr):
                expr = expr.toarray()
            return float(np.mean(expr > 0))  # Fraction of cells expressing

        cd8_score = _mean_expr(cd8_markers)
        cd4_score = _mean_expr(cd4_markers)
        naive_score = _mean_expr(naive_markers)
        effector_score = _mean_expr(effector_markers)
        exhausted_score = _mean_expr(exhausted_markers)
        treg_score = _mean_expr(treg_markers)
        prolif_score = _mean_expr(proliferating_markers)

        # Decision tree
        label = "T cell"
        if cd8_score > cd4_score:
            if prolif_score > 0.3:
                label = "CD8 Proliferating"
            elif exhausted_score > 0.3:
                label = "CD8 Exhausted"
            elif effector_score > 0.3:
                label = "CD8 Effector"
            elif naive_score > 0.3:
                label = "CD8 Naive"
            else:
                label = "CD8 Memory"
        else:
            if treg_score > 0.2:
                label = "CD4 Treg"
            elif prolif_score > 0.3:
                label = "CD4 Proliferating"
            elif naive_score > 0.3:
                label = "CD4 Naive"
            else:
                label = "CD4 Memory"

        cluster_labels[cl] = label
        logger.info(f"  Cluster {cl}: {label} "
                    f"(CD8={cd8_score:.2f}, CD4={cd4_score:.2f}, "
                    f"Naive={naive_score:.2f}, Eff={effector_score:.2f}, "
                    f"Exh={exhausted_score:.2f}, Treg={treg_score:.2f})")

    adata_t.obs["tcell_subtype"] = adata_t.obs[cluster_key].map(cluster_labels)

    return adata_t


def plot_tcell_results(adata_t: sc.AnnData, cfg: dict, save_dir: str) -> None:
    """Generate T cell subclustering plots."""
    save_path = Path(save_dir)
    save_path.mkdir(parents=True, exist_ok=True)

    # UMAP by subtype
    fig, axes = plt.subplots(1, 3, figsize=(24, 7))
    sc.pl.umap(adata_t, color="tcell_subtype", ax=axes[0], show=False,
               title="T cell subtype")
    sc.pl.umap(adata_t, color="HIVstatus", ax=axes[1], show=False,
               title="HIV status")
    sc.pl.umap(adata_t, color="leiden_1.0", ax=axes[2], show=False,
               title="Leiden clusters")
    plt.tight_layout()
    fig.savefig(save_path / "tcell_umap_overview.png", bbox_inches="tight", dpi=150)
    plt.close(fig)

    # Dot plot of T cell markers
    all_markers = []
    for lineage in ["cd8_markers", "cd4_markers"]:
        for subtype, genes in cfg["tcell"].get(lineage, {}).items():
            for g in genes:
                if g in adata_t.var_names or (adata_t.raw is not None and g in adata_t.raw.var_names):
                    all_markers.append(g)
    all_markers = list(dict.fromkeys(all_markers))  # Deduplicate

    if all_markers:
        fig, ax = plt.subplots(figsize=(max(14, len(all_markers) * 0.6), 8))
        sc.pl.dotplot(adata_t, var_names=all_markers, groupby="tcell_subtype",
                      standard_scale="var", show=False, ax=ax)
        plt.tight_layout()
        fig.savefig(save_path / "tcell_dotplot_subtypes.png",
                    bbox_inches="tight", dpi=150)
        plt.close(fig)

    # Per-subtype marker genes
    sc.tl.rank_genes_groups(adata_t, groupby="tcell_subtype", method="wilcoxon",
                            use_raw=True)
    fig = sc.pl.rank_genes_groups_dotplot(
        adata_t, n_genes=5, show=False, return_fig=True
    )
    if fig is not None:
        fig.savefig(save_path / "tcell_de_markers.png", bbox_inches="tight", dpi=150)
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="T cell subclustering")
    parser.add_argument("--config", default="pipeline/config.yaml")
    parser.add_argument("--arm", default="primary")
    parser.add_argument("--skip-scvi", action="store_true",
                        help="Skip scVI re-integration (use PCA + Harmony instead)")
    args = parser.parse_args()

    cfg = load_config(args.config)
    setup_paths(cfg)
    set_plotting_defaults()

    # Load annotated data
    in_path = Path(cfg["paths"]["annotation_output"]) / f"annotated_{args.arm}.h5ad"
    logger.info(f"Loading annotated data from {in_path}")
    adata = load_adata(str(in_path))

    # Step 1: Subset T cells
    adata_t = subset_tcells(adata)
    del adata  # Free memory

    # Step 2: Re-process (HVG on T cell subset)
    adata_t = reprocess_tcells(adata_t, cfg)

    # Step 3: Re-integrate
    if not args.skip_scvi:
        try:
            adata_t = integrate_tcells_scvi(adata_t, cfg)
        except Exception as e:
            logger.error(f"scVI failed for T cells: {e}")
            logger.info("Falling back to PCA-based approach")
            sc.pp.pca(adata_t, n_comps=30)
            import harmonypy
            harmony_out = harmonypy.run_harmony(
                adata_t.obsm["X_pca"], adata_t.obs, "Batch"
            )
            adata_t.obsm["X_scVI"] = harmony_out.Z_corr.T
    else:
        sc.pp.pca(adata_t, n_comps=30)
        import harmonypy
        harmony_out = harmonypy.run_harmony(
            adata_t.obsm["X_pca"], adata_t.obs, "Batch"
        )
        adata_t.obsm["X_scVI"] = harmony_out.Z_corr.T

    # Step 4: Cluster + UMAP
    adata_t = cluster_tcells(adata_t, cfg)

    # Step 5: Score and annotate subtypes
    adata_t = score_tcell_subtypes(adata_t, cfg)
    adata_t = annotate_tcell_subtypes(adata_t, cluster_key=f"leiden_{cfg['tcell']['resolution']}")

    # Step 6: Plots
    fig_dir = Path(cfg["paths"]["figures"]) / "tcell"
    plot_tcell_results(adata_t, cfg, str(fig_dir))

    # Save
    out_path = Path(cfg["paths"]["tcell_output"]) / f"tcell_{args.arm}.h5ad"
    save_adata(adata_t, str(out_path))

    logger.info(f"\nDone. T cell data: {out_path}")
    logger.info(f"T cell subtypes:")
    for st, count in adata_t.obs["tcell_subtype"].value_counts().items():
        logger.info(f"  {st}: {count}")


if __name__ == "__main__":
    main()
