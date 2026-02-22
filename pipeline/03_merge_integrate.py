#!/usr/bin/env python3
"""
Phase 3: Merge QC-filtered libraries + scVI integration.

Steps:
  1. Concatenate per-library AnnData objects (inner join on genes)
  2. Add sample/patient/batch metadata
  3. Store raw counts in adata.layers["counts"]
  4. Normalize + log1p (for visualization/HVG only)
  5. Batch-aware HVG selection (seurat_v3 flavor)
  6. scVI integration (GPU) with batch_key="Batch"
  7. Harmony fallback
  8. Neighbor graph + multi-resolution Leiden + UMAP
  9. Integration quality assessment with scib

Usage:
    python pipeline/03_merge_integrate.py [--config pipeline/config.yaml]
    python pipeline/03_merge_integrate.py --arm primary   # Default: unsorted only
    python pipeline/03_merge_integrate.py --arm focused   # Sorted only
    python pipeline/03_merge_integrate.py --harmony       # Run Harmony as well
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.sparse import issparse

from utils import (
    load_config, load_metadata, get_gex_samples, filter_samples_by_arm,
    map_clinical_metadata, load_clinical_metadata, setup_logging,
    setup_paths, save_adata, load_adata, set_plotting_defaults,
    plot_umap_batch,
)

logger = setup_logging("03_merge_integrate")


def load_qc_libraries(sample_names: list[str], cfg: dict) -> list[sc.AnnData]:
    """Load QC-filtered h5ad files for each sample."""
    qc_dir = Path(cfg["paths"]["qc_output"])
    adata_list = []
    skipped = []

    for sample_name in sample_names:
        path = qc_dir / f"{sample_name}_qc.h5ad"
        if not path.exists():
            logger.warning(f"QC output not found, skipping: {path}")
            skipped.append(sample_name)
            continue
        adata = sc.read_h5ad(str(path))
        adata.obs["sample_id"] = sample_name
        adata_list.append(adata)
        logger.info(f"  Loaded {sample_name}: {adata.n_obs} cells x {adata.n_vars} genes")

    if skipped:
        logger.warning(f"Skipped {len(skipped)} samples with missing QC output")

    return adata_list


def merge_libraries(
    adata_list: list[sc.AnnData],
    meta: pd.DataFrame,
    clinical: pd.DataFrame | None = None,
) -> sc.AnnData:
    """
    Concatenate libraries and add metadata.
    Uses inner join to keep only genes present in ALL libraries.
    """
    logger.info(f"Merging {len(adata_list)} libraries (inner join on genes)...")

    adata = ad.concat(adata_list, join="inner", label="sample_id",
                      keys=[a.obs["sample_id"].iloc[0] for a in adata_list],
                      index_unique="-")
    adata.obs_names_make_unique()

    logger.info(f"Merged: {adata.n_obs} cells x {adata.n_vars} genes")

    # Add metadata from the sample table
    # Build a sample-level lookup
    gex_meta = meta[meta["locus"] == "5primeGEX"].drop_duplicates(subset="sampleName")
    sample_lookup = gex_meta.set_index("sampleName")[
        ["patientID", "Batch", "sortingCT", "tissueType", "HIVstatus"]
    ].to_dict("index")

    for col in ["patientID", "Batch", "sortingCT", "tissueType", "HIVstatus"]:
        adata.obs[col] = adata.obs["sample_id"].map(
            lambda s, c=col: sample_lookup.get(s, {}).get(c, "unknown")
        )

    # Add clinical metadata if available
    if clinical is not None:
        obs = adata.obs.copy()
        obs["ptid"] = obs["patientID"].str.extract(r"008_(\d+)").astype(float)
        clinical_c = clinical.copy()
        clinical_c["ptid"] = clinical_c["ptid"].astype(float)
        obs = obs.merge(clinical_c, on="ptid", how="left", suffixes=("", "_clin"))
        obs.index = adata.obs.index
        # Transfer clinical columns
        for col in ["best_resp", "bl_plasma_grp", "best_resp_new"]:
            if col in obs.columns:
                adata.obs[col] = obs[col].values

    # Store raw counts
    adata.layers["counts"] = adata.X.copy()
    logger.info("Stored raw counts in adata.layers['counts']")

    return adata


def normalize_and_hvg(adata: sc.AnnData, cfg: dict) -> sc.AnnData:
    """
    Normalize, log-transform, and select HVGs.
    scVI will use raw counts from layers["counts"], not the normalized data.
    """
    int_cfg = cfg["integration"]

    # Normalize for visualization and HVG selection
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # Store normalized in .raw for plotting

    # Batch-aware HVG selection using raw counts
    logger.info(f"Selecting {int_cfg['n_top_genes']} HVGs "
                f"(flavor={int_cfg['hvg_flavor']}, batch-aware)...")
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=int_cfg["n_top_genes"],
        flavor=int_cfg["hvg_flavor"],
        batch_key="sample_id",
        layer="counts",
    )

    # Exclude mito/ribo/hb genes from HVGs
    if "mt" in adata.var.columns:
        adata.var.loc[adata.var["mt"], "highly_variable"] = False
    if "ribo" in adata.var.columns:
        adata.var.loc[adata.var["ribo"], "highly_variable"] = False
    if "hb" in adata.var.columns:
        adata.var.loc[adata.var["hb"], "highly_variable"] = False

    n_hvg = adata.var["highly_variable"].sum()
    logger.info(f"Selected {n_hvg} HVGs (after excluding MT/ribo/HB genes)")

    return adata


def run_scvi(adata: sc.AnnData, cfg: dict) -> sc.AnnData:
    """Run scVI integration. Requires GPU."""
    import scvi as scvi_module

    scvi_cfg = cfg["integration"]["scvi"]
    logger.info("Setting up scVI model...")

    scvi_module.model.SCVI.setup_anndata(
        adata,
        layer="counts",
        batch_key=scvi_cfg["batch_key"],
        continuous_covariate_keys=scvi_cfg["continuous_covariates"],
    )

    model = scvi_module.model.SCVI(
        adata,
        n_latent=scvi_cfg["n_latent"],
        n_layers=scvi_cfg["n_layers"],
    )

    logger.info(f"Training scVI (max_epochs={scvi_cfg['max_epochs']}, "
                f"early_stopping={scvi_cfg['early_stopping']})...")
    model.train(
        max_epochs=scvi_cfg["max_epochs"],
        early_stopping=scvi_cfg["early_stopping"],
    )

    # Get latent representation
    adata.obsm["X_scVI"] = model.get_latent_representation()
    logger.info(f"scVI latent space: {adata.obsm['X_scVI'].shape}")

    # Save trained model
    model_dir = Path(cfg["paths"]["models"]) / "scvi_model"
    model_dir.mkdir(parents=True, exist_ok=True)
    model.save(str(model_dir), overwrite=True)
    logger.info(f"Saved scVI model to {model_dir}")

    return adata


def run_harmony(adata: sc.AnnData, cfg: dict) -> sc.AnnData:
    """Run Harmony integration as fallback/comparison."""
    import harmonypy

    h_cfg = cfg["integration"]["harmony"]
    logger.info("Running Harmony integration...")

    # PCA first
    sc.pp.pca(adata, n_comps=h_cfg["n_pcs"])

    harmony_out = harmonypy.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        h_cfg["batch_key"],
    )
    adata.obsm["X_harmony"] = harmony_out.Z_corr.T

    logger.info(f"Harmony embedding: {adata.obsm['X_harmony'].shape}")
    return adata


def cluster_and_umap(adata: sc.AnnData, cfg: dict, use_rep: str = "X_scVI") -> sc.AnnData:
    """Compute neighbor graph, multi-resolution Leiden, and UMAP."""
    int_cfg = cfg["integration"]

    logger.info(f"Computing neighbors (use_rep={use_rep}, n_neighbors={int_cfg['n_neighbors']})...")
    sc.pp.neighbors(adata, use_rep=use_rep, n_neighbors=int_cfg["n_neighbors"])

    # Multi-resolution Leiden clustering
    for res in int_cfg["leiden_resolutions"]:
        key = f"leiden_{res}"
        sc.tl.leiden(adata, resolution=res, key_added=key)
        n_clusters = adata.obs[key].nunique()
        logger.info(f"  Leiden res={res}: {n_clusters} clusters")

    # UMAP
    logger.info(f"Computing UMAP (min_dist={int_cfg['umap_min_dist']})...")
    sc.tl.umap(adata, min_dist=int_cfg["umap_min_dist"])

    return adata


def assess_integration(adata: sc.AnnData, cfg: dict, rep_key: str = "X_scVI") -> dict:
    """
    Assess integration quality using scib metrics.
    Returns dict of metric name -> value.
    """
    try:
        import scib

        logger.info("Computing integration quality metrics (scib)...")
        batch_key = cfg["integration"]["scvi"]["batch_key"]
        label_key = f"leiden_{cfg['integration']['default_resolution']}"

        metrics = {}

        # Batch mixing metrics
        try:
            metrics["batch_ASW"] = scib.metrics.silhouette_batch(
                adata, batch_key=batch_key, label_key=label_key, embed=rep_key
            )
        except Exception as e:
            logger.warning(f"batch_ASW failed: {e}")

        try:
            metrics["graph_connectivity"] = scib.metrics.graph_connectivity(
                adata, label_key=label_key
            )
        except Exception as e:
            logger.warning(f"graph_connectivity failed: {e}")

        # Bio conservation
        try:
            metrics["cell_type_ASW"] = scib.metrics.silhouette(
                adata, label_key=label_key, embed=rep_key
            )
        except Exception as e:
            logger.warning(f"cell_type_ASW failed: {e}")

        try:
            metrics["NMI"] = scib.metrics.nmi(adata, label_key, label_key)
        except Exception as e:
            logger.warning(f"NMI failed: {e}")

        logger.info("Integration metrics:")
        for name, val in metrics.items():
            logger.info(f"  {name}: {val:.4f}")

        return metrics

    except ImportError:
        logger.warning("scib not installed, skipping integration assessment")
        return {}


def main():
    parser = argparse.ArgumentParser(description="Merge + scVI integration")
    parser.add_argument("--config", default="pipeline/config.yaml")
    parser.add_argument("--arm", default="primary",
                        choices=["primary", "focused"],
                        help="Analysis arm (default: primary = unsorted only)")
    parser.add_argument("--harmony", action="store_true",
                        help="Also run Harmony for comparison")
    parser.add_argument("--skip-scvi", action="store_true",
                        help="Skip scVI (e.g., if no GPU)")
    args = parser.parse_args()

    cfg = load_config(args.config)
    setup_paths(cfg)
    set_plotting_defaults()

    meta = load_metadata(cfg)

    # Filter to analysis arm
    arm_meta = filter_samples_by_arm(get_gex_samples(meta), args.arm, cfg)
    # Exclude missing batch
    arm_meta = arm_meta[arm_meta["Batch"] != cfg["study"]["missing_batch"]]
    sample_names = arm_meta["sampleName"].unique()
    logger.info(f"Analysis arm: {args.arm} ({len(sample_names)} samples)")

    # Load clinical metadata
    clinical = None
    try:
        clinical = load_clinical_metadata(cfg)
    except Exception as e:
        logger.warning(f"Could not load clinical metadata: {e}")

    # Step 1: Load QC-filtered libraries
    logger.info("Loading QC-filtered libraries...")
    adata_list = load_qc_libraries(sample_names.tolist(), cfg)
    if not adata_list:
        logger.error("No libraries loaded. Check QC output directory.")
        return

    # Step 2: Merge
    adata = merge_libraries(adata_list, meta, clinical)

    # Step 3: Normalize + HVG
    adata = normalize_and_hvg(adata, cfg)

    # Step 4: Integration
    use_rep = "X_scVI"

    if not args.skip_scvi:
        try:
            adata = run_scvi(adata, cfg)
        except Exception as e:
            logger.error(f"scVI failed: {e}")
            logger.info("Falling back to Harmony...")
            args.harmony = True
            args.skip_scvi = True

    if args.harmony or args.skip_scvi:
        adata = run_harmony(adata, cfg)
        if args.skip_scvi:
            use_rep = "X_harmony"

    # Step 5: Clustering + UMAP
    adata = cluster_and_umap(adata, cfg, use_rep=use_rep)

    # Step 6: Integration assessment
    metrics = assess_integration(adata, cfg, rep_key=use_rep)
    if metrics:
        metrics_df = pd.DataFrame([metrics])
        metrics_path = Path(cfg["paths"]["integration_output"]) / f"integration_metrics_{args.arm}.csv"
        metrics_df.to_csv(metrics_path, index=False)

    # Save
    out_dir = Path(cfg["paths"]["integration_output"])
    out_path = out_dir / f"integrated_{args.arm}.h5ad"
    save_adata(adata, str(out_path))

    # Generate plots
    fig_dir = Path(cfg["paths"]["figures"]) / "integration"
    fig_dir.mkdir(parents=True, exist_ok=True)
    plot_umap_batch(adata, save_path=str(fig_dir / f"umap_batch_{args.arm}.png"))

    logger.info(f"\nDone. Output: {out_path}")
    logger.info(f"Final dataset: {adata.n_obs} cells x {adata.n_vars} genes")


if __name__ == "__main__":
    main()
