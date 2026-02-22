#!/usr/bin/env python3
"""
Phase 6: VDJ (TCR) integration with scirpy.

Steps:
  1. Read VDJ contigs per sample (Cell Ranger filtered_contig_annotations.csv)
  2. Merge VDJ data into T cell AnnData
  3. Define clonotypes by CDR3 amino acid identity
  4. Cross-reference doublets (cells with >1 productive TRB)
  5. Clonal expansion analysis
  6. Diversity metrics (Shannon, Simpson) per sample
  7. V(D)J gene usage analysis
  8. Clonotype network for cross-patient sharing

Usage:
    python pipeline/06_vdj_integration.py [--config pipeline/config.yaml]
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scirpy as ir
import matplotlib.pyplot as plt

from utils import (
    load_config, load_metadata, get_vdj_samples, setup_logging,
    setup_paths, save_adata, load_adata, set_plotting_defaults,
)

logger = setup_logging("06_vdj_integration")


def load_vdj_data(meta: pd.DataFrame, cfg: dict) -> dict:
    """
    Load VDJ contig annotations per sample.
    Returns dict of {sample_name: AnnData with IR data}.
    """
    vdj_meta = get_vdj_samples(meta)
    # Exclude missing batch
    vdj_meta = vdj_meta[vdj_meta["Batch"] != cfg["study"]["missing_batch"]]
    sample_names = vdj_meta["sampleName"].unique()

    vdj_data = {}
    for sample_name in sample_names:
        cr_base = Path(cfg["paths"]["cellranger_output"]) / sample_name / "outs"
        contigs_path = (cr_base / "per_sample_outs" / sample_name /
                        "vdj_t" / "filtered_contig_annotations.csv")

        # Try alternative paths
        if not contigs_path.exists():
            contigs_path = cr_base / "vdj_t" / "filtered_contig_annotations.csv"
        if not contigs_path.exists():
            contigs_path = (Path(cfg["paths"]["cellranger_output"]) /
                            sample_name / "vdj_t" / "filtered_contig_annotations.csv")

        if contigs_path.exists():
            try:
                adata_vdj = ir.io.read_10x_vdj(str(contigs_path))
                adata_vdj.obs["sample_id"] = sample_name
                vdj_data[sample_name] = adata_vdj
                logger.info(f"  Loaded VDJ for {sample_name}: "
                            f"{adata_vdj.n_obs} cells with contigs")
            except Exception as e:
                logger.warning(f"  Failed to load VDJ for {sample_name}: {e}")
        else:
            logger.warning(f"  VDJ contigs not found for {sample_name}: {contigs_path}")

    logger.info(f"Loaded VDJ data for {len(vdj_data)}/{len(sample_names)} samples")
    return vdj_data


def merge_vdj_with_gex(adata_t: sc.AnnData, vdj_data: dict) -> sc.AnnData:
    """
    Merge VDJ data into the T cell AnnData.
    Uses scirpy's merge_with_ir to add IR columns to .obs.
    """
    # Concatenate all VDJ anndata objects
    if not vdj_data:
        logger.error("No VDJ data to merge")
        return adata_t

    vdj_list = list(vdj_data.values())
    vdj_combined = vdj_list[0].concatenate(vdj_list[1:]) if len(vdj_list) > 1 else vdj_list[0]

    logger.info(f"Merging VDJ ({vdj_combined.n_obs} cells) with GEX ({adata_t.n_obs} cells)")

    ir.pp.merge_with_ir(adata_t, vdj_combined)

    # Count cells with TCR data
    has_tcr = adata_t.obs["IR_VDJ_1_junction_aa"].notna()
    n_tcr = has_tcr.sum()
    logger.info(f"Cells with TCR data: {n_tcr} ({n_tcr/adata_t.n_obs*100:.1f}%)")

    return adata_t


def flag_multi_chain_doublets(adata_t: sc.AnnData) -> sc.AnnData:
    """
    Flag cells with >1 productive TRB chain as likely doublets.
    Cross-reference with Scrublet doublet scores.
    """
    # Check for dual IR chains
    has_dual_trb = adata_t.obs["IR_VDJ_2_junction_aa"].notna()
    adata_t.obs["multi_trb"] = has_dual_trb

    n_dual = has_dual_trb.sum()
    logger.info(f"Cells with >1 TRB chain (likely doublets): {n_dual}")

    # Cross-reference with Scrublet
    if "doublet_score" in adata_t.obs.columns:
        dual_scores = adata_t.obs.loc[has_dual_trb, "doublet_score"]
        if len(dual_scores) > 0:
            logger.info(f"  Mean doublet score for multi-TRB cells: {dual_scores.mean():.3f}")
            logger.info(f"  Mean doublet score for single-TRB cells: "
                        f"{adata_t.obs.loc[~has_dual_trb, 'doublet_score'].mean():.3f}")

    return adata_t


def define_clonotypes(adata_t: sc.AnnData, cfg: dict) -> sc.AnnData:
    """
    Define clonotypes by CDR3 amino acid identity.
    Uses strict identity matching (CDR3α + CDR3β).
    """
    vdj_cfg = cfg["vdj"]

    logger.info(f"Defining clonotypes (metric={vdj_cfg['metric']}, "
                f"sequence={vdj_cfg['sequence']})...")

    ir.pp.ir_dist(adata_t, metric=vdj_cfg["metric"], sequence=vdj_cfg["sequence"])
    ir.tl.define_clonotypes(
        adata_t,
        receptor_arms=vdj_cfg["receptor_arms"],
        dual_ir=vdj_cfg["dual_ir"],
    )

    n_clonotypes = adata_t.obs["clone_id"].nunique()
    logger.info(f"Defined {n_clonotypes} unique clonotypes")

    return adata_t


def analyze_clonal_expansion(adata_t: sc.AnnData, cfg: dict) -> sc.AnnData:
    """Compute clonal expansion categories and diversity metrics."""
    vdj_cfg = cfg["vdj"]

    # Expansion categories
    ir.tl.clonal_expansion(
        adata_t,
        clip_at=vdj_cfg["expansion_clip_at"],
        expanded_in="sample_id",
    )

    logger.info("Clonal expansion distribution:")
    if "clonal_expansion" in adata_t.obs.columns:
        for cat, count in adata_t.obs["clonal_expansion"].value_counts().items():
            logger.info(f"  {cat}: {count}")

    # Diversity per sample
    for metric in vdj_cfg["diversity_metrics"]:
        try:
            ir.tl.alpha_diversity(adata_t, metric=metric, groupby="sample_id")
            logger.info(f"Computed {metric} diversity per sample")
        except Exception as e:
            logger.warning(f"Failed to compute {metric} diversity: {e}")

    return adata_t


def analyze_vgene_usage(adata_t: sc.AnnData) -> pd.DataFrame | None:
    """Analyze V gene usage across T cell subtypes."""
    try:
        vgene_col = "IR_VDJ_1_v_call"
        if vgene_col not in adata_t.obs.columns:
            logger.warning("V gene column not found")
            return None

        subtype_col = "tcell_subtype" if "tcell_subtype" in adata_t.obs.columns else "cell_type"

        result = ir.tl.group_abundance(
            adata_t,
            groupby=vgene_col,
            target_col=subtype_col,
        )
        logger.info(f"V gene usage analysis: {len(result)} V genes across subtypes")
        return result
    except Exception as e:
        logger.warning(f"V gene usage analysis failed: {e}")
        return None


def plot_vdj_results(adata_t: sc.AnnData, cfg: dict, save_dir: str) -> None:
    """Generate VDJ analysis plots."""
    save_path = Path(save_dir)
    save_path.mkdir(parents=True, exist_ok=True)

    # Clonal expansion on UMAP
    if "clonal_expansion" in adata_t.obs.columns:
        fig, axes = plt.subplots(1, 2, figsize=(16, 6))
        sc.pl.umap(adata_t, color="clonal_expansion", ax=axes[0], show=False,
                   title="Clonal expansion")
        if "tcell_subtype" in adata_t.obs.columns:
            sc.pl.umap(adata_t, color="tcell_subtype", ax=axes[1], show=False,
                       title="T cell subtype")
        plt.tight_layout()
        fig.savefig(save_path / "vdj_umap_expansion.png", bbox_inches="tight", dpi=150)
        plt.close(fig)

    # Clonal expansion by HIV status
    if "clonal_expansion" in adata_t.obs.columns and "HIVstatus" in adata_t.obs.columns:
        try:
            fig, ax = plt.subplots(figsize=(8, 5))
            ct = pd.crosstab(
                adata_t.obs["HIVstatus"],
                adata_t.obs["clonal_expansion"],
                normalize="index",
            )
            ct.plot(kind="bar", stacked=True, ax=ax)
            ax.set_ylabel("Fraction of cells")
            ax.set_title("Clonal expansion by HIV status")
            ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left")
            plt.tight_layout()
            fig.savefig(save_path / "vdj_expansion_by_hiv.png",
                        bbox_inches="tight", dpi=150)
            plt.close(fig)
        except Exception as e:
            logger.warning(f"Could not generate expansion by HIV plot: {e}")

    # TCR chain pairing
    has_tra = adata_t.obs["IR_VJ_1_junction_aa"].notna()
    has_trb = adata_t.obs["IR_VDJ_1_junction_aa"].notna()
    has_both = has_tra & has_trb
    has_any = has_tra | has_trb

    logger.info(f"TCR chain pairing:")
    logger.info(f"  TRA only: {(has_tra & ~has_trb).sum()}")
    logger.info(f"  TRB only: {(has_trb & ~has_tra).sum()}")
    logger.info(f"  Paired (TRA+TRB): {has_both.sum()}")
    logger.info(f"  Any chain: {has_any.sum()}")
    logger.info(f"  Pairing rate: {has_both.sum()/has_any.sum()*100:.1f}% (of cells with any chain)")

    # Diversity per sample barplot
    diversity_col = None
    for col in adata_t.obs.columns:
        if "shannon" in col.lower() or "diversity" in col.lower():
            diversity_col = col
            break

    if diversity_col:
        try:
            div_per_sample = adata_t.obs.groupby("sample_id")[diversity_col].first().dropna()
            fig, ax = plt.subplots(figsize=(max(10, len(div_per_sample) * 0.4), 5))
            div_per_sample.sort_values().plot(kind="bar", ax=ax)
            ax.set_ylabel("Shannon diversity")
            ax.set_title("TCR diversity per sample")
            plt.xticks(rotation=90)
            plt.tight_layout()
            fig.savefig(save_path / "vdj_diversity_per_sample.png",
                        bbox_inches="tight", dpi=150)
            plt.close(fig)
        except Exception as e:
            logger.warning(f"Diversity plot failed: {e}")


def main():
    parser = argparse.ArgumentParser(description="VDJ integration with scirpy")
    parser.add_argument("--config", default="pipeline/config.yaml")
    parser.add_argument("--arm", default="primary")
    args = parser.parse_args()

    cfg = load_config(args.config)
    setup_paths(cfg)
    set_plotting_defaults()

    meta = load_metadata(cfg)

    # Load T cell data
    t_path = Path(cfg["paths"]["tcell_output"]) / f"tcell_{args.arm}.h5ad"
    logger.info(f"Loading T cell data from {t_path}")
    adata_t = load_adata(str(t_path))

    # Step 1: Load VDJ data
    logger.info("Loading VDJ data per sample...")
    vdj_data = load_vdj_data(meta, cfg)

    if not vdj_data:
        logger.error("No VDJ data loaded. Check Cell Ranger output paths.")
        return

    # Step 2: Merge VDJ with GEX
    adata_t = merge_vdj_with_gex(adata_t, vdj_data)

    # Step 3: Flag multi-chain doublets
    adata_t = flag_multi_chain_doublets(adata_t)

    # Step 4: Define clonotypes
    adata_t = define_clonotypes(adata_t, cfg)

    # Step 5: Clonal expansion + diversity
    adata_t = analyze_clonal_expansion(adata_t, cfg)

    # Step 6: V gene usage
    vgene_df = analyze_vgene_usage(adata_t)
    if vgene_df is not None:
        vgene_path = Path(cfg["paths"]["vdj_output"]) / f"vgene_usage_{args.arm}.csv"
        vgene_path.parent.mkdir(parents=True, exist_ok=True)
        vgene_df.to_csv(vgene_path, index=False)

    # Step 7: Plots
    fig_dir = Path(cfg["paths"]["figures"]) / "vdj"
    plot_vdj_results(adata_t, cfg, str(fig_dir))

    # Save
    out_path = Path(cfg["paths"]["vdj_output"]) / f"tcell_vdj_{args.arm}.h5ad"
    save_adata(adata_t, str(out_path))

    logger.info(f"\nDone. VDJ-integrated T cell data: {out_path}")


if __name__ == "__main__":
    main()
