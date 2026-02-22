"""
Shared utility functions for the scRNA-seq + VDJ pipeline.
"""

import os
import logging
from pathlib import Path

import yaml
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

def load_config(config_path: str = "pipeline/config.yaml") -> dict:
    """Load pipeline configuration from YAML."""
    with open(config_path) as f:
        return yaml.safe_load(f)


def setup_paths(cfg: dict) -> None:
    """Create all output directories from config."""
    for key, path in cfg["paths"].items():
        if key in ("metadata", "clinical_metadata", "fastq_base",
                    "gex_reference", "vdj_reference"):
            continue
        Path(path).mkdir(parents=True, exist_ok=True)


def setup_logging(name: str, log_dir: str = "output/logs") -> logging.Logger:
    """Configure logging for a pipeline step."""
    Path(log_dir).mkdir(parents=True, exist_ok=True)
    log = logging.getLogger(name)
    log.setLevel(logging.INFO)
    if not log.handlers:
        fh = logging.FileHandler(f"{log_dir}/{name}.log")
        fh.setFormatter(logging.Formatter(
            "%(asctime)s [%(levelname)s] %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        ))
        log.addHandler(fh)
        sh = logging.StreamHandler()
        sh.setFormatter(logging.Formatter("[%(levelname)s] %(message)s"))
        log.addHandler(sh)
    return log

# ---------------------------------------------------------------------------
# Metadata
# ---------------------------------------------------------------------------

def load_metadata(cfg: dict) -> pd.DataFrame:
    """Load and clean sample metadata, applying known corrections."""
    meta = pd.read_csv(cfg["paths"]["metadata"])

    # Apply HIV status overrides (008_252 VDJ row mislabeled)
    for patient, status in cfg["study"].get("hiv_status_overrides", {}).items():
        mask = meta["patientID"] == patient
        n_corrected = mask.sum() - (meta.loc[mask, "HIVstatus"] == status).sum()
        if n_corrected > 0:
            logger.info(f"Correcting HIV status for {patient}: {n_corrected} rows -> {status}")
        meta.loc[mask, "HIVstatus"] = status

    return meta


def load_clinical_metadata(cfg: dict) -> pd.DataFrame:
    """Load clinical outcome metadata."""
    return pd.read_csv(cfg["paths"]["clinical_metadata"])


def get_gex_samples(meta: pd.DataFrame) -> pd.DataFrame:
    """Filter metadata to GEX (gene expression) rows only."""
    return meta[meta["locus"] == "5primeGEX"].copy()


def get_vdj_samples(meta: pd.DataFrame) -> pd.DataFrame:
    """Filter metadata to VDJ rows only."""
    return meta[meta["locus"] == "5primeVDJ"].copy()


def get_sample_arm(meta_row: pd.Series, cfg: dict) -> str:
    """Classify a sample into analysis arm: primary, focused, or validation."""
    sorting = meta_row["sortingCT"]
    patient = meta_row["patientID"]
    validation_patients = cfg["study"]["analysis_arms"]["validation"]["patients"]

    if sorting in ("CD4pTCells", "CD8pTCells"):
        if patient in validation_patients:
            return "validation"
        return "focused"
    return "primary"


def filter_samples_by_arm(meta: pd.DataFrame, arm: str, cfg: dict) -> pd.DataFrame:
    """Filter metadata to a specific analysis arm."""
    arms = cfg["study"]["analysis_arms"]
    if arm == "primary":
        return meta[meta["sortingCT"].isin(arms["primary"]["sorting_filter"])].copy()
    elif arm == "focused":
        return meta[meta["sortingCT"].isin(arms["focused"]["sorting_filter"])].copy()
    elif arm == "validation":
        return meta[meta["patientID"].isin(arms["validation"]["patients"])].copy()
    else:
        raise ValueError(f"Unknown analysis arm: {arm}")


def map_clinical_metadata(meta: pd.DataFrame, clinical: pd.DataFrame) -> pd.DataFrame:
    """
    Merge clinical metadata with sample metadata.
    Extracts numeric patient ID from '008_XXX' format to match clinical ptid.
    """
    meta = meta.copy()
    # Extract numeric portion: '008_216' -> 216
    meta["ptid"] = meta["patientID"].str.extract(r"008_(\d+)").astype(float)
    clinical = clinical.copy()
    clinical["ptid"] = clinical["ptid"].astype(float)
    merged = meta.merge(clinical, on="ptid", how="left")
    return merged

# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------

def find_cellranger_output(sample_name: str, cfg: dict) -> dict:
    """Locate Cell Ranger multi output directories for a sample."""
    base = Path(cfg["paths"]["cellranger_output"]) / sample_name / "outs"
    return {
        "raw_h5": base / "multi" / "count" / "raw_feature_bc_matrix.h5",
        "filtered_h5": base / "multi" / "count" / "filtered_feature_bc_matrix.h5",
        "metrics": base / "multi" / "count" / "summary.csv",
        "vdj_contigs": base / "per_sample_outs" / sample_name / "vdj_t" / "filtered_contig_annotations.csv",
        "vdj_airr": base / "per_sample_outs" / sample_name / "vdj_t" / "airr_rearrangement.tsv",
    }


def save_adata(adata, path: str, **kwargs) -> None:
    """Save AnnData with directory creation."""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path, **kwargs)
    logger.info(f"Saved AnnData ({adata.n_obs} cells x {adata.n_vars} genes) -> {path}")


def load_adata(path: str) -> "sc.AnnData":
    """Load AnnData from h5ad."""
    adata = sc.read_h5ad(path)
    logger.info(f"Loaded AnnData ({adata.n_obs} cells x {adata.n_vars} genes) <- {path}")
    return adata

# ---------------------------------------------------------------------------
# QC helpers
# ---------------------------------------------------------------------------

def calculate_qc_metrics(adata: sc.AnnData) -> sc.AnnData:
    """Calculate standard QC metrics: mito, ribo, hemoglobin."""
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    adata.var["ribo"] = adata.var_names.str.match(r"^RP[SL]\d")
    adata.var["hb"] = adata.var_names.str.match(r"^HB[^P]")

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"],
        percent_top=None, log1p=True, inplace=True
    )
    return adata


def mad_outlier(values: np.ndarray, nmads: float = 3, direction: str = "both") -> np.ndarray:
    """
    Identify outliers using Median Absolute Deviation.

    Parameters
    ----------
    values : array-like
        Metric values (e.g., log1p_total_counts).
    nmads : float
        Number of MADs for threshold.
    direction : str
        'both', 'upper', or 'lower'.

    Returns
    -------
    Boolean mask where True = outlier.
    """
    values = np.asarray(values, dtype=float)
    med = np.nanmedian(values)
    mad = np.nanmedian(np.abs(values - med))
    # Avoid zero MAD for constant data
    if mad == 0:
        mad = np.nanstd(values)
    if mad == 0:
        return np.zeros(len(values), dtype=bool)

    lower = med - nmads * 1.4826 * mad
    upper = med + nmads * 1.4826 * mad

    if direction == "both":
        return (values < lower) | (values > upper)
    elif direction == "upper":
        return values > upper
    elif direction == "lower":
        return values < lower
    else:
        raise ValueError(f"Unknown direction: {direction}")

# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def set_plotting_defaults():
    """Configure matplotlib/scanpy plotting defaults."""
    sc.settings.set_figure_params(
        dpi=100, dpi_save=300, frameon=False,
        figsize=(6, 4), facecolor="white"
    )
    plt.rcParams.update({
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
    })


def plot_qc_violins(adata: sc.AnnData, sample_col: str = "sample_id",
                    save_path: str | None = None) -> None:
    """Plot per-sample violin plots of QC metrics."""
    metrics = ["n_genes_by_counts", "total_counts", "pct_counts_mt",
               "pct_counts_ribo", "pct_counts_hb"]
    fig, axes = plt.subplots(len(metrics), 1, figsize=(max(12, adata.obs[sample_col].nunique() * 0.4), 4 * len(metrics)))
    for ax, metric in zip(axes, metrics):
        if metric not in adata.obs.columns:
            continue
        sc.pl.violin(adata, metric, groupby=sample_col, rotation=90,
                     ax=ax, show=False, stripplot=False)
        ax.set_title(metric)
    plt.tight_layout()
    if save_path:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


def plot_filtering_summary(pre_counts: dict, post_counts: dict,
                           save_path: str | None = None) -> None:
    """Bar chart of cells retained per sample after QC."""
    df = pd.DataFrame({
        "sample": list(pre_counts.keys()),
        "before": list(pre_counts.values()),
        "after": [post_counts.get(s, 0) for s in pre_counts],
    })
    df["removed"] = df["before"] - df["after"]
    df["pct_retained"] = (df["after"] / df["before"] * 100).round(1)

    fig, ax = plt.subplots(figsize=(max(10, len(df) * 0.5), 5))
    x = np.arange(len(df))
    ax.bar(x, df["before"], label="Before QC", alpha=0.5, color="steelblue")
    ax.bar(x, df["after"], label="After QC", alpha=0.8, color="darkorange")
    ax.set_xticks(x)
    ax.set_xticklabels(df["sample"], rotation=90)
    ax.set_ylabel("Number of cells")
    ax.set_title("Cells retained after QC filtering")
    ax.legend()

    # Annotate % retained
    for i, row in df.iterrows():
        ax.text(i, row["after"] + 20, f"{row['pct_retained']}%",
                ha="center", va="bottom", fontsize=7)

    plt.tight_layout()
    if save_path:
        Path(save_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()


def plot_umap_batch(adata: sc.AnnData, batch_key: str = "Batch",
                    save_path: str | None = None) -> None:
    """UMAP colored by batch for integration assessment."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    sc.pl.umap(adata, color=batch_key, ax=axes[0], show=False, title=f"UMAP by {batch_key}")
    sc.pl.umap(adata, color="leiden_0.8", ax=axes[1], show=False, title="UMAP by cluster")
    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, bbox_inches="tight")
        plt.close(fig)
    else:
        plt.show()
