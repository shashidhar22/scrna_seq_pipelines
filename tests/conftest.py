"""
Shared test fixtures for the scRNA-seq + VDJ pipeline test suite.
"""

import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import scanpy as sc
import anndata as ad
from scipy.sparse import csr_matrix


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

REPO_ROOT = Path(__file__).resolve().parents[1]
PIPELINE_DIR = REPO_ROOT / "pipeline"


@pytest.fixture()
def tmp_dir(tmp_path):
    """Return a temporary directory for test outputs."""
    return tmp_path


# ---------------------------------------------------------------------------
# Config fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def minimal_config(tmp_path):
    """Minimal pipeline config dict (no file I/O needed)."""
    return {
        "paths": {
            "metadata": str(tmp_path / "meta.csv"),
            "clinical_metadata": str(tmp_path / "clinical.csv"),
            "fastq_base": str(tmp_path / "fastq"),
            "gex_reference": "/ref/gex",
            "vdj_reference": "/ref/vdj",
            "cellranger_output": str(tmp_path / "cellranger"),
            "cellbender_output": str(tmp_path / "cellbender"),
            "qc_output": str(tmp_path / "qc"),
            "integration_output": str(tmp_path / "integration"),
            "annotation_output": str(tmp_path / "annotation"),
            "tcell_output": str(tmp_path / "tcell"),
            "vdj_output": str(tmp_path / "vdj"),
            "de_output": str(tmp_path / "de"),
            "trajectory_output": str(tmp_path / "trajectory"),
            "figures": str(tmp_path / "figures"),
            "models": str(tmp_path / "models"),
        },
        "study": {
            "study_id": "TEST",
            "expected_cells": 3000,
            "chemistry": "auto",
            "missing_batch": "AAAMMK2M5",
            "resequenced_samples": [
                {"sample": "008_216_V09", "batches": ["HVFNLDRXX", "AAAMMK2M5"]},
            ],
            "hiv_status_overrides": {"008_252": "negative"},
            "analysis_arms": {
                "primary": {
                    "description": "Unsorted",
                    "sorting_filter": ["Unsorted"],
                },
                "focused": {
                    "description": "Sorted",
                    "sorting_filter": ["CD4pTCells", "CD8pTCells"],
                },
                "validation": {
                    "description": "Both",
                    "patients": ["008_216", "008_217", "008_220"],
                },
            },
        },
        "cellbender": {
            "fpr": 0.01,
            "total_droplets_included": 20000,
            "epochs": 150,
            "learning_rate": 1e-4,
            "max_epochs_fallback": 300,
            "convergence_check": True,
        },
        "qc": {
            "mad_threshold": 3,
            "max_pct_hb": 5.0,
            "min_cells_per_gene": 3,
            "expected_doublet_rate": 0.046,
            "min_cells_per_library": 200,
        },
        "integration": {
            "n_top_genes": 4000,
            "hvg_flavor": "seurat_v3",
            "scvi": {
                "batch_key": "Batch",
                "n_latent": 30,
                "n_layers": 2,
                "max_epochs": 400,
                "early_stopping": True,
                "continuous_covariates": ["pct_counts_mt"],
            },
            "harmony": {"n_pcs": 50, "batch_key": "Batch"},
            "n_neighbors": 30,
            "leiden_resolutions": [0.3, 0.5, 0.8, 1.0, 1.5],
            "default_resolution": 0.8,
            "umap_min_dist": 0.3,
        },
        "annotation": {
            "celltypist_models": ["Immune_All_Low.pkl", "Immune_All_High.pkl"],
            "majority_voting": True,
            "markers": {
                "T_cells": ["CD3D", "CD3E"],
                "CD8_T": ["CD8A", "CD8B"],
            },
        },
        "tcell": {
            "n_top_genes": 3000,
            "scvi": {"n_latent": 20, "max_epochs": 300, "batch_key": "Batch"},
            "n_neighbors": 20,
            "resolution": 1.0,
            "cd8_markers": {
                "Naive": ["CCR7", "SELL", "LEF1"],
                "Effector": ["GZMB", "PRF1"],
                "Exhausted": ["PDCD1", "LAG3", "HAVCR2", "TOX"],
            },
            "cd4_markers": {
                "Naive": ["CCR7"],
                "Treg": ["FOXP3", "IL2RA"],
            },
        },
        "vdj": {
            "metric": "identity",
            "sequence": "aa",
            "receptor_arms": "all",
            "dual_ir": "primary_only",
            "expansion_clip_at": 4,
            "diversity_metrics": ["shannon", "simpson"],
        },
        "de": {
            "min_cells_per_pseudobulk": 10,
            "min_counts_per_pseudobulk": 1000,
            "sample_col": "patientID",
            "contrasts": [
                {
                    "name": "HIV_status",
                    "column": "HIVstatus",
                    "reference": "negative",
                    "target": "positive",
                },
            ],
        },
        "trajectory": {
            "n_dcs": 15,
            "root_cell_type": "Naive",
        },
        "compute": {
            "n_jobs": 4,
            "gpu": True,
            "random_seed": 42,
        },
    }


# ---------------------------------------------------------------------------
# Metadata fixtures
# ---------------------------------------------------------------------------

@pytest.fixture()
def sample_metadata_df():
    """A small metadata DataFrame mirroring KSTME_metadata.csv structure."""
    rows = [
        # Unsorted GEX + VDJ pairs
        {"repoName": "R1", "sampleName": "008_216_V01", "patientID": "008_216",
         "visitCode": "V01", "expected_cells": 3000, "chemistry": "fiveprime",
         "nucliecAcid": "cDNA", "locus": "5primeGEX", "platform": "10XGenomics",
         "indices": "A5", "sortingCT": "Unsorted", "tissueType": "PBMC",
         "VDJType": "NA", "Batch": "HTVWKDRXX", "HIVstatus": "negative",
         "repoID": "R1", "Note": ""},
        {"repoName": "R2", "sampleName": "008_216_V01", "patientID": "008_216",
         "visitCode": "V01", "expected_cells": 3000, "chemistry": "fiveprime",
         "nucliecAcid": "cDNA", "locus": "5primeVDJ", "platform": "10XGenomics",
         "indices": "E5", "sortingCT": "Unsorted", "tissueType": "PBMC",
         "VDJType": "T-cell", "Batch": "HTVWKDRXX", "HIVstatus": "negative",
         "repoID": "R2", "Note": ""},
        # Another patient (HIV+)
        {"repoName": "R3", "sampleName": "008_220_V01", "patientID": "008_220",
         "visitCode": "V01", "expected_cells": 3000, "chemistry": "fiveprime",
         "nucliecAcid": "cDNA", "locus": "5primeGEX", "platform": "10XGenomics",
         "indices": "D7", "sortingCT": "Unsorted", "tissueType": "PBMC",
         "VDJType": "NA", "Batch": "HVFNLDRXX", "HIVstatus": "positive",
         "repoID": "R3", "Note": ""},
        # Sorted CD4 sample (focused arm)
        {"repoName": "R4", "sampleName": "008_216_V02", "patientID": "008_216",
         "visitCode": "V02", "expected_cells": 3000, "chemistry": "fiveprime",
         "nucliecAcid": "cDNA", "locus": "5primeGEX", "platform": "10XGenomics",
         "indices": "A2", "sortingCT": "CD4pTCells", "tissueType": "CD4 sorted PBMC",
         "VDJType": "NA", "Batch": "HGMLLDRXX", "HIVstatus": "negative",
         "repoID": "R4", "Note": ""},
        # Sample in missing batch
        {"repoName": "R5", "sampleName": "008_216_V09", "patientID": "008_216",
         "visitCode": "V09", "expected_cells": 3000, "chemistry": "fiveprime",
         "nucliecAcid": "cDNA", "locus": "5primeGEX", "platform": "10XGenomics",
         "indices": "A7", "sortingCT": "Unsorted", "tissueType": "PBMC",
         "VDJType": "NA", "Batch": "AAAMMK2M5", "HIVstatus": "negative",
         "repoID": "R5", "Note": ""},
        # Patient with HIV status override (mislabeled VDJ row)
        {"repoName": "R6", "sampleName": "008_252_V09", "patientID": "008_252",
         "visitCode": "V09", "expected_cells": 3000, "chemistry": "fiveprime",
         "nucliecAcid": "cDNA", "locus": "5primeVDJ", "platform": "10XGenomics",
         "indices": "D6", "sortingCT": "Unsorted", "tissueType": "PBMC",
         "VDJType": "T-cell", "Batch": "HMVGWDRX2", "HIVstatus": "positive",
         "repoID": "R6", "Note": "Sample was mislabeled 225."},
        {"repoName": "R7", "sampleName": "008_252_V01", "patientID": "008_252",
         "visitCode": "V01", "expected_cells": 3000, "chemistry": "fiveprime",
         "nucliecAcid": "cDNA", "locus": "5primeGEX", "platform": "10XGenomics",
         "indices": "E8", "sortingCT": "Unsorted", "tissueType": "PBMC",
         "VDJType": "NA", "Batch": "HJTLFDRX2", "HIVstatus": "negative",
         "repoID": "R7", "Note": ""},
    ]
    return pd.DataFrame(rows)


@pytest.fixture()
def sample_metadata_csv(tmp_path, sample_metadata_df):
    """Write sample metadata to a CSV and return path."""
    path = tmp_path / "meta.csv"
    sample_metadata_df.to_csv(path, index=False)
    return str(path)


@pytest.fixture()
def clinical_metadata_df():
    """Small clinical metadata DataFrame."""
    return pd.DataFrame({
        "ptid": [216, 220, 252, 999],
        "best_resp": ["CR", "PR", "SD", "PD"],
        "bl_plasma_grp": ["low", "high", "low", "high"],
        "bl_plasma_grp_adj": ["low", "high", "low", "high"],
        "best_resp_new": ["CR", "PR", "SD", "PD"],
    })


@pytest.fixture()
def clinical_metadata_csv(tmp_path, clinical_metadata_df):
    """Write clinical metadata to a CSV and return path."""
    path = tmp_path / "clinical.csv"
    clinical_metadata_df.to_csv(path, index=False)
    return str(path)


@pytest.fixture()
def config_with_metadata(minimal_config, sample_metadata_csv, clinical_metadata_csv):
    """Config pointing to real CSV files."""
    cfg = minimal_config.copy()
    cfg["paths"]["metadata"] = sample_metadata_csv
    cfg["paths"]["clinical_metadata"] = clinical_metadata_csv
    return cfg


# ---------------------------------------------------------------------------
# AnnData fixtures
# ---------------------------------------------------------------------------

def _make_adata(n_obs=200, n_vars=100, n_genes_prefix=None, seed=42):
    """Create a synthetic AnnData for testing."""
    rng = np.random.default_rng(seed)
    counts = rng.poisson(lam=3, size=(n_obs, n_vars)).astype(np.float32)
    X = csr_matrix(counts)

    gene_names = [f"Gene{i}" for i in range(n_vars)]
    if n_genes_prefix:
        for prefix, idxs in n_genes_prefix.items():
            for idx in idxs:
                if idx < n_vars:
                    gene_names[idx] = f"{prefix}{idx}"

    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)])
    var = pd.DataFrame(index=gene_names)
    return ad.AnnData(X=X, obs=obs, var=var)


@pytest.fixture()
def simple_adata():
    """Simple 200-cell, 100-gene AnnData."""
    return _make_adata(200, 100)


@pytest.fixture()
def adata_with_qc_genes():
    """AnnData with MT-, RP, HB gene prefixes for QC testing."""
    n_obs, n_vars = 500, 120
    rng = np.random.default_rng(42)
    counts = rng.poisson(lam=5, size=(n_obs, n_vars)).astype(np.float32)

    gene_names = [f"Gene{i}" for i in range(n_vars)]
    # Add MT genes (indices 0-4)
    for i in range(5):
        gene_names[i] = f"MT-GENE{i}"
    # Add ribo genes (indices 5-9)
    for i in range(5, 10):
        gene_names[i] = f"RPS{i}"
    # Add HB genes (indices 10-12)
    gene_names[10] = "HBA1"
    gene_names[11] = "HBA2"
    gene_names[12] = "HBB"
    # Add canonical markers for annotation tests
    gene_names[13] = "CD3D"
    gene_names[14] = "CD3E"
    gene_names[15] = "CD4"
    gene_names[16] = "CD8A"
    gene_names[17] = "CD8B"
    gene_names[18] = "MS4A1"
    gene_names[19] = "CD14"
    gene_names[20] = "NKG7"
    gene_names[21] = "PPBP"
    gene_names[22] = "IRF7"
    gene_names[23] = "GNLY"
    gene_names[24] = "FOXP3"
    gene_names[25] = "IL2RA"
    gene_names[26] = "CCR7"
    gene_names[27] = "SELL"
    gene_names[28] = "LEF1"
    gene_names[29] = "GZMB"
    gene_names[30] = "PRF1"
    gene_names[31] = "PDCD1"
    gene_names[32] = "LAG3"
    gene_names[33] = "HAVCR2"
    gene_names[34] = "TOX"
    gene_names[35] = "MKI67"
    gene_names[36] = "TOP2A"
    gene_names[37] = "IL7R"

    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(n_obs)])
    var = pd.DataFrame(index=gene_names)
    adata = ad.AnnData(X=csr_matrix(counts), obs=obs, var=var)
    return adata


@pytest.fixture()
def adata_with_metadata(adata_with_qc_genes):
    """AnnData with sample/batch/patient metadata in .obs."""
    adata = adata_with_qc_genes.copy()
    n = adata.n_obs
    rng = np.random.default_rng(42)

    # Assign cells to samples
    samples = ["008_216_V01", "008_220_V01", "008_252_V01"]
    adata.obs["sample_id"] = rng.choice(samples, size=n)
    adata.obs["patientID"] = adata.obs["sample_id"].str.extract(r"(008_\d+)")[0].values
    adata.obs["Batch"] = adata.obs["sample_id"].map({
        "008_216_V01": "HTVWKDRXX",
        "008_220_V01": "HVFNLDRXX",
        "008_252_V01": "HJTLFDRX2",
    })
    adata.obs["HIVstatus"] = adata.obs["patientID"].map({
        "008_216": "negative",
        "008_220": "positive",
        "008_252": "negative",
    })
    adata.obs["sortingCT"] = "Unsorted"
    adata.obs["tissueType"] = "PBMC"

    # Store counts layer
    adata.layers["counts"] = adata.X.copy()

    return adata


@pytest.fixture()
def adata_annotated(adata_with_metadata):
    """AnnData with cell type annotations for downstream tests."""
    adata = adata_with_metadata.copy()
    n = adata.n_obs
    rng = np.random.default_rng(42)

    cell_types = ["CD4 T", "CD8 T", "T cells", "NK", "B cells", "Monocytes"]
    adata.obs["cell_type"] = rng.choice(cell_types, size=n, p=[0.25, 0.25, 0.1, 0.1, 0.15, 0.15])
    adata.obs["cell_subtype"] = adata.obs["cell_type"]

    # Add cluster labels
    adata.obs["leiden_0.8"] = pd.Categorical(rng.choice([str(i) for i in range(8)], size=n))
    adata.obs["leiden_1.0"] = pd.Categorical(rng.choice([str(i) for i in range(12)], size=n))

    return adata


@pytest.fixture()
def adata_tcell(adata_annotated):
    """AnnData subset to T cells with tcell_subtype annotations."""
    adata = adata_annotated.copy()
    t_mask = adata.obs["cell_type"].isin(["CD4 T", "CD8 T", "T cells"])
    adata_t = adata[t_mask].copy()
    n = adata_t.n_obs
    rng = np.random.default_rng(42)

    subtypes = ["CD8 Naive", "CD8 Effector", "CD8 Exhausted", "CD4 Naive", "CD4 Treg", "CD4 Memory"]
    adata_t.obs["tcell_subtype"] = rng.choice(subtypes, size=n)

    # Add pseudotime-like data
    adata_t.obs["dpt_pseudotime"] = rng.uniform(0, 1, size=n)

    # Add VDJ-like columns
    clone_ids = [f"clone_{i}" for i in range(20)] + [None] * 5
    adata_t.obs["clone_id"] = rng.choice(clone_ids, size=n)
    adata_t.obs["clonal_expansion"] = rng.choice(["1", "2", "3", "4+"], size=n, p=[0.6, 0.2, 0.1, 0.1])
    adata_t.obs["IR_VDJ_1_junction_aa"] = rng.choice(["CASSLAPGATNEKLFF", None], size=n, p=[0.7, 0.3])
    adata_t.obs["IR_VDJ_2_junction_aa"] = rng.choice(["CASSLGG", None], size=n, p=[0.05, 0.95])
    adata_t.obs["doublet_score"] = rng.uniform(0, 0.5, size=n)

    return adata_t
