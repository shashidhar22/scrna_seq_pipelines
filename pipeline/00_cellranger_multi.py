#!/usr/bin/env python3
"""
Phase 0: Generate Cell Ranger multi config CSVs from metadata.

Parses KSTME_metadata.csv, pairs GEX + VDJ libraries per sample,
resolves FASTQ paths across sequencing batches, and handles the
AAAMMK2M5 re-sequencing logic.

Usage:
    python pipeline/00_cellranger_multi.py [--config pipeline/config.yaml]
    python pipeline/00_cellranger_multi.py --verify  # Check Cell Ranger outputs
"""

import argparse
import csv
import sys
from pathlib import Path

import pandas as pd

from utils import load_config, load_metadata, get_gex_samples, get_vdj_samples, setup_logging

logger = setup_logging("00_cellranger_multi")


def resolve_fastq_paths(repo_id: str, batch: str, fastq_base: str) -> list[str]:
    """
    Find FASTQ directory for a given repoID and batch.

    The mkfastq output structure is:
      {fastq_base}/{batch}/outs/fastq_path/{batch}/{repoID}/
    or sometimes directly:
      {fastq_base}/{batch}/{repoID}/

    Returns list of valid FASTQ directories found.
    """
    candidates = [
        Path(fastq_base) / batch / "outs" / "fastq_path" / batch / repo_id,
        Path(fastq_base) / batch / repo_id,
        # Some mkfastq outputs use MAKE_FASTQS_CS path
        Path(fastq_base) / batch / "MAKE_FASTQS_CS" / "MAKE_FASTQS" / "BCL2FASTQ_WITH_SAMPLESHEET" / "fork0" / "files" / "fastq_path" / batch / repo_id,
    ]
    return [str(p) for p in candidates if p.is_dir()]


def build_sample_table(meta: pd.DataFrame, cfg: dict) -> pd.DataFrame:
    """
    Build a table of GEX-VDJ paired samples with resolved FASTQ paths.

    For re-sequenced samples (appearing in both AAAMMK2M5 and HVFNLDRXX),
    merge FASTQ paths from all available flowcells.
    """
    fastq_base = cfg["paths"]["fastq_base"]
    missing_batch = cfg["study"]["missing_batch"]
    reseq_samples = {s["sample"]: s["batches"]
                     for s in cfg["study"].get("resequenced_samples", [])}

    rows = []
    # Group by sampleName to pair GEX + VDJ
    for sample_name, group in meta.groupby("sampleName"):
        gex_rows = group[group["locus"] == "5primeGEX"]
        vdj_rows = group[group["locus"] == "5primeVDJ"]

        sample_info = group.iloc[0]
        patient_id = sample_info["patientID"]
        expected_cells = sample_info["expected_cells"]

        # Collect GEX FASTQ paths
        gex_fastqs = []
        for _, row in gex_rows.iterrows():
            batch = row["Batch"]
            repo_id = row["repoID"]
            if batch == missing_batch:
                # Skip missing batch - FASTQs should be in the paired batch
                logger.warning(
                    f"Sample {sample_name}: skipping missing batch {missing_batch} "
                    f"(repo: {repo_id}). Data expected in alternate batch."
                )
                continue
            paths = resolve_fastq_paths(repo_id, batch, fastq_base)
            if paths:
                gex_fastqs.extend(paths)
            else:
                logger.warning(f"No FASTQ dir found for GEX {repo_id} in batch {batch}")

        # Collect VDJ FASTQ paths
        vdj_fastqs = []
        for _, row in vdj_rows.iterrows():
            batch = row["Batch"]
            repo_id = row["repoID"]
            if batch == missing_batch:
                logger.warning(
                    f"Sample {sample_name}: skipping missing batch {missing_batch} "
                    f"(repo: {repo_id}). Data expected in alternate batch."
                )
                continue
            paths = resolve_fastq_paths(repo_id, batch, fastq_base)
            if paths:
                vdj_fastqs.extend(paths)
            else:
                logger.warning(f"No FASTQ dir found for VDJ {repo_id} in batch {batch}")

        rows.append({
            "sample_name": sample_name,
            "patient_id": patient_id,
            "expected_cells": expected_cells,
            "gex_fastqs": gex_fastqs,
            "gex_fastq_ids": [r["repoID"] for _, r in gex_rows.iterrows()
                              if r["Batch"] != missing_batch],
            "vdj_fastqs": vdj_fastqs,
            "vdj_fastq_ids": [r["repoID"] for _, r in vdj_rows.iterrows()
                              if r["Batch"] != missing_batch],
            "has_gex": len(gex_fastqs) > 0,
            "has_vdj": len(vdj_fastqs) > 0,
        })

    return pd.DataFrame(rows)


def write_multi_config(sample_row: pd.Series, cfg: dict, output_dir: str) -> str:
    """
    Write a Cell Ranger multi config CSV for a single sample.

    Config format:
        [gene-expression]
        reference,/path/to/gex/ref
        chemistry,auto
        expect-cells,3000
        [vdj]
        reference,/path/to/vdj/ref
        [libraries]
        fastq_id,fastqs,feature_types,subsample_rate
        ...
    """
    sample_name = sample_row["sample_name"]
    out_path = Path(output_dir) / f"{sample_name}_config.csv"

    gex_ref = cfg["paths"]["gex_reference"]
    vdj_ref = cfg["paths"]["vdj_reference"]
    expected = sample_row["expected_cells"]

    lines = []

    # GEX section
    if sample_row["has_gex"]:
        lines.extend([
            "[gene-expression]",
            f"reference,{gex_ref}",
            "chemistry,auto",
            f"expect-cells,{expected}",
        ])

    # VDJ section
    if sample_row["has_vdj"]:
        lines.extend([
            "[vdj]",
            f"reference,{vdj_ref}",
        ])

    # Libraries section
    lines.append("[libraries]")
    lines.append("fastq_id,fastqs,feature_types,subsample_rate")

    # GEX libraries
    for fq_id, fq_path in zip(sample_row["gex_fastq_ids"], sample_row["gex_fastqs"]):
        lines.append(f"{fq_id},{fq_path},Gene Expression,")

    # VDJ libraries
    for fq_id, fq_path in zip(sample_row["vdj_fastq_ids"], sample_row["vdj_fastqs"]):
        lines.append(f"{fq_id},{fq_path},VDJ-T,")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n")
    logger.info(f"Wrote config: {out_path}")
    return str(out_path)


def verify_cellranger_outputs(cfg: dict, meta: pd.DataFrame) -> pd.DataFrame:
    """
    Check Cell Ranger multi outputs and flag problematic samples.

    Reads metrics_summary.csv for each sample and flags:
    - <1,000 cells recovered
    - <500 median genes per cell
    - <50% sequencing saturation
    """
    gex_samples = get_gex_samples(meta)
    # Deduplicate: one row per unique sampleName
    sample_names = gex_samples["sampleName"].unique()

    results = []
    for sample_name in sample_names:
        cr_dir = Path(cfg["paths"]["cellranger_output"]) / sample_name / "outs"
        metrics_path = cr_dir / "multi" / "count" / "summary.csv"

        if not metrics_path.exists():
            # Try alternative path
            metrics_path = cr_dir / "metrics_summary.csv"

        if not metrics_path.exists():
            results.append({
                "sample": sample_name,
                "status": "MISSING",
                "cells": None,
                "median_genes": None,
                "saturation": None,
            })
            continue

        try:
            metrics = pd.read_csv(metrics_path)
            # Cell Ranger metrics columns vary by version; try common names
            cells = None
            median_genes = None
            saturation = None

            for col in metrics.columns:
                col_lower = col.lower().strip()
                val = metrics[col].iloc[0]
                # Parse percentage/comma-formatted strings
                if isinstance(val, str):
                    val = val.replace(",", "").replace("%", "")
                    try:
                        val = float(val)
                    except ValueError:
                        continue

                if "estimated" in col_lower and "cell" in col_lower:
                    cells = val
                elif "median genes" in col_lower:
                    median_genes = val
                elif "sequencing saturation" in col_lower:
                    saturation = val

            flags = []
            if cells is not None and cells < 1000:
                flags.append("LOW_CELLS")
            if median_genes is not None and median_genes < 500:
                flags.append("LOW_GENES")
            if saturation is not None and saturation < 50:
                flags.append("LOW_SATURATION")

            results.append({
                "sample": sample_name,
                "status": "FLAG" if flags else "OK",
                "flags": ",".join(flags) if flags else "",
                "cells": cells,
                "median_genes": median_genes,
                "saturation": saturation,
            })
        except Exception as e:
            results.append({
                "sample": sample_name,
                "status": "ERROR",
                "cells": None,
                "median_genes": None,
                "saturation": None,
                "error": str(e),
            })

    df = pd.DataFrame(results)
    return df


def main():
    parser = argparse.ArgumentParser(description="Generate Cell Ranger multi configs")
    parser.add_argument("--config", default="pipeline/config.yaml")
    parser.add_argument("--verify", action="store_true",
                        help="Verify Cell Ranger outputs instead of generating configs")
    parser.add_argument("--output-dir", default=None,
                        help="Output directory for config files (default: cellranger_output)")
    args = parser.parse_args()

    cfg = load_config(args.config)
    meta = load_metadata(cfg)

    if args.verify:
        logger.info("Verifying Cell Ranger multi outputs...")
        results = verify_cellranger_outputs(cfg, meta)
        out_path = Path(cfg["paths"]["cellranger_output"]) / "verification_summary.csv"
        out_path.parent.mkdir(parents=True, exist_ok=True)
        results.to_csv(out_path, index=False)
        logger.info(f"Verification summary: {out_path}")

        # Report
        n_ok = (results["status"] == "OK").sum()
        n_flag = (results["status"] == "FLAG").sum()
        n_missing = (results["status"] == "MISSING").sum()
        logger.info(f"Results: {n_ok} OK, {n_flag} flagged, {n_missing} missing")
        if n_flag > 0:
            logger.warning("Flagged samples:")
            for _, row in results[results["status"] == "FLAG"].iterrows():
                logger.warning(f"  {row['sample']}: {row.get('flags', '')}")
        return

    # Generate configs
    output_dir = args.output_dir or Path(cfg["paths"]["cellranger_output"]) / "configs"
    logger.info("Building sample table from metadata...")
    sample_table = build_sample_table(meta, cfg)

    logger.info(f"Found {len(sample_table)} unique samples")
    logger.info(f"  With GEX: {sample_table['has_gex'].sum()}")
    logger.info(f"  With VDJ: {sample_table['has_vdj'].sum()}")

    # Write individual config CSVs
    config_paths = []
    for _, row in sample_table.iterrows():
        if not row["has_gex"]:
            logger.warning(f"Skipping {row['sample_name']}: no GEX FASTQs found")
            continue
        path = write_multi_config(row, cfg, output_dir)
        config_paths.append({"sample_name": row["sample_name"], "config_file": path})

    # Write master config table
    config_df = pd.DataFrame(config_paths)
    master_path = Path(output_dir) / "multi_config.tsv"
    config_df.to_csv(master_path, sep="\t", index=False)
    logger.info(f"Wrote master config table: {master_path} ({len(config_df)} samples)")

    # Report samples missing FASTQ data
    missing_gex = sample_table[~sample_table["has_gex"]]
    if len(missing_gex) > 0:
        logger.warning(f"{len(missing_gex)} samples missing GEX FASTQs:")
        for _, row in missing_gex.iterrows():
            logger.warning(f"  {row['sample_name']} (patient: {row['patient_id']})")


if __name__ == "__main__":
    main()
