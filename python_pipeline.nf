#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * =============================================================================
 * scRNA-seq + VDJ Python Pipeline — AWS Batch
 * =============================================================================
 *
 * DSL2 workflow wrapping all 9 Python pipeline phases (00–08).
 * Designed for the AWS Batch executor with Spot instances.
 *
 * Usage:
 *   nextflow run python_pipeline.nf \
 *       -c nextflow.aws.config \
 *       -work-dir s3://kstme-scrna/nextflow-work \
 *       --arm primary \
 *       -with-report report.html \
 *       -resume
 *
 * Scatter/gather strategy:
 *   Phases 00-02: parallelised per-sample (scatter)
 *   Phase 03:     gathers all QC outputs (single job)
 *   Phases 04-08: serial single jobs (07 & 08 run in parallel)
 *
 * Staging strategy:
 *   Python scripts read input paths from config.yaml (e.g., output/cellbender/).
 *   Nextflow stages inputs as flat files in the work directory. Each process
 *   creates the directory structure the Python scripts expect and symlinks
 *   staged files into place, so no Python code changes are needed.
 * =============================================================================
 */


// ---------------------------------------------------------------------------
// Phase 00 — Cell Ranger Multi
// ---------------------------------------------------------------------------
// The sample channel is built from metadata in the workflow block (Groovy).
// Cell Ranger config CSVs are generated inline since the Python script
// (00_cellranger_multi.py) needs local FASTQ directory access for path
// resolution, which isn't available until after S3 staging.
//
// The shell: directive separates Nextflow (!{var}) from bash ($var) cleanly.
// ---------------------------------------------------------------------------

process phase00_cellranger_multi {
    label 'high_mem_cpu'
    container "${params.cpu_image}"
    publishDir "${params.results}/cellranger/", mode: 'copy'

    input:
    tuple val(sample_name), val(expected_cells),
          val(all_repo_ids), val(all_loci), path("fastqs/*")

    output:
    tuple val(sample_name), path("${sample_name}/"), emit: cellranger_out

    shell:
    '''
    # Generate Cell Ranger multi config CSV
    {
    echo "[gene-expression]"
    echo "reference,/opt/references/refdata-gex-GRCh38-2024-A"
    echo "chemistry,auto"
    echo "expect-cells,!{expected_cells}"

    # Check for VDJ libraries
    IFS=',' read -ra LOCI <<< "!{all_loci}"
    HAS_VDJ=false
    for l in "${LOCI[@]}"; do
        if [[ "$l" == "5primeVDJ" ]]; then HAS_VDJ=true; break; fi
    done
    if [[ "$HAS_VDJ" == "true" ]]; then
        echo "[vdj]"
        echo "reference,/opt/references/refdata-cellranger-vdj-GRCh38"
    fi

    echo "[libraries]"
    echo "fastq_id,fastqs,feature_types,subsample_rate"

    IFS=',' read -ra REPOS <<< "!{all_repo_ids}"
    N=${#REPOS[@]}
    for ((i=0; i<N; i++)); do
        if [[ "${LOCI[i]}" == "5primeGEX" ]]; then
            echo "${REPOS[i]},fastqs/${REPOS[i]},Gene Expression,"
        elif [[ "${LOCI[i]}" == "5primeVDJ" ]]; then
            echo "${REPOS[i]},fastqs/${REPOS[i]},VDJ-T,"
        fi
    done
    } > !{sample_name}_config.csv

    cellranger multi \
        --id=!{sample_name} \
        --csv=!{sample_name}_config.csv \
        --localcores !{task.cpus} \
        --localmem $(( !{task.memory.toGiga()} - 10 ))
    '''
}


// ---------------------------------------------------------------------------
// Phase 01 — CellBender (GPU)
// ---------------------------------------------------------------------------
// The Python script reads Cell Ranger output from cfg["paths"]["cellranger_output"]
// We symlink the staged directory to that path.
// ---------------------------------------------------------------------------

process phase01_cellbender {
    label 'gpu'
    container "${params.gpu_image}"
    publishDir "${params.results}/cellbender/", mode: 'copy'

    input:
    tuple val(sample_name), path(cellranger_dir)
    path config
    path metadata

    output:
    tuple val(sample_name), path("output/cellbender/${sample_name}_cellbender.h5ad"), emit: cellbender_out
    path "output/cellbender/cellbender_summary.csv", optional: true

    script:
    """
    # Stage inputs at paths the Python script expects (per config.yaml)
    mkdir -p output/cellranger data
    ln -sf "\$(pwd)/${cellranger_dir}" "output/cellranger/${sample_name}"
    ln -sf "\$(pwd)/${metadata}" "data/kstme_sc_metadata.csv"

    python /app/pipeline/01_cellbender.py \\
        --config ${config} \\
        --sample ${sample_name}
    """
}


// ---------------------------------------------------------------------------
// Phase 02 — QC Filtering
// ---------------------------------------------------------------------------
// The Python script reads CellBender output from cfg["paths"]["cellbender_output"]
// ---------------------------------------------------------------------------

process phase02_qc {
    label 'low_mem_cpu'
    container "${params.cpu_image}"
    publishDir "${params.results}/qc/", mode: 'copy'

    input:
    tuple val(sample_name), path(cellbender_h5ad)
    path config
    path metadata

    output:
    tuple val(sample_name), path("output/qc/${sample_name}_qc.h5ad"), emit: qc_out
    path "output/figures/qc/*.png", optional: true

    script:
    """
    # Stage inputs at paths the Python script expects
    mkdir -p output/cellbender data
    ln -sf "\$(pwd)/${cellbender_h5ad}" "output/cellbender/${sample_name}_cellbender.h5ad"
    ln -sf "\$(pwd)/${metadata}" "data/kstme_sc_metadata.csv"

    python /app/pipeline/02_qc_filtering.py \\
        --config ${config} \\
        --sample ${sample_name}
    """
}


// ---------------------------------------------------------------------------
// Phase 03 — Merge + scVI Integration (GPU, gathers all samples)
// ---------------------------------------------------------------------------
// The Python script reads QC h5ad files from cfg["paths"]["qc_output"]
// ---------------------------------------------------------------------------

process phase03_integration {
    label 'gpu_large'
    container "${params.gpu_image}"
    publishDir "${params.results}/integration/", mode: 'copy'

    input:
    path qc_files       // all per-sample QC h5ad files (collected)
    path config
    path metadata
    path clinical_metadata
    val arm

    output:
    path "output/integration/integrated_${arm}.h5ad", emit: integrated
    path "output/integration/integration_metrics_${arm}.csv", optional: true
    path "output/figures/integration/*.png", optional: true

    script:
    """
    # Stage QC h5ad files at the paths the Python script expects
    mkdir -p output/qc data
    for f in ${qc_files}; do
        ln -sf "\$(pwd)/\$f" "output/qc/\$f"
    done

    # Stage metadata at config-expected paths
    ln -sf "\$(pwd)/${metadata}" "data/kstme_sc_metadata.csv"
    ln -sf "\$(pwd)/${clinical_metadata}" "data/KSTME_updated_clinical_metadata.csv"

    python /app/pipeline/03_merge_integrate.py \\
        --config ${config} \\
        --arm ${arm}
    """
}


// ---------------------------------------------------------------------------
// Phase 04 — Cell Type Annotation
// ---------------------------------------------------------------------------
// The Python script reads from cfg["paths"]["integration_output"]
// ---------------------------------------------------------------------------

process phase04_annotation {
    label 'mid_mem_cpu'
    container "${params.cpu_image}"
    publishDir "${params.results}/annotation/", mode: 'copy'

    input:
    path integrated_h5ad
    path config
    val arm

    output:
    path "output/annotation/annotated_${arm}.h5ad", emit: annotated
    path "output/annotation/cluster_markers_${arm}.csv", optional: true
    path "output/figures/annotation/*.png", optional: true

    script:
    """
    # Stage input at the path the Python script expects
    mkdir -p output/integration
    ln -sf "\$(pwd)/${integrated_h5ad}" "output/integration/integrated_${arm}.h5ad"

    python /app/pipeline/04_cluster_annotate.py \\
        --config ${config} \\
        --arm ${arm}
    """
}


// ---------------------------------------------------------------------------
// Phase 05 — T Cell Subclustering (GPU)
// ---------------------------------------------------------------------------
// The Python script reads from cfg["paths"]["annotation_output"]
// ---------------------------------------------------------------------------

process phase05_tcell {
    label 'gpu'
    container "${params.gpu_image}"
    publishDir "${params.results}/tcell/", mode: 'copy'

    input:
    path annotated_h5ad
    path config
    val arm

    output:
    path "output/tcell/tcell_${arm}.h5ad", emit: tcell
    path "output/figures/tcell/*.png", optional: true

    script:
    """
    # Stage input at the path the Python script expects
    mkdir -p output/annotation
    ln -sf "\$(pwd)/${annotated_h5ad}" "output/annotation/annotated_${arm}.h5ad"

    python /app/pipeline/05_tcell_subcluster.py \\
        --config ${config} \\
        --arm ${arm}
    """
}


// ---------------------------------------------------------------------------
// Phase 06 — VDJ Integration
// ---------------------------------------------------------------------------
// The Python script reads T cell data from cfg["paths"]["tcell_output"]
// and VDJ contigs from cfg["paths"]["cellranger_output"]
// ---------------------------------------------------------------------------

process phase06_vdj {
    label 'low_mem_cpu'
    container "${params.cpu_image}"
    publishDir "${params.results}/vdj/", mode: 'copy'

    input:
    path tcell_h5ad
    path cellranger_dirs   // all Cell Ranger output dirs (collected)
    path config
    path metadata
    val arm

    output:
    path "output/vdj/tcell_vdj_${arm}.h5ad", emit: vdj
    path "output/vdj/vgene_usage_${arm}.csv", optional: true
    path "output/figures/vdj/*.png", optional: true

    script:
    """
    # Stage T cell h5ad at the path the Python script expects
    mkdir -p output/tcell output/cellranger data
    ln -sf "\$(pwd)/${tcell_h5ad}" "output/tcell/tcell_${arm}.h5ad"

    # Stage Cell Ranger output dirs (for VDJ contigs)
    for d in ${cellranger_dirs}; do
        ln -sf "\$(pwd)/\$d" "output/cellranger/\$d"
    done

    # Stage metadata
    ln -sf "\$(pwd)/${metadata}" "data/kstme_sc_metadata.csv"

    python /app/pipeline/06_vdj_integration.py \\
        --config ${config} \\
        --arm ${arm}
    """
}


// ---------------------------------------------------------------------------
// Phase 07 — Differential Expression
// ---------------------------------------------------------------------------
// The Python script reads from cfg["paths"]["annotation_output"]
// ---------------------------------------------------------------------------

process phase07_de {
    label 'low_mem_cpu'
    container "${params.cpu_image}"
    publishDir "${params.results}/differential_expression/", mode: 'copy'

    input:
    path annotated_h5ad
    path config
    val arm

    output:
    path "output/differential_expression/de_results_*_${arm}.csv", optional: true, emit: de_results
    path "output/differential_expression/cell_proportions_*_${arm}.csv", optional: true
    path "output/figures/de/*.png", optional: true

    script:
    """
    # Stage input at the path the Python script expects
    mkdir -p output/annotation
    ln -sf "\$(pwd)/${annotated_h5ad}" "output/annotation/annotated_${arm}.h5ad"

    python /app/pipeline/07_differential_expression.py \\
        --config ${config} \\
        --arm ${arm}
    """
}


// ---------------------------------------------------------------------------
// Phase 08 — Trajectory Analysis
// ---------------------------------------------------------------------------
// The Python script reads from cfg["paths"]["vdj_output"]
// ---------------------------------------------------------------------------

process phase08_trajectory {
    label 'low_mem_cpu'
    container "${params.cpu_image}"
    publishDir "${params.results}/trajectory/", mode: 'copy'

    input:
    path vdj_h5ad
    path config
    val arm

    output:
    path "output/trajectory/trajectory_${arm}.h5ad", optional: true, emit: trajectory
    path "output/trajectory/clonotype_trajectories_${arm}.csv", optional: true
    path "output/figures/trajectory/*.png", optional: true

    script:
    """
    # Stage input at the path the Python script expects
    mkdir -p output/vdj
    ln -sf "\$(pwd)/${vdj_h5ad}" "output/vdj/tcell_vdj_${arm}.h5ad"

    python /app/pipeline/08_trajectory.py \\
        --config ${config} \\
        --arm ${arm}
    """
}


// ===========================================================================
// Workflow
// ===========================================================================

workflow {

    // -----------------------------------------------------------------------
    // Input channels (files on S3, staged by Nextflow)
    // -----------------------------------------------------------------------
    config_ch          = Channel.fromPath(params.config)
    metadata_ch        = Channel.fromPath(params.metadata)
    clinical_meta_ch   = Channel.fromPath(params.clinical_metadata)

    // -----------------------------------------------------------------------
    // Build sample channel from metadata CSV
    // -----------------------------------------------------------------------
    // Parse metadata, group by sampleName. Each sample gets a tuple with
    // all its library info (GEX + VDJ) and corresponding FASTQ directories.
    // The AAAMMK2M5 batch is excluded (FASTQs missing; data in HVFNLDRXX).
    //
    // Channel emits: (sampleName, expected_cells, repo_ids_csv, loci_csv, [fastq_dirs])

    sample_ch = Channel.fromPath(params.metadata)
        .splitCsv(header: true)
        .filter { it.Batch != 'AAAMMK2M5' }
        .map { row ->
            tuple(
                row.sampleName,
                row.expected_cells,
                row.repoID,
                row.locus,
                file("${params.fastq_base}/${row.Batch}/${row.repoID}")
            )
        }
        .groupTuple(by: 0)
        // After groupTuple: (sampleName, [expected_cells...], [repoIDs...], [loci...], [dirs...])
        .map { sampleName, expected_list, repo_ids, loci, fastq_dirs ->
            tuple(
                sampleName,
                expected_list[0],
                repo_ids.join(','),
                loci.join(','),
                fastq_dirs
            )
        }

    // -----------------------------------------------------------------------
    // Phase 00: Cell Ranger multi (per-sample, scatter)
    // -----------------------------------------------------------------------
    phase00_cellranger_multi(sample_ch)

    // -----------------------------------------------------------------------
    // Phase 01: CellBender (per-sample, GPU)
    // -----------------------------------------------------------------------
    phase01_cellbender(
        phase00_cellranger_multi.out.cellranger_out,
        config_ch,
        metadata_ch
    )

    // -----------------------------------------------------------------------
    // Phase 02: QC (per-sample, CPU)
    // -----------------------------------------------------------------------
    phase02_qc(
        phase01_cellbender.out.cellbender_out,
        config_ch,
        metadata_ch
    )

    // -----------------------------------------------------------------------
    // Phase 03: Integration (gather all QC outputs, single GPU job)
    // -----------------------------------------------------------------------
    all_qc_h5ads = phase02_qc.out.qc_out
        .map { sample_name, h5ad -> h5ad }
        .collect()

    phase03_integration(
        all_qc_h5ads,
        config_ch,
        metadata_ch,
        clinical_meta_ch,
        params.arm
    )

    // -----------------------------------------------------------------------
    // Phase 04: Annotation (single CPU job)
    // -----------------------------------------------------------------------
    phase04_annotation(
        phase03_integration.out.integrated,
        config_ch,
        params.arm
    )

    // -----------------------------------------------------------------------
    // Phase 05: T Cell Subclustering (single GPU job)
    // -----------------------------------------------------------------------
    phase05_tcell(
        phase04_annotation.out.annotated,
        config_ch,
        params.arm
    )

    // -----------------------------------------------------------------------
    // Phase 06: VDJ Integration (needs T cell data + Cell Ranger VDJ contigs)
    // -----------------------------------------------------------------------
    all_cr_dirs = phase00_cellranger_multi.out.cellranger_out
        .map { sample_name, cr_dir -> cr_dir }
        .collect()

    phase06_vdj(
        phase05_tcell.out.tcell,
        all_cr_dirs,
        config_ch,
        metadata_ch,
        params.arm
    )

    // -----------------------------------------------------------------------
    // Phases 07 & 08 run in parallel (different input dependencies)
    // -----------------------------------------------------------------------

    // Phase 07: DE (depends on Phase 04 annotated data)
    phase07_de(
        phase04_annotation.out.annotated,
        config_ch,
        params.arm
    )

    // Phase 08: Trajectory (depends on Phase 06 VDJ data)
    phase08_trajectory(
        phase06_vdj.out.vdj,
        config_ch,
        params.arm
    )
}
