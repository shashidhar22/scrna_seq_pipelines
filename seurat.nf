#!/usr/bin/env nextflow

process seuratRead10X {
  container "$params.sctools"
  errorStrategy 'retry'
  label 'mid_mem'
  cache false
  publishDir "$params.output.folder/seurat/object", mode : "copy"
  input:
    each path(count_path), path(vdj_path)
    path script_path
  output:
    path "*.h5Seurat", emit: seurat_object
  script:
  def sample_name = count_path.getSimpleName().replaceAll(/_counts/, "")
  def seurat_extension = count_path.getExtension()
  if(seurat_extension == 'h5Seurat')
    """
    export pdir=\$PWD
    echo \$pdir
    tar -zxvf ${vdj_path}
    Rscript scripts/seuratReadData.R -m ${count_path} -n ${sample_name}
    """
  else
    """
    export pdir=\$PWD
    echo \$pdir
    tar -zxvf ${count_path}
    tar -zxvf ${vdj_path}
    Rscript scripts/seuratReadData.R -m ./Counts/ -n ${sample_name}
    """
}


process seuratMergeH5Obj {
  container "$params.sctools"
  errorStrategy 'retry'
  label 'mid_mem'
  publishDir "$params.output.folder/seurat/object", mode : "copy"
  input:
    path seurat_objects
    val study_id
    path script_path
  output:
    path "${study_id}_merged_object.h5Seurat", emit: merged_object
  script:
    """
    Rscript scripts/seuratMergeData.R -n ${study_id}
    """
}