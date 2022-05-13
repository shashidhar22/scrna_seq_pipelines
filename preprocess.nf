#!/usr/bin/env nextflow

process preProcess10XGEX {
  container "$params.sctools"
  errorStrategy 'retry'
  label 'high_mem'
  cache false
  publishDir "$params.output.folder/preprocess/object", mode : "copy", pattern : "h5Seurat"
  publishDir "$params.output.folder/preprocess/figures", mode : "copy", pattern : "pdf"
  input:
    each path(count_path)
    path script_path
  output:
    path "*.h5Seurat", emit: singlecell_object
    path "*.pdf", emit: singlecell_plot
  script:
  def sample_name = count_path.getSimpleName().replaceAll(/_counts/, "")
  """
  export pdir=\$PWD
  echo \$pdir
  tar -zxvf ${count_path}
  Rscript scripts/preProcessData.R -m ./Counts/ -n ${sample_name}
  """
}
