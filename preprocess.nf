#!/usr/bin/env nextflow

process preProcess10XGEX {
  container "$params.sctools"
  errorStrategy 'retry'
  label 'high_mem'
  cache true
  publishDir "$params.output.folder/preprocess/object", mode : "copy", pattern : "*.h5Seurat"
  publishDir "$params.output.folder/preprocess/figures/${sample_name}", mode : "copy", pattern : "*.pdf"
  input:
    each path(count_path)
    path script_path
    val mode
  output:
    path "${sample_name}.h5Seurat", emit: singlecell_object
    path "*.pdf", emit: singlecell_plot
  script:
  sample_name = count_path.getSimpleName().replaceAll(/_counts/, "")
  if (mode == "h5")
  """
  export pdir=\$PWD
  echo \$pdir
  export TMPDIR=$workDir
  Rscript scripts/preProcessData.R -m ${count_path} -n ${sample_name}
  """
  else if (mode == "10X")
  """
  export pdir=\$PWD
  echo \$pdir
  export TMPDIR=$workDir
  tar -zxvf ${count_path}
  Rscript scripts/preProcessData.R -m ./Counts/ -n ${sample_name}
  """
}
