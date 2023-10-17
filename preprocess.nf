#!/usr/bin/env nextflow

process preProcess10XGEX {
  module "R/4.2.2-foss-2022b"
  errorStrategy 'retry'
  label 'high_mem'
  cache true
  publishDir "$params.output.folder/preprocess/object/sce", mode : "copy", pattern : "*.rda"
  publishDir "$params.output.folder/preprocess/figures/${sample_name}", mode : "copy", pattern : "*.pdf"
  input:
    each path(count_path)
    path script_path
    val mode
  output:
    path "${sample_name}.rda", emit: singlecell_object
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
  #tar -zxvf ${count_path}
  Rscript scripts/preProcessData.R -m ${count_path} -n ${sample_name}
  """
}

process convertFormats {
  module "R/4.2.2-foss-2022b"
  errorStrategy 'retry'
  label 'high_mem'
  cache true
  publishDir "$params.output.folder/preprocess/object/$params.meta.format", mode : "copy", pattern : "*.h5Seurat"
  input:
    each path(count_path)
  output:
    path "${sample_name}.h5Seurat", emit: singlecell_object
    //path "*.pdf", emit: singlecell_plot
  script:
  sample_name = count_path.getSimpleName()
  """
  Rscript $params.input.script_path/formatConverter.R -c ${count_path} -f SCE \
    -t $params.meta.format -n ${sample_name} 1> stdout.log 2> stderr.log
  """
}

process integrateTCRs {
  module "R/4.2.2-foss-2022b"
  errorStrategy 'retry'
  label 'high_mem'
  cache true
  publishDir "$params.output.folder/preprocess/object/TCRmerged", mode : "copy", pattern : "*.h5Seurat"
  input:
    each path(count_path)
    path vdj_path
  output:
    path "${sample_name}_Tmerged.h5Seurat", emit: singlecell_object
    //path "*.pdf", emit: singlecell_plot
  script:
  sample_name = count_path.getSimpleName()
  """
  Rscript $params.input.script_path/mergeTCRData.R -c ${count_path} \
    -t ${sample_name}_contig_annotations.csv -n ${sample_name} 1> stdout.log 2> stderr.log
  """
}
