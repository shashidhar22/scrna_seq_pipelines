nextflow.enable.dsl=2
include {runCellranger ; prepCellranger ; cellrangerMulti} from './cellranger'
include {seuratRead10X ; seuratMergeH5Obj} from './seurat'
include {preProcess10XGEX} from './preprocess'

workflow CellRangerMulti {
  metadata = Channel.fromPath(params.input.metadata)
  gex_reference = Channel.fromPath(params.input.gex_reference)
  vdj_reference = Channel.fromPath(params.input.vdj_reference)
  fastq_paths = Channel.fromPath(params.input.fastq_path)
  study_id = Channel.from(params.input.study_id)
  script_path = Channel.fromPath(params.input.script_path)
  mode = Channel.from(params.input.process_mode)
  main:
    prepCellranger(metadata.collect{"$it"}, gex_reference, vdj_reference, fastq_paths.collect{"$it"}, study_id, script_path)
    cellrangerMulti(prepCellranger.out.config_files, prepCellranger.out.multi_config.splitCsv(header: true, quote: '\"', sep: '\t'), gex_reference, vdj_reference, fastq_paths.collect{"$it"}, script_path)
    preProcess10XGEX(cellrangerMulti.out.count_list, script_path, mode)
}

workflow PreprocessSCDataset {
  count_list = Channel.fromPath(params.input.count_list)
  script_path = Channel.fromPath(params.input.script_path)
  mode = Channel.from(params.input.process_mode)
  study_id = Channel.from(params.input.study_id)
  main:
    preProcess10XGEX(count_list, script_path, mode)
    seuratMergeH5Obj(preProcess10XGEX.out.singlecell_object.toSortedList(), study_id, script_path)

}


workflow KS10X011622 {
  metadata = Channel.fromPath(params.input.metadata)
  gex_reference = Channel.fromPath(params.input.gex_reference)
  vdj_reference = Channel.fromPath(params.input.vdj_reference)
  fastq_paths = Channel.fromPath(params.input.fastq_path)
  study_id = Channel.from(params.input.study_id)
  script_path = Channel.fromPath(params.input.script_path)
  main:
    prepCellranger(metadata.collect{"$it"}, gex_reference, vdj_reference, fastq_paths.collect{"$it"}, study_id, script_path)
    cellrangerMulti(prepCellranger.out.config_files, prepCellranger.out.multi_config.splitCsv(header: true, quote: '\"', sep: '\t'), gex_reference, vdj_reference, fastq_paths, script_path)
    seuratRead10X(cellrangerMulti.out.count_list.merge(cellrangerMulti.out.vdj_list), script_path)
    seuratMergeH5Obj(seuratRead10X.out.seurat_object.toSortedList().collect(), study_id, script_path)
}

workflow KS10X012622_HHV8 {
  metadata = Channel.fromPath(params.input.metadata)
  gex_reference = Channel.fromPath(params.input.gex_reference)
  vdj_reference = Channel.fromPath(params.input.vdj_reference)
  fastq_paths = Channel.fromPath(params.input.fastq_path)
  study_id = Channel.from(params.input.study_id)
  script_path = Channel.fromPath(params.input.script_path)
  main:
    runCellranger(metadata.splitCsv(header: true, quote: '\"', sep: ',').filter{ it.locus == "5primeGEX"}, fastq_paths.collect{"$it"}, gex_reference, vdj_reference, study_id, script_path)

}

workflow RC032822_merge {
  h5paths = Channel.fromPath(params.input.h5paths)
  study_id = Channel.from(params.input.study_id)
  script_path = Channel.fromPath(params.input.script_path)

  main:
    seuratMergeH5Obj(h5paths.toSortedList(), study_id, script_path)
}
