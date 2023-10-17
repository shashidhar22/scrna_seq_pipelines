#!/usr/bin/env nextflow

/*
========================================================================================
                         Warren lab Cellranger pipeline
========================================================================================
 Warren lab Cellranger pipeline
----------------------------------------------------------------------------------------
*/
nextflow.enable.dsl=2

def helpMessage() {
    log.info"""
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run cellranger.nf -c '10Xcellranger.config'
    The it is recommended to use a config file to run this pipeline, as mentioned in the README. 
    Alternatively, the following parameters need to be provided to run this workflow.
    Alternate arguments: 
      --metadata                       Path to metadata file for single cell cellranger run.
      --gex_reference                  Path to cellranger GEX reference folder
      --vdj_reference                  Path to cellranger VDJ reference folder
      --gex_fastq_paths                Path to 10X GEX FASTQ files obtained from cellranger mkfastq, demux, or Illumina bcl2fastq
      --vdj_fastq_paths                Path to 10X VDJ FASTQ files obtained from cellranger mkfastq, demux, or Illumina bcl2fastq
      --study_id                       Unique study ID for cellranger run
      --mode                           Cellranger aggr normalization mode [mapped, none]
      --run_gex                        Boolean value to indicate whether cellranger count and aggr should be run on the dataset
      --run_vdj                        Boolean value to indicate whether cellranger VDJ should be run on the dataset
      --fastq_type                     Type of demultiplexing used to generate the FASTQ reads ['mkfastq', 'bcl2fastq', 'demux'] 
      --output_folder                  Path to output folder to store results
    """.stripIndent()
}

//Print help and exit
if (params.help){
    helpMessage()
    exit 0
}

process runCellranger {
  container "$params.sctools"
  label 'high_mem'
  publishDir "$params.output.folder/summaries/", mode : "copy", pattern : "*.html"
  publishDir "$params.output.folder/matrices/", mode : "copy", pattern : "*.h5"
  publishDir "$params.output.folder/${meta_row.Batch}/", mode : "copy", pattern : "${meta_row.sampleName}"
  input:
    each meta_row
    path fastq_path
    path gex_reference
    path vdj_reference
    val study_id
  output:
    path "${meta_row.sampleName}", type: dir, emit: cellranger_out
    path "*.html"
    path "*.h5"
  script:
    sample = meta_row.sampleName
    library = meta_row.repoID
    batch = meta_row.Batch
    if ("$meta_row.locus" == "5primeGEX") 
      """
      cellranger count --id=${sample} --transcriptome=${gex_reference} \
        --fastqs=${batch}/outs/fastq_path/${batch} \
        --sample=${library} --expect-cells=$meta_row.expected_cells \
        --localcores=6 --localmem=120 --chemistry=SC5P-R2

      #tar -czf ${sample}_counts.tar.gz ./${sample}/*  
      cp ${sample}/outs/web_summary.html ${sample}.html 
      cp ${sample}/outs/filtered_feature_bc_matrix.h5 ${sample}_filterd.h5
      """
    else if ("$meta_row.locus" == "5primeVDJ" )
      """
      cellranger vdj --id=$meta_row.sampleName --reference=${vdj_reference} \
        --fastqs=$meta_row.Batch/outs/fastq_path/$meta_row.Batch \
        --sample=$meta_row.repoName --localcores=6 --localmem=120

      #tar -czf "$meta_row.sampleName"_vdj.tar.gz ./$meta_row.sampleName/*  
      """
}




process prepCellranger {
  module "R"
  label 'low_mem'
  errorStrategy 'retry'
  cache true
  publishDir "$params.output.folder/config/multi_config/", mode : "copy", pattern : "*_config.csv"
  publishDir "$params.output.folder/config/", mode : "copy", pattern : "multi_config.tsv"
  input:
    path metadata
    path gex_reference
    path vdj_reference
    path fastq_path
    val study_id
    path script_path
  output:
    path "multi_config.tsv", emit: multi_config
    path "*_config.csv", emit: config_files
  script:
  """
  /app/software/R/4.2.2-foss-2021b/bin/Rscript ${script_path}/prepCellranger.R -m ${metadata.join(',')} -g ${gex_reference} -v ${vdj_reference} -s ${study_id} -f ${fastq_path.join(',')} 
  """
}

process cellrangerMulti {
  container "$params.sctools"
  label 'high_mem'
  errorStrategy 'retry'
  publishDir "$params.output.folder/cellranger/gex/", mode : "copy", pattern : "${multi_sample.sample_name}"
  publishDir "$params.output.folder/cellranger/vdj/", mode : "copy", pattern : "${multi_sample.sample_name}"
  publishDir "$params.output.folder/cellranger/summaries/", mode : "copy", pattern : "${multi_sample.sample_name}.html"
  publishDir "$params.output.folder/cellranger/vdj/airr/$multi_sample.sample_name/", mode : "copy", pattern : "${multi_sample.sample_name}_airr.tsv"
  publishDir "$params.output.folder/cellranger/vdj/screpertoire/$multi_sample.sample_name/", mode : "copy", pattern : "${multi_sample.sample_name}_contig_annotations.csv"
  publishDir "$params.output.folder/cellranger/vdj/contigs/$multi_sample.sample_name/", mode : "copy", pattern : "${multi_sample.sample_name}_consensus.fa"
  input:
    path config_files
    each multi_sample
    path gex_reference
    path vdj_reference
    path fastq_path
    path script_path
  output:
    val "${sample}", emit: sample_list
    path "${multi_sample.sample_name}", type: 'dir', emit: count_list
    path "${multi_sample.sample_name}", type: 'dir', emit: vdj_list
    path "${multi_sample.sample_name}_config.csv", emit: multi_config_files
    path "${multi_sample.sample_name}.html", emit: summary_files
    path "${multi_sample.sample_name}_contig_annotations.csv", emit: tcr_file
    path "${multi_sample.sample_name}_airr.tsv", emit: airr_file
    path "${multi_sample.sample_name}_consensus.fa", emit: tcr_contigs

  script:
  def memory = "$task.memory" =~  /\d+/
  def sample = "$multi_sample.sample_name"
  """
  cellranger multi --id=${sample} --csv=$multi_sample.config_file  --localcores 12 --localmem 300 1> ${sample}_cellranger.log  2>&1
  cp -r ${sample}/outs/per_sample_outs/${sample}/count Counts
  #tar -czf ${sample}_counts.tar.gz ./Counts/*
  mkdir VDJ
  cp -r ${sample}/outs/per_sample_outs/${sample}/vdj_*/* VDJ/
  #tar -czf ${sample}_vdj.tar.gz ./VDJ/*
  cp ${sample}/outs/per_sample_outs/${sample}/web_summary.html ${sample}.html 
  cp ${sample}/outs/per_sample_outs/${sample}/vdj_t/filtered_contig_annotations.csv ${sample}_contig_annotations.csv
  cp ${sample}/outs/per_sample_outs/${sample}/vdj_t/airr_rearrangement.tsv ${sample}_airr.tsv
  cp ${sample}/outs/per_sample_outs/${sample}/vdj_t/consensus.fasta ${sample}_consensus.fa
  """

}

workflow {
  metadata = Channel.fromPath(params.input.metadata)
  gex_reference = Channel.fromPath(params.input.gex_reference)
  vdj_reference = Channel.fromPath(params.input.vdj_reference)
  fastq_paths = Channel.fromPath(params.input.fastq_path)
  study_id = Channel.from(params.input.study_id)

  main:
    prepCellranger(metadata.collect{"$it"}, gex_reference, vdj_reference, fastq_paths.collect{"$it"}, study_id)
    cellrangerMulti(prepCellranger.out.config_files, prepCellranger.out.multi_config.splitCsv(header: true, quote: '\"', sep: '\t'), gex_reference, vdj_reference, fastq_paths)
}
/*workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()

    sendMail(to: $params.input.email , subject: 'Cellranger pipeline execution: ${workflow.success}', body: msg)
}*/