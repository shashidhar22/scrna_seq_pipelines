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
  module 'CellRanger/6.1.1'
  label 'mid_mem'
  publishDir "$params.output.folder/summaries/", mode : "copy", pattern : "*.html"
  publishDir "$params.output.folder/matrices/", mode : "copy", pattern : "*.h5"
  publishDir "$params.output.folder/${meta_row.Batch}/${meta_row.sampleName}/", mode : "copy", pattern : "*.tar.gz"
  input:
    each meta_row
    path fastq_path
    path gex_reference
    path vdj_reference
    val study_id
  output:
    path "*.tar.gz", emit: cellranger_out
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

      tar -czf ${sample}_counts.tar.gz ./${sample}/*  
      cp ${sample}/outs/web_summary.html ${sample}.html 
      cp ${sample}/outs/filtered_feature_bc_matrix.h5 ${sample}_filterd.h5
      """
    else if ("$meta_row.locus" == "5primeVDJ" )
      """
      cellranger vdj --id=$meta_row.sampleName --reference=${vdj_reference} \
        --fastqs=$meta_row.Batch/outs/fastq_path/$meta_row.Batch \
        --sample=$meta_row.repoName --localcores=6 --localmem=120

      tar -czf "$meta_row.sampleName"_vdj.tar.gz ./$meta_row.sampleName/*  
      """
}




process prepCellranger {
  module 'R/4.1.0-foss-2020b'
  label 'low_mem'
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
  Rscript scripts/prepCellranger.R -m ${metadata.join(',')} -g ${gex_reference} -v ${vdj_reference} -s ${study_id} -f ${fastq_path.join(',')} 
  """
}

process cellrangerMulti {
  module 'R/4.1.0-foss-2020b'
  module 'CellRanger/6.1.1'
  label 'high_mem'
  publishDir "$params.output.folder/$multi_sample.sample_name/", mode : "copy", pattern : "*counts.tar.gz"
  publishDir "$params.output.folder/$multi_sample.sample_name/", mode : "copy", pattern : "*vdj.tar.gz"
  publishDir "$params.output.folder/$multi_sample.sample_name/", mode : "copy", pattern : "${multi_sample.sample_name}.html"
  input:
    path config_files
    each multi_sample
    path gex_reference
    path vdj_reference
    path fastq_path
    path script_path
  output:
    val "${sample}", emit: sample_list
    path "${multi_sample.sample_name}_counts.tar.gz", emit: count_list
    path "${multi_sample.sample_name}_vdj.tar.gz", emit: vdj_list
    path "${multi_sample.sample_name}_config.csv", emit: multi_config_files
    path "${multi_sample.sample_name}.html", emit: summary_files
  script:
  def memory = "$task.memory" =~  /\d+/
  def sample = "$multi_sample.sample_name"
  """
  cellranger multi --id=${sample} --csv=$multi_sample.config_file --localcores $task.cpus --localmem ${memory[0]}
  cp -r ${sample}/outs/per_sample_outs/${sample}/count Counts
  tar -czf ${sample}_counts.tar.gz ./Counts/*
  mkdir VDJ
  cp -r ${sample}/outs/per_sample_outs/${sample}/vdj_*/* VDJ/
  tar -czf ${sample}_vdj.tar.gz ./VDJ/*
  cp ${sample}/outs/per_sample_outs/${sample}/web_summary.html ${sample}.html 
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
/*process TENX_GEX_MAP {
  echo false
  label 'low_mem'
  publishDir "$params.output.folder/Metadata", mode : 'copy'
  module 'R/4.0.2-foss-2019b'
  
  input:
    path meta_file from gex_metadata
    val fastq_path from gex_fastq_paths
    val study from gex_study_id
    val map from map_gex

  output:
    path "GEX_h5_samplesheet.csv" into gex_h5sheet
    path "GEX_samplesheet.csv" into gex_samplesheet
  
  when:
    map == true
        
  script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    # Metadata file provide should have the 12 mandatory fields mentioned in the README 
    data <- readr::read_csv("${meta_file}")
    data <- data %>% 
            dplyr::filter(platform == '10XGenomics' & str_detect(locus, "GEX"))
    data <- data %>% 
            dplyr::mutate(library= "GEX")
    data <- data %>% 
            dplyr::mutate(fastqType = "$params.count.fastq_type")
    # Account for folder structure differences between mkfastq and demux
    if ("$params.count.fastq_type" == "mkfastq" | "$params.count.fastq_type" == "bcl2fastq") {
      data <- data %>% 
              dplyr::mutate(fastqSample = repoName,
                            fastqPath = "${fastq_path}")
    } else if ("$params.count.fastq_type" == "demux") {
      data <- data %>% 
              dplyr::mutate(fastqSample = paste(repoName, "demux_id", sep="/"), 
                            fastqFolder = "${fastq_path}",
                            fastqPath =  paste(fastqFolder, fastqSample, sep="/"))
    }
    data <- data %>% mutate(indices = paste("SI-GA-", indices, sep=""))
    # Split VDJ and GEX libraries if present
    data_gex <- data %>% filter(library == "GEX")
    data_gex <- data_gex %>% mutate(molecule_h5 = paste("$params.output.folder", "Counts", 
                                    sampleName, "outs/molecule_info.h5", sep="/"))
    # Create samplesheet for mkfastq, and H5 and analysis sheets for VDJ analysis
    data_gex_h5 <- data_gex %>% 
                   dplyr::rename(library_id = sampleName)
    data_gex_h5 <- data_gex_h5 %>% 
                   dplyr::select(library_id, molecule_h5, everything())
    data_gex_ss <- data_gex %>%  
                   dplyr::select(repoName, indices) %>% 
                   dplyr::rename(Sample = repoName, Index = indices) 
    data_gex_ss <- data_gex_ss %>% 
                   dplyr::mutate(Lane = '1-2') %>% 
                   dplyr::select(Lane, Sample, Index)
    # Write metadata files
    readr::write_csv(data_gex_h5, "GEX_h5_samplesheet.csv")
    readr::write_csv(data_gex_ss, "GEX_samplesheet.csv")
    """
}

process TENX_VDJ_MAP {
  echo false
  label 'low_mem'
  publishDir "$params.output.folder/Metadata", mode : 'copy'
  module 'R/4.0.2-foss-2019b'
  
  input:
    path meta_file from vdj_metadata
    val fastq_path from vdj_fastq_paths
    val study from vdj_study_id
    val map from map_vdj

  output:
    path "VDJ_samplesheet.csv" into vdj_samplesheet
    path "VDJ_analysis_samplesheet.csv" into vdj_analysissheet
  
  when:
    map == true
        
  script:
    """
    #!/usr/bin/env Rscript
    library(tidyverse)
    # Metadata file provide should have the 12 mandatory fields mentioned in the README 
    data <- readr::read_csv('${meta_file}')
    data <- data %>% 
            dplyr::filter(platform == '10XGenomics' & str_detect(locus, "VDJ"))
    data <- data %>% 
            dplyr::mutate(library = "VDJ", 
                          fastqType = "$params.count.fastq_type")
    # Account for folder structure differences between mkfastq and demux
    if ("$params.count.fastq_type" == "mkfastq" | "$params.count.fastq_type" == "bcl2fastq") {
      data <- data %>% 
              dplyr::mutate(fastqSample = repoName, 
                            fastqPath = "${fastq_path}")
    } else if ("$params.count.fastq_type" == "demux") {
      data <- data %>% 
              dplyr::mutate(fastqSample = paste(repoName, "demux_id", sep="/"), 
                            fastqFolder = "${fastq_path}", 
                            fastqPath =  paste(fastqFolder, fastqSample, sep="/"))
    }
    data <- data %>% 
            dplyr::mutate(indices = paste("SI-GA-", indices, sep=""))
    # Split VDJ and GEX libraries if present
    data_vdj <- data %>% 
                dplyr::filter(library == "VDJ")
    data_vdj <- data %>% 
                dplyr::mutate(library_id = paste(sampleName, VDJType, sep='_'))
    data_vdj <- data_vdj %>% 
                dplyr::mutate(vdj_sequences = paste("$params.output.folder", 
                                                    "VDJ", 
                                                    library_id, 
                                                    "outs/airr_rearrangement.tsv", 
                                                    sep="/"))
    # Create samplesheet for mkfastq, and H5 and analysis sheets for VDJ analysis
    data_vdj_ss <- data_vdj %>% 
                   dplyr::select(repoName, indices) %>% 
                   dplyr::rename(Sample = repoName, Index = indices) 
    data_vdj_ss <- data_vdj_ss %>% 
                   dplyr::mutate(Lane = '1-2') %>% 
                   dplyr::select(Lane, Sample, Index)
    # Write metadata files
    write_csv(data_vdj_ss, "VDJ_samplesheet.csv")
    write_csv(data_vdj, "VDJ_analysis_samplesheet.csv")
    """
}

//Copy GEX output channels into channels for counting and aggreagtion
gex_h5sheet.into {count_gex_h5sheet; aggr_gex_h5sheet}

process TENX_COUNT {
  echo false 
  publishDir "$params.output.folder/Counts" , mode : 'copy'
  label 'mid_mem'  
  module 'CellRanger/4.0.0'

  input:
    each sample from count_gex_h5sheet.splitCsv(header: true, quote: '\"')
    val run_count from count_gex
  
  output:
    path "${sample.library_id}" into count_path
    val task.exitStatus into count_status

  when:
    run_count == true

  script:
    memory = "$task.memory" =~  /\d+/
    if("$params.count.fastq_type" == "mkfastq" | "$params.count.fastq_type" == "bcl2fastq")
        """
        cellranger count --id=$sample.library_id --transcriptome=$params.input.gex_reference \
          --fastqs=$sample.fastqPath --sample=$sample.fastqSample --expect-cells=$sample.expected_cells \
          --localcores=$task.cpus --localmem=${memory[0]}
        """
    else if("$params.count.fastq_type" == "demux")
        """
        cellranger count --id=$sample.library_id --transcriptome=$params.input.gex_reference \
          --fastqs=$sample.fastqPath --expect-cells=$sample.expected_cells \
          --chemistry=$sample.chemistry --localcores=$task.cpus --localmem=${memory[0]}
        """
}

process TENX_AGGR {
  echo false
  publishDir "$params.output.folder/Counts" , mode : 'copy'
  label 'mid_mem'
  module 'CellRanger/4.0.0'

  input:
    path samplesheet from aggr_gex_h5sheet
    each mode from modes
    val status from count_status.collect()
    val run_aggr from aggr_gex
  
  output:
    path "Aggregate_${mode}_normalized" into aggr_path

  when:
    run_aggr == true

  script:
    """
    cellranger aggr --id=Aggregate_${mode}_normalized --csv=${samplesheet} --normalize=${mode} 
    """
}

process TENX_MATRIX {
  echo false
  publishDir "$params.output.folder/Counts/${aggr_out}/out/filtered_feature_bc_matrix" , mode : 'copy'
  module 'CellRanger/4.0.0'
  label 'low_mem'

  input:
    path aggr_out from aggr_path
    val run_matrix from mat_gex

  output:
    path "Filtered_expression_matrix.csv" into mtx_path

  when:
    run_matrix == true

  script:
    """
    cellranger mat2csv ${aggr_out}/outs/filtered_feature_bc_matrix Filtered_expression_matrix.csv
    """
}

process TENX_VDJ {
  echo false
  publishDir "$params.output.folder/VDJ" , mode : 'copy'
  module 'CellRanger/4.0.0'
  label 'mid_mem'

  input:
    each sample from vdj_analysissheet.splitCsv(header: true, quote: '\"') 
    val analyze from ana_vdj

  output:
    path "${sample.library_id}"

  when:
    analyze == true

  script:
    memory = "$task.memory" =~  /\d+/
    if("$params.count.fastq_type" == "mkfastq" | "$params.count.fastq_type" == "bcl2fastq")
      """
      cellranger vdj --id=$sample.library_id --reference=$params.input.vdj_reference --fastqs=$sample.fastqPath \
        --sample=$sample.fastqSample --localcores=$task.cpus --localmem=${memory[0]}
      """
    else if("$params.count.fastq_type" == "demux")
      """
      cellranger vdj --id=$sample.library_id --reference=$params.input.vdj_reference --fastqs=$sample.fastqPath \
        --localcores=$task.cpus --localmem=${memory[0]}
      """
}

*/
