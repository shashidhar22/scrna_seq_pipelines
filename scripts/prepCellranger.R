library(tidyverse)     
library(lubridate)
library(readxl)
library(optparse)
# Accept arguments
option_list <- list(
    make_option(c("-m", "--meta_path"), 
                type="character",  
                help="Path to metadata files. When multiple sudies are to be processed, pass a comma separated list of metadata files",
                dest="mlist"),
    make_option(c("-g", "--gref_path"), 
                type="character",  
                help="Path to GEX cellranger reference",
                dest="gref"),
    make_option(c("-v", "--vref_path"), 
                type="character",  
                help="Path to VDJ cellranger reference",
                dest="vref"),
    make_option(c("-s", "--study_id"), 
                type="character",  
                help="Study name, ouput folder will be named using this field",
                dest="sid"),
    make_option(c("-f", "--fastq_path"), 
                type="character",  
                help="Path to mkfastq output folders. When multiple sudies are to be processed, pass a comma separated list of metadata files",
                dest="flist")

)
parser <- parse_args(OptionParser(option_list=option_list), print_help_and_exit = TRUE)
meta_list <- strsplit(parser$mlist, ",")
fastq_list <- strsplit(parser$flist, ",")

mtable <- meta_list %>%
    purrr::map(~readr::read_csv(.x, na = c("", "NA"),)) %>%
    dplyr::bind_rows()

getPaths <- function(run_name) {
    fastq_path <- paste(run_name, "outs/fastq_path",run_name, sep = "/")
    multi_path <- paste(run_name, "outs/fastq_path", run_name,sep = "/")
    rtable <- list.dirs(fastq_path, recursive = FALSE, full.names = FALSE) %>% 
        tibble::as_tibble_col(column_name = "repoName") %>% 
        dplyr::mutate(run = run_name, 
            path = fastq_path, 
            fastqs = tools::file_path_as_absolute(fastq_path))
}


ftable <- fastq_list %>%
    unlist() %>% 
    purrr::map(getPaths) %>%
    bind_rows()

ctable <- dplyr::left_join(mtable, ftable, by = "repoName") %>%
    dplyr::mutate(feature_types = case_when(locus == "5primeGEX" & is.na(VDJType) ~ "Gene Expression",
        locus == "5primeVDJ" & VDJType == "T-cell" ~ "VDJ-T"),
        subsample_rate = "") %>%
    dplyr::select(repoName, sampleName, feature_types, fastqs, subsample_rate) %>%
    dplyr::rename(fastq_id = repoName)

prepMulti <- function(ctable, gref, vref, ecell) {
    sample <- ctable %>% pull(sampleName) %>% unique()
    ctable <- ctable %>%
        select(fastq_id, fastqs, feature_types, subsample_rate)
    out_path <- paste(sample,  "config.csv", sep = "_")
    file_conn <- file(out_path)
    line_vector <- c("[gene-expression]",
        paste("reference", tools::file_path_as_absolute(gref), sep=","),
        paste("chemistry", 'auto', sep=","), 
        paste("expect-cells", ecell, sep=","),
        "[vdj]",
        paste("reference", tools::file_path_as_absolute(vref), sep=","),
        "[libraries]",
        paste(colnames(ctable), collapse=","))
    for (i in seq_along(ctable$fastq_id)) {
        row <- ctable[i, ]
        line_vector <- c(line_vector, paste(unname(unlist((row))), collapse = ","))
    }
    writeLines(line_vector, file_conn)
    close(file_conn)    

    ftable <- tibble(sample_name = c(sample),
        config_file = c(out_path))
}

ftable <- ctable %>%
    group_by(sampleName) %>%
    group_split() %>%
    purrr::map(~prepMulti(.x, parser$gref, parser$vref, 3000)) %>%
    bind_rows()
write_tsv( ftable, "multi_config.tsv")
write_tsv( ctable, "run_data.tsv")
