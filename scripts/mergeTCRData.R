library(optparse)
library(tidyverse)
library(Seurat)     
library(SeuratDisk)
library(SingleCellExperiment)
#library(monocle3)
library(optparse)
library(tidyverse)
library(scRepertoire)
# Accept arguments
option_list <- list(
    make_option(c("-c", "--count_path"), 
                type="character",  
                help="Path to count file",
                dest="count_path"),
    make_option(c("-t", "--tcr_path"), 
                type="character",  
                help="Path to TCR data",
                dest="tcr_path"),    
    make_option(c("-n", "--sample_name"), 
            type="character",  
            help="Sample name",
            dest="sample_name")        
)
parser <- parse_args(OptionParser(option_list=option_list), print_help_and_exit = TRUE)

# Load and combine TCR data
tcr_table <- readr::read_csv(parser$tcr_path)
sample_name <- str_remove(parser$sample_name, "_GEX")
sample_id <- str_extract(parser$sample_name, "GEX")
tcr_table <- combineTCR(list(tcr_table), samples = c(sample_name),
    ID = c(sample_name), cells = "T-AB", removeNA =TRUE, filterMulti = TRUE)

# Load and merge gene expression data with the TCR data
gex_object <- LoadH5Seurat(parser$count_path)
gex_colnames <- gex_object@meta.data |>
    as_tibble() |>
    mutate(screp_barcode = str_c(Sample, Barcode, sep = "_")) |>
    pull(screp_barcode)
rownames(gex_object@meta.data) <- gex_colnames
gex_object <- RenameCells(gex_object, new.names = gex_colnames)
gex_object <- combineExpression(tcr_table, gex_object,
    cloneCall = "gene", 
    proportion = FALSE, 
    cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

# Define output file name and store result
output_path <- str_c(parser$sample_name, "_Tmerged.h5Seurat", sep = "")
SaveH5Seurat(gex_object, output_path)