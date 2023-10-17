library(optparse)
library(tidyverse)
library(Seurat)     
library(SeuratDisk)
library(SingleCellExperiment)
#library(monocle3)
library(optparse)
library(tidyverse)
# Accept arguments
option_list <- list(
    make_option(c("-c", "--count_path"), 
                type="character",  
                help="Path to count file",
                dest="count_path"),
    make_option(c("-f", "--from"), 
            type="character",  
            help="Format of input file",
            dest="from"),
    make_option(c("-t", "--to"), 
            type="character",  
            help="Format of output file",
            dest="to"),
    make_option(c("-n", "--sample_name"), 
            type="character",  
            help="Sample name",
            dest="sample_name")        
)
parser <- parse_args(OptionParser(option_list=option_list), print_help_and_exit = TRUE)

if (parser$from == "SCE" & parser$to == "Seurat") {
    load(parser$count_path)
    seu_object <- SeuratObject::as.Seurat(sc_object, counts = "counts", 
        data = "logcounts") 
    seu_object <- RenameAssays(seu_object, originalexp = 'RNA')
    seu_object$orig.ident <- seu_object$Sample
    Idents(seu_object) <- "orig.idents"
    output_path <- stringr::str_c(parser$sample_name, "h5Seurat", sep = ".")
    SaveH5Seurat(seu_object, output_path)
}

if (parser$from == "SCE" & parser$to == "monocle3") {
    load(parser$count_path)
    rowData(sc_object)$gene_short_name <- rowData(sc_object)$Symbol
    cds_object <- new_cell_data_set(counts(sc_object), 
        cell_metadata = colData(sc_object), gene_metadata = rowData(sc_object))
    cds_object <- preprocess_cds(cds_object, num_dim = 100)
    cds_object <- reduce_dimension(cds_object, reduction_method = "UMAP")
    cds_object <- reduce_dimension(cds_object, reduction_method = "tSNE")
    save(cds_object, file = str_c(file_name, "_cds.rda", sep =""))
}