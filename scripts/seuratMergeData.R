library(optparse)
library(tidyverse)
library(Seurat)     
library(SeuratDisk)
# Accept arguments
option_list <- list(
    make_option(c("-n", "--study_id"), 
                type="character",  
                help="Study ID",
                dest="study_id")
)
parser <- parse_args(OptionParser(option_list=option_list), print_help_and_exit = TRUE)

readCountObject <- function(h5_path) {
    print(h5_path)
    h5_object <- SeuratDisk::LoadH5Seurat(h5_path)
    h5_object <- NormalizeData(h5_object)
    h5_object <- FindVariableFeatures(h5_object, selection.method = "vst", nfeatures = 2000)
    return(h5_object)
}

runPCA <- function(h5_object, features) {
    print(h5_object)
    h5_object <- Seurat::ScaleData(h5_object, features = features, 
        verbose = FALSE)
    h5_object <- Seurat::RunPCA(h5_object, features = features, verbose = FALSE)
    return(h5_object)
}   

h5_list <- base::list.files(path = "./",
        pattern = "h5Seurat",
        full.names = TRUE) %>% 
    purrr::map(readCountObject) 

features <- Seurat::SelectIntegrationFeatures(object.list = h5_list)
h5_list <- h5_list %>% 
    purrr::map(~runPCA(.x, features = features ))

anchors <- Seurat::FindIntegrationAnchors(object.list = h5_list, anchor.features = features, 
    reference = c(1,2), reduction = "rpca", dims = 1:50)

intergrated_object <- Seurat::IntegrateData(anchorset = anchors, dims = 1:50)
intergrated_object <- Seurat::ScaleData(intergrated_object, verbose = FALSE)
intergrated_object <- Seurat::RunPCA(intergrated_object, verbose = FALSE)
intergrated_object <- Seurat::RunUMAP(intergrated_object, dims = 1:50, reduction = "pca")

output_path <- base::paste(parser$study_id, "merged_object.h5Seurat", sep = "_")
SeuratDisk::SaveH5Seurat(object = intergrated_object, filename = output_path, overwrite = TRUE)


