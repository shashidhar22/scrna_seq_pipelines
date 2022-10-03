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
# Load in Seurat object and run SCTranform to scale and normalize each samples
readCountObject <- function(h5_path) {
    h5_object <- SeuratDisk::LoadH5Seurat(h5_path)
    h5_object <- subset(x = h5_object, subset = discard == FALSE)
    h5_object <- Seurat::SCTransform(h5_object, assay = "RNA", 
        vars.to.regress = "subsets_Mito_percent", verbose = FALSE)
    return(h5_object)
}
    # Scale data and run PCA on each seurat object
# Note: When integrating public datasets, this step can fail on low quality
# runs. Adjusting the npcs can be one solution in such cases
runPCA <- function(h5_object, features) {
    h5_object <- ScaleData(h5_object, features = features, verbose = FALSE)
    h5_object <- RunPCA(h5_object, features = features, verbose = FALSE, npcs = 30)
    return(h5_object)
}
# Set SCT as the default assay for integration
setDefaultAssay <- function(h5_object) {
    DefaultAssay(h5_object) <- "SCT"
    return(h5_object)
}

h5_list <- base::list.files(path = "./",
        pattern = "h5Seurat",
        full.names = TRUE) %>% 
    purrr::map(readCountObject) 

features <- Seurat::SelectIntegrationFeatures(object.list = h5_list)
h5_list <- h5_list %>% 
    purrr::map(~runPCA(.x, features))

anchors <- Seurat::FindIntegrationAnchors(object.list = h5_list, 
    anchor.features = features, reduction = "rpca", l2.norm = FALSE)

intergrated_object <- Seurat::IntegrateData(anchorset = anchors)
intergrated_object <- Seurat::ScaleData(intergrated_object, verbose = FALSE)
intergrated_object <- Seurat::RunPCA(intergrated_object, verbose = FALSE)
intergrated_object <- Seurat::RunUMAP(intergrated_object, dims = 1:30, 
    reduction = "pca")

output_path <- base::paste(parser$study_id, "merged_object.h5Seurat", sep = "_")
SeuratDisk::SaveH5Seurat(object = intergrated_object, 
    filename = output_path, overwrite = TRUE)


