library(Seurat)  
library(SeuratDisk)   
library(optparse)
library(tidyverse)
library(SingleCellExperiment)
library(scuttle)
library(scrna)
library(scater)
library(uwot)
library(rtracklayer)
library(DropletUtils)
# Accept arguments
option_list <- list(
    make_option(c("-m", "--count_path"), 
        type="character",  
        help="Path to 10X count data",
        dest="count_path"),
    make_option(c("-n", "--name"),
        type="character",
        help="Sample name",
        dest="sample_name")
)
parser <- parse_args(OptionParser(option_list=option_list), print_help_and_exit = TRUE)
## Load Seurat object
if (str_detect(parser$count_path, "h5Seurat")) {
    count_object <- SeuratDisk::LoadH5Seurat(parser$count_path)
    sc_object <- DropletUtils::read10xCounts(matrix_path)
} else {
    count_path <- paste(parser$count_path, "/sample_feature_bc_matrix", sep = "")
    count_data <- Seurat::Read10X(count_path)
    count_object <- Seurat::CreateSeuratObject(counts = count_data, min.cells = 100, min.features = 500)

}
count_object$sample_name <- parser$sample_name 
ptid <- str_extract(parser$sample_name, "^\\d+_\\d+|^\\w+_\\d+|^\\w+")
count_object$patient_id <- ptid
count_object <- Seurat::PercentageFeatureSet(count_object, 
    pattern = "^MT-", col.name = "percent.mt")
count_object <- Seurat::SCTransform(count_object, 
    vars.to.regress = "percent.mt", verbose = FALSE)
count_object <- Seurat::RunPCA(count_object, verbose = FALSE)
count_object <- Seurat::RunUMAP(count_object, dims = 1:50, verbose = FALSE)
count_object <- Seurat::FindNeighbors(count_object, dims = 1:50, verbose = FALSE)
count_object <- Seurat::FindClusters(count_object, dims = 1:50, verbose = FALSE)

## Load VDJ data
vdj_path <- paste(parser$vdj_path, "/filtered_contig_annotations.csv", sep = "")
vdj_data <- read.csv(vdj_path, header = TRUE)
vdj_object <- scRepertoire::combineTCR(vdj_data, sample = parser$sample_name,
    filterMulti = TRUE)
vdj_object <- vdj_object[[1]] %>% 
    as_tibble() %>% 
    separate(TCR1, into = c("TRAV", "TRAJ", "TRAC"), sep = "\\.", remove = FALSE) %>%
    separate(TCR2, into = c("TRBV", "TRBD", "TRBJ", "TRBC"), sep = "\\.", remove = FALSE)

## Merge Seurat and VDJ objects
count_object$barcodes <- rownames(count_object@meta.data)
cell_names <- count_object@meta.data %>% 
    as_tibble() %>% 
    mutate(barcodes = str_remove(barcodes, "_\\d+$"), 
        barcodes = str_c(sample_name, barcodes, sep = "_")) %>% 
    pull(barcodes) %>% 
    as.vector()
rownames(count_object@meta.data) <- cell_names
count_object <- RenameCells(count_object, new.names = cell_names)
count_object <- combineExpression(vdj_object, count_object, 
    cloneCall = "strict", proportion = TRUE)

out_path <- paste(parser$sample_name, "h5Seurat", sep = ".")
SeuratDisk::SaveH5Seurat(count_object, out_path, overwrite = TRUE)
