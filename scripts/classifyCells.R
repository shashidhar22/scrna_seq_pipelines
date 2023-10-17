library(optparse)
library(tidyverse)
library(monocle3)
library(garnett)
# Accept arguments
option_list <- list(
    make_option(c("-c", "--count_path"), 
                type="character",  
                help="Path to count file",
                dest="count_path"),
    make_option(c("-m", "--marker_path"), 
            type="character",  
            help="Path to marker file",
            dest="from")
)
parser <- parse_args(OptionParser(option_list=option_list), print_help_and_exit = TRUE)
count_object <- SeuratDisk::LoadH5Seurat(count_path)
metadata <- read_csv("../../metadata/KS_ICMH_metadata.csv")

icmh_samples <- metadata %>% 
    filter(str_detect(locus, "5primeGEX")) %>% 
    mutate(orig_sample = str_remove(repoID, "GEX_|_GEX|_GE|_KS_"), orig_sample = str_replace_all(orig_sample, "-", "_"), orig_sample = str_remove(orig_sample, "_KS_")) %>% 
    pull(orig_sample)

icmh_seurat_pbmc <- subset(count_object, subset = Sample %in% icmh_samples & !str_detect(Sample, "scREP|SC")) 
icmh_seurat_screp <- subset(count_object, subset = Sample %in% icmh_samples & str_detect(Sample, "scREP|SC"))
screp_batch <- metadata %>% 
    dplyr::filter(str_detect(locus, "5primeGEX")) %>% 
    mutate(orig_sample = str_remove(repoID, "GEX_|_GEX|_GE|_KS_"), orig_sample = str_replace_all(orig_sample, "-", "_"), orig_sample = str_remove(orig_sample, "_KS_")) %>% 
    dplyr::filter(str_detect(orig_sample, "scREP|SC")) %>% pull(Batch)1

names(screp_batch) <- metadata %>% 
    dplyr::filter(str_detect(locus, "5primeGEX")) %>% 
    mutate(orig_sample = str_remove(repoID, "GEX_|_GEX|_GE|_KS_"), orig_sample = str_replace_all(orig_sample, "-", "_"), orig_sample = str_remove(orig_sample, "_KS_")) %>%  
    dplyr::filter(str_detect(orig_sample, "scREP|SC")) %>% pull(orig_sample)

colData(icmh_cds_screp)$Batch <- screp_batch[colData(icmh_cds_screp)$Sample]

icmh_cds_screp <- as.cell_data_set(icmh_seurat_screp)
colData(icmh_cds_screp)$Batch <- screp_batch[colData(icmh_cds_screp)$Sample]
icmh_cds_screp <- estimate_size_factors(icmh_cds_screp)
icmh_cds_screp <- preprocess_cds(icmh_cds_screp, num_dim = 100)
icmh_cds_screp <- align_cds(icmh_cds_screp, num_dim = 100, alignment_group = "Batch")
icmh_cds_screp <- reduce_dimension(icmh_cds_screp, reduction_method="UMAP")

icmh_cds_screp <- cluster_cells(icmh_cds_screp, resolution=1e-5)
screp_classifier <- train_cell_classifier(cds = icmh_cds_screp,
    marker_file = "garnettModel.yaml",
    db=org.Hs.eg.db,
    cds_gene_id_type = "ENSEMBL",
    num_unknown = 500,
    marker_file_gene_id_type = "SYMBOL")




icmh_cds_screp <- classify_cells(icmh_cds_screp, screp_classifier, db = org.Hs.eg.db, cluster_extend = TRUE, cds_gene_id_type = "ENSEMBL")
icmh_cds_screp <- reduce_dimension(icmh_cds_screp, reduction_method="UMAP", preprocess_method="PCA")
icmh_cds_screp <- reduce_dimension(icmh_cds_screp, reduction_method="tSNE", preprocess_method="PCA")
icmh_cds_scsubset <- icmh_cds_screp[, pData(icmh_cds_screp)$cluster_ext_type != "Unknown"]
icmh_cds_vdj <- icmh_cds_scsubset

icmh_cds_screp_plot <- plot_cells(icmh_cds_screp, color_cells_by="partition", group_cells_by="partition")
ggsave("KS_scREP_clusters.png", icmh_cds_screp_plot, height = 10, width = 10)

icmh_cds_screp_plot <- plot_cells(icmh_cds_scsubset, color_cells_by="cluster_ext_type", group_cells_by="partition")
ggsave("KS_scREP_celltypes.png", icmh_cds_screp_plot, height = 10, width = 10)

readAIRR <- function(file_path) {
    repertoire_id <- str_remove(basename(file_path), "_airr.tsv")
    airr_table <- read_tsv(file_path) %>% 
        mutate(repertoire_id = repertoire_id)
    return(airr_table)
}

vdj_table <- list.files("airr/", pattern = "_airr.tsv", full.names = TRUE, recursive = TRUE) %>% 
    map(readAIRR) %>% 
    bind_rows()



vdj_map <- vdj_table %>%  mutate(sample_name = str_c(cell_id, repertoire_id, sep = "_")) %>% dplyr::select(sample_name, junction_aa) %>% pull(junction_aa)
names(vdj_map) <- vdj_table  %>% mutate(sample_name = str_c(cell_id, repertoire_id, sep = "_")) %>% dplyr::select(sample_name, junction_aa) %>% pull(sample_name)


vdj_counts <- vdj_table %>%  mutate(sample_name = str_c(cell_id, repertoire_id, sep = "_")) %>% dplyr::select(sample_name, duplicate_count) %>% pull(duplicate_count)
names(vdj_counts) <- vdj_table  %>% mutate(sample_name = str_c(cell_id, repertoire_id, sep = "_")) %>% dplyr::select(sample_name, junction_aa) %>% pull(sample_name)

pData(icmh_cds_vdj)$junction_aa <- vdj_map[str_c(pData(icmh_cds_vdj)$Barcode, pData(icmh_cds_vdj)$Sample, sep = "_")]
pData(icmh_cds_vdj)$duplicate_count <- vdj_counts[str_c(pData(icmh_cds_vdj)$Barcode, pData(icmh_cds_vdj)$Sample, sep = "_")]


cpal <- c("B cells" = '#a6cee3', "Macrophage" ='#1f78b4', "NK Cell" = '#b2df8a',
    "T Cells" = '#33a02c', "Tcm" = '#fb9a99', "Teff" = '#e31a1c', 
    "Tem" = '#fdbf6f', "Th1" = '#ff7f00', "Th2" = '#cab2d6',
    "Trm" = '#6a3d9a', "Tscm" = '#ffff99')

icmh_cds_screp_plot <- plot_cells(icmh_cds_scsubset, color_cells_by="cluster_ext_type", group_cells_by="partition", label_cell_groups=FALSE) + theme(legend.position =  "bottom") + 
    scale_color_manual(values = cpal)
ggsave("KS_scREP_celltypes.pdf", icmh_cds_screp_plot, height = 10, width = 10)

icmh_cds_vdj_subset <- icmh_cds_vdj[, !is.na(pData(icmh_cds_vdj)$junction_aa)]

icmh_tcell_screp_plot <- plot_cells(icmh_cds_vdj_subset, color_cells_by="Sample", group_cells_by="partition", cell_size = 1) + theme(legend.position =  "bottom") + xlim(-9.7652, 9.9598) + ylim(-10.6565, 7.6317) 
ggsave("KS_scREP_tcells.pdf", icmh_tcell_screp_plot, height = 10, width = 10)
