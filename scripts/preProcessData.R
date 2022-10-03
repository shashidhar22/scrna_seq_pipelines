library(Seurat)  
library(SeuratDisk)   
library(optparse)
library(tidyverse)
library(SingleCellExperiment)
library(scuttle)
library(scran)
library(scater)
library(uwot)
library(rtracklayer)
library(DropletUtils)
library(GenomeInfoDb)
library(org.Hs.eg.db)
library(scDblFinder)
library(robustbase)
library(AnnotationHub)
# Set seed
set.seed(12357)
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
if (str_detect(parser$count_path, "h5Seurat")) {
    sc_object <- LoadH5Seurat(parser$count_path) 
    sc_object <- Seurat::as.SingleCellExperiment(sc_object)
} else {
    matrix_path <- base::paste(parser$count_path, "/sample_filtered_feature_bc_matrix", sep = "")
    sc_object <- DropletUtils::read10xCounts(matrix_path)
    sc_object$Sample <- parser$sample_name
}
if (ncol(rowData(sc_object)) != 0) {
    rownames(sc_object) <- rowData(sc_object)$ID
    rowData(sc_object)$GENEID <- rowData(sc_object)$ID
} else {
    nomenclature <- rownames(sc_object) %>% 
        head(1)
    annotation_hub <- AnnotationHub(cache = "./")
    anno_db <- query(annotation_hub, pattern = c("Homo Sapiens", "EnsDb"))[['AH104864']]
    if (!str_detect(nomenclature, "ENSG")) {
        genes <- rownames(sc_object)
        gene_map <- select(anno_db, genes, c("GENEID", "GENENAME"), "SYMBOL")
        gene_list <- gene_map$SYMBOL %>% 
            unique()
        sc_object <- sc_object[gene_list]
        ensembl_dict <- gene_map$GENEID
        names(ensembl_dict) <- gene_map$SYMBOL
        rowData(sc_object)$Symbol <- rownames(sc_object)
        rowData(sc_object)$GENEID <- ensembl_dict[rownames(sc_object)]
    } else {
        genes <- rownames(sc_object)
        gene_map <- select(anno_db, genes, c("SYMBOL", "GENENAME"), "GENEID") %>% 
            as_tibble() %>% 
            dplyr::filter(!is.na(SYMBOL) & SYMBOL != "")
        gene_list <- gene_map %>% 
            pull(GENEID) %>% 
            unique()
        sc_object <- sc_object[gene_list]
        ensembl_dict <- gene_map$SYMBOL
        names(ensembl_dict) <- gene_map$GENEID
        rowData(sc_object)$GENEID <- rownames(sc_object)
        rowData(sc_object)$Symbol <- ensembl_dict[rownames(sc_object)]
    }
    rownames(sc_object) <- rowData(sc_object)$GENEID
    rowData(sc_object)$ID <- rowData(sc_object)$GENEID
}



## Plot Barcode rank plots
set.seed(12357)
barcode_rank_table <- DropletUtils::barcodeRanks(counts(sc_object))
knee <- metadata(barcode_rank_table)$knee
inflection <- metadata(barcode_rank_table)$inflection
barcode_rank_plot <- barcode_rank_table %>% 
    as_tibble() %>% 
    ggplot2::ggplot(aes(x = rank, y = total)) + 
    ggplot2::geom_point() +
    ggplot2::geom_line(aes(x = rank, y = fitted), color = "red") +
    ggplot2::ylab("Total") +
    ggplot2::xlab("Rank") + 
    ggplot2::geom_hline(yintercept = knee, color = "green", linetype = "dotted") + 
    ggplot2::geom_hline(yintercept = inflection, color = "blue", linetype = "dotted") + 
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10() + 
    ggplot2::theme_classic(base_size = 28) + 
    ggplot2::theme(legend.position = c(1, 1)) 
barcode_rank_file <- paste(parser$sample_name, "barcode_rank.pdf", sep = "_")
ggplot2::ggsave(barcode_rank_file, plot = barcode_rank_plot, height = 7,
    width = 7, units = "in", device = "pdf") 

# Per cell QC metrics
all_locations <- MatrixGenerics::rowRanges(sc_object)
mito_locations <- rownames(sc_object) %>% str_detect("MT|mt")
names(mito_locations) <- rownames(sc_object)
qc_table <- scuttle::perCellQCMetrics(sc_object, subsets=list(Mito=mito_locations))
sc_object <- scuttle::addPerCellQC(sc_object, subsets = list(Mito = mito_locations))


#Discard reasons
qc_reasons <- scuttle::perCellQCFilters(qc_table, sub.fields = c("subsets_Mito_percent"))
discard_reasons <- qc_reasons %>% 
    as_tibble() %>% 
    summarise(low_lib_size = sum(low_lib_size),
        low_n_features = sum(low_n_features),
        high_subsets_Mito_percent = sum(high_subsets_Mito_percent),
        discard = sum(discard))
# Outlier detection
qc_stats <- cbind(log10(qc_table$sum), log10(qc_table$detected),
    qc_table$subsets_Mito_percent, qc_table$altexps_ERCC_percent)
outlying <- adjOutlyingness(fullRank(qc_stats), only.outlyingness = TRUE)
# Record outlier and discard cell information
discard_reasons <- qc_reasons %>% 
    as_tibble() %>% 
    summarise(low_lib_size = sum(low_lib_size),
        low_n_features = sum(low_n_features),
        high_subsets_Mito_percent = sum(high_subsets_Mito_percent),
        discard = sum(discard), outlier = sum(isOutlier(outlying, type = "higher")))
qc_table$outlyingness <- isOutlier(outlying, type = "higher")
qc_table$discard <- qc_reasons$discard
sc_object$outlier <- isOutlier(outlying, type = "higher")
sc_object$discard <- qc_reasons$discard
# Plot outlier counts
discard_plot <- discard_reasons %>% 
    mutate(sample_name = parser$sample_name) %>% 
    pivot_longer(-sample_name, names_to = "Discard reason", values_to = "Number of cells") %>% 
    ggplot(aes(y= `Discard reason`, x = `Number of cells`)) +
    geom_bar(stat = "identity") +
    scale_y_discrete(labels = c("low_lib_size" = "Low library size", 
        "low_n_features" = "Low feature counts",
        "high_subsets_Mito_percent" = "High mitochondrial content",
        "discard" = "Total discard count", "outlier" = "Outlier cells")) +
    theme_classic(base_size = 28) 
discard_output <- paste(parser$sample_name, "discard.pdf", sep = "_")
ggsave(discard_output, plot = discard_plot, width = 10, height = 4, device = "pdf", units = "in")
# Plot QC characteristics
total_count <- qc_table %>% 
    as_tibble() %>% 
    ggplot(aes(x = discard, y = sum)) + 
    geom_boxplot() + 
    theme_classic(base_size = 28) +
    xlab("Discard") + 
    ylab("Total library size")
total_features <- qc_table %>%
    as_tibble() %>% 
    ggplot(aes(x = discard, y = detected)) + 
    geom_boxplot() + 
    theme_classic(base_size = 28) +
    xlab("Discard") + 
    ylab("Total feature counts")
total_mito <- qc_table %>% 
    as_tibble() %>% 
    ggplot(aes(x = discard, y = subsets_Mito_percent)) + 
    geom_boxplot() + 
    theme_classic(base_size = 28) +
    xlab("Discard") + 
    ylab("Total percent\nmitochondrial genes")
cell_metrics_path <- paste(parser$sample_name, "cell_metrics.pdf", sep = "_")
cell_metrics_plot <- total_count + total_features + total_mito 
ggsave(cell_metrics_path, plot = cell_metrics_plot, width = 14, height = 5, device = "pdf", units = "in")

#Identify doublets from single cell data
cluster_object <- scran::quickCluster(sc_object, d = 10)
sc_object <- scran::computeSumFactors(sc_object, cluster = cluster_object, 
    min.mean = 0.1)
sc_object <- scuttle::logNormCounts(sc_object)
top_variation <- scran::getTopHVGs(sc_object, n = 2000)
sc_object <- scran::fixedPCA(sc_object, subset.row = top_variation)
set.seed(12357)
sc_object <- scater::runTSNE(sc_object, dimred = "PCA")
sc_clusters <- clusterCells(sc_object, use.dimred="PCA", full=TRUE)
colLabels(sc_object) <- sc_clusters$clusters
if (length(unique(sc_clusters$clusters)) > 2) { 
    sc_doublets <- findDoubletClusters(sc_object)
    chosen_doubles <-  rownames(sc_doublets)[isOutlier(sc_doublets$num.de, type="lower", log=TRUE)]
    if (!is.na(max(chosen_doubles))) {
        sc_markers <- scran::findMarkers(sc_object, direction = "up")
        combineTables <- function(table) {
            table_rows <- rownames(table)
            table <- table %>% 
                as_tibble() %>% 
                mutate(genes = table_rows)
            return(table)}
        if( length(sc_markers[chosen_doubles]) == 0) {
            print("No outlier gene expression detected between doublets")
        } else {
            doublet_markers <- sc_markers[chosen_doubles] %>% 
                lapply(combineTables) %>% 
                bind_rows()
            chosen_markers <- doublet_markers %>% 
                dplyr::filter(Top <= 10) %>% 
                pull(genes) %>% 
                unique()
            doublet_heatmap_path <- paste(parser$sample_name, "doublet_heatmap.pdf", sep = "_")
            pdf(doublet_heatmap_path, width = 10, height = 10)
            doublet_heatmap <- scater::plotHeatmap(sc_object, order_columns_by = "label", 
                features = chosen_markers, center = TRUE, symmetric = TRUE, 
                zlim = c(-5,5))
            doublet_heatmap
            dev.off() 
    }} 
    # Identify doublets from RNA-seq data
    sc_object <- scDblFinder(sc_object, clusters=colLabels(sc_object))
    print(table(sc_object$scDblFinder.class))
} else {
    print("Single cluster detected")
}

if ("Barcode" %in% colnames(colData(sc_object))) {
    colnames(sc_object) <- sc_object$Barcode
} 


seurat_object <- as.Seurat(sc_object)
# Set feature metadata, AKA rowData. Super intuitive, right?
try(seurat_object[["RNA"]][[]] <- as.data.frame(rowData(sc_object)), silent = TRUE)
Idents(seurat_object) <- parser$sample_name
seurat_object$ident <- parser$sample_name
seurat_object <- seurat_object %>% 
    subset(subset = outlier == FALSE)
seurat_object_path <- paste(parser$sample_name, ".h5Seurat", sep = "")
SeuratDisk::SaveH5Seurat(seurat_object, seurat_object_path, overwrite=TRUE)
