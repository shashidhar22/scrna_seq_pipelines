library(tidyverse)     
library(lubridate)
library(readxl)
library(patchwork)
library(optparse)
library(scater)
library(DropletUtils)
library(LymphoSeq2)
# Accept arguments
option_list <- list(
    make_option(c("-m", "--meta_path"), 
                type="character",  
                help="Path to metadata file",
                dest="mpath")
)
parser <- parse_args(OptionParser(option_list=option_list), 
                     print_help_and_exit = TRUE)
# Get all the paths for the count matrices and the VDJ rearrangement files
cpaths <- Sys.glob("*/outs/count/filtered_feature_bc_matrix")
vpaths <- Sys.glob("*/outs/vdj_t/airr_rearrangement.tsv")
# Get the study table from the VDJ data for each sample
getVDJ <- function(vdj_path) {
  vdj <- LymphoSeq2::readImmunoSeq(vdj_path) 
  sample <- str_split(vdj_path, "/")[[1]][1]
  vdj <- vdj %>% 
         mutate(repertoire_id = sample)
  return(vdj)
}
# Extract the quality information for the count data include gene and read
# counts and percent mitochondrial content
getQualities <- function(cds_path) {
  sce <- DropletUtils::read10xCounts(cds_path)
  sample <- str_split(cds_path, "/")[[1]][1]
  is_mito <- grepl(rowData(sce)$Symbol, 
                   pattern= "^MT-")
  qc_df <- scater::perCellQCMetrics(sce, 
                                    subsets=list(mitochondrial= is_mito))
  discard <- scater::quickPerCellQC(qc_df, 
                                    percent_subsets=c("subsets_mitochondrial_percent"), 
                                    nmad=3) 
  qc_tibble <- qc_df %>% 
               tibble::as_tibble() %>% 
               tibble::rowid_to_column()
  discard_tibble <- discard %>% 
                    tibble::as_tibble() %>% 
                    tibble::rowid_to_column() %>% 
                    dplyr::select(rowid, discard)
  qc_tibble <- dplyr::left_join(qc_tibble, discard_tibble, by = "rowid") %>% 
               dplyr::mutate(sample = sample)
}
# Generate knee plots for sample to understand the data threshold for Seurat
# data loading
getKneePlot <- function(cds_path) {
  sce <- DropletUtils::read10xCounts(cds_path)
  sample <- str_split(cds_path, "/")[[1]][1]
  sce_matrix <- assay(sce)
  barcode_ranks <- barcodeRanks(sce_matrix)
  barcode_table <- barcode_ranks %>% 
                   as_tibble()
  barcode_points <-  tibble(type = c("inflection", "knee"),
                            value = c(barcode_ranks@metadata$inflection, 
                                      barcode_ranks@metadata$knee))
  knee_plot <- ggplot(barcode_table, aes(x = rank, y = total)) +
               geom_point() +
               geom_hline(aes(yintercept = value,
                              colour = type), 
                          barcode_points,
                          linetype = 2) +
               scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                           labels = trans_format("log10", math_format(10^.x))) +
               scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                           labels = trans_format("log10", math_format(10^.x))) +
               labs(x = "Log Rank",
                    y = "Log counts") +
               theme_classic(base_size = 26)
   out_plot <- paste(sample, "knee", "plot.pdf", sep="_")
   ggsave(out_plot, knee_plot, 
          width = 16, height = 9, 
          units = "in", device = "pdf")
}
# Extract the quality information for also samples and combine across all
# samples and merge with metadata file
qtable <- cpaths %>% 
          purrr::map(getQualities) %>%
          bind_rows()
meta_gex <- read_csv(parser$mpath) %>% filter(str_detect(locus, "GEX"))
qtable <- qtable %>% left_join(meta_gex, by =  c("sample" = "sampleName"))
qmtable <- qtable %>% 
           group_by(sample) %>% 
           summarize(mean_sum = mean(sum),
                     mean_gene = mean(detected),
                     mean_mito = mean(subsets_mitochondrial_percent),
                     tissueType = unique(tissueType),
                     Batch = unique(Batch),
                     patientID = unique(patientID)) %>%
           ungroup()
# Plot percent mitochondrial content distribution
mito_plot <- ggplot(qtable, aes(x=sample, y=subsets_mitochondrial_percent, color=tissueType)) + 
             geom_boxplot() + 
             theme_classic(base_size = 16) +
             theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
             facet_wrap(~Batch, scale="free") +
             labs(x = "Sample Name", y = "Percent Mitochondrial Content", color = "Tissue type") + 
             ggtitle("Per cell mitochondrial content distribution")
# Plot read count distribtuion
read_plot <- ggplot(qtable, aes(x=sample, y=sum, color=tissueType)) + 
             geom_boxplot() + 
             theme_classic(base_size = 16) +
             theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
             facet_wrap(~Batch, scale="free") +
             labs(x = "Sample Name", y = "Total Read Count", color = "Tissue type") + 
             ggtitle("Per cell read count distribution")
# Plot gene count distribution             
gene_plot <- ggplot(qtable, aes(x=sample, y=detected, color=tissueType)) + 
             geom_boxplot() + 
             theme_classic(base_size = 16) +
             theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
             facet_wrap(~Batch, scale="free") +
             labs(x = "Sample Name", y = "Genes Detected", color = "Tissue type") + 
             ggtitle("Genes per cell distribution")
# Plot read count vs gene counts
read_gene <- ggplot(qmtable, aes(x=mean_sum, y = mean_gene, shape =tissueType, color = patientID)) +
             geom_point() +
             theme_classic(base_size = 16) +
             facet_wrap(~Batch, scale="free") +
             xlim(0, 10000) + 
             ylim(0, 2500) +
             labs(x = "Total read count", y = "Genes detected", shape = "Tissue type", color = "Patient ID")
# Plot knee plot
read_matrix <- cpaths %>% 
               purrr::map(getKneePlot)
# Save plots
ggsave("Gene_counts.pdf", gene_plot, width = 16, height = 9, units = "in", device = "pdf")
ggsave("Read_counts.pdf", read_plot, width = 16, height = 9, units = "in", device = "pdf")
ggsave("Mito_content.pdf", mito_plot, width = 16, height = 9, units = "in", device = "pdf")
ggsave("Read_gene.pdf", read_gene, width = 18, height = 6, units = "in", device = "pdf")
# Combine VDJ table from all samples
vtable <- vpaths %>% 
          purrr::map(getVDJ) %>%
          bind_rows()
# Separate the TRA and TRB sequences
tra_table <- vtable %>%
             filter(str_detect(v_call, "TRA"))

trb_table <- vtable %>%
             filter(str_detect(v_call, "TRB"))

tra_itable <- clonality(tra_table) %>%
              left_join(meta_gex, by = c("repertoire_id" = "sampleName"))
trb_itable <- clonality(trb_table) %>%
              left_join(meta_gex, by = c("repertoire_id" = "sampleName"))
# Plot Gini and counts for TRA table
tra_plot <- ggplot(tra_itable, aes(x=unique_productive_sequences, y = gini_coefficient, shape = tissueType, color = patientID)) +
             geom_point(size = 4) +
             theme_classic(base_size = 16) +
             facet_wrap(~Batch, scale="free") +
             xlim(0, 10000) + 
             ylim(0, 1) +
             labs(x = "Unique productive sequences", y = "Gini coefficient", shape = "Tissue type", color = "Patient ID") +
             ggtitle("TRA")
# Plot Gini and counts for TRB table
trb_plot <- ggplot(trb_itable, aes(x=unique_productive_sequences, y = gini_coefficient, shape = tissueType, color = patientID)) +
             geom_point(size = 4) +
             theme_classic(base_size = 16) +
             facet_wrap(~Batch, scale="free") +
             xlim(0, 10000) + 
             ylim(0, 1) +
             labs(x = "Unique productive sequences", y = "Gini coefficient", shape = "Tissue type", color = "Patient ID") +
             ggtitle("TRB")
# Combine and save plots
vdj_plot <- tra_plot + trb_plot + plot_layout(guides = 'collect')
ggsave("VDJ_plot.pdf", vdj_plot, width = 24, height = 6, units = "in", device = "pdf")
# Plot top ten sequences for TRA sequences
top_plot_tra <- tra_table %>%
                producitveSeq() %>% 
                topSeqsPlot(top = 10) + 
                ggtitle("TRA")
# Plot top ten sequences for TRB sequences
top_plot_trb <- trb_table %>%
                producitveSeq() %>% 
                topSeqsPlot(top = 10) +
                ggtitle("TRB")
# Combine and save plots
top_plot <- top_plot_tra + top_plot_trb
ggsave("VDJ_top_plot.pdf", vdj_plot, width = 12, height = 6, units = "in", device = "pdf")
