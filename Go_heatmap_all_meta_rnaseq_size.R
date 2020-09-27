library(pheatmap)
library(dplyr)
library(stringr)
library(tidyverse)
library(writexl)

# For selected Go_terms with multiple sample types
common_go_terms <- readxl::read_xlsx(
  "./Selected_terms_with_all_information_ma.xlsx", sheet = 3) %>%
  as.data.frame() %>% 
  filter(row_number() < 14)

meta_tbl <- common_go_terms %>% 
  separate_rows(Sample_Type, Log2_Fold_Enrichment, sep = ", ")

unique_types <- unique(meta_tbl$Sample_Type)
unique_go_terms <- unique(meta_tbl$Go_terms)

tbl_for_heatmap <- expand_grid(Go_terms = unique_go_terms, 
                               Sample_Type = unique_types, 
                               Log2_Fold_Enrichment = 0) %>% as.data.frame()

tbl_for_heatmap <- tbl_for_heatmap %>% reshape(timevar = "Sample_Type",
          idvar = "Go_terms",
          direction = "wide") 

colnames(tbl_for_heatmap) <- gsub("Log2_Fold_Enrichment.", "", colnames(tbl_for_heatmap))
rownames(tbl_for_heatmap) <- tbl_for_heatmap$Go_terms ; tbl_for_heatmap$Go_terms = NULL

new_col_order <- c("CX_E17.5_cluster2", "CB_E17.5_cluster1", "CB_E17.5_cluster2", "CB_E17.5_cluster3",
                   "HIP_E17.5_cluster1", "HIP_E17.5_cluster2", "HIP_E17.5_cluster3", "CX_P7_cluster3", 
                   "CB_P7_cluster1", "CB_P7_cluster3", "CB_P7_cluster4", "HIP_P7_cluster1", "HIP_P7_cluster3",
                    "CX_P35_cluster3", "CX_P35_cluster4", "CB_P35_cluster1", "CB_P35_cluster2",
                   "CB_P35_cluster3", "HIP_P35_cluster1", "HIP_P35_cluster2", "HIP_P35_cluster3")

new_col_order_by_region <- c("CX_E17.5_cluster2", "CX_P7_cluster3",  "CX_P35_cluster3", 
                             "CX_P35_cluster4","CB_E17.5_cluster1", "CB_E17.5_cluster2", "CB_E17.5_cluster3",
                             "CB_P7_cluster1", "CB_P7_cluster3", "CB_P7_cluster4", "CB_P35_cluster1", 
                             "CB_P35_cluster2","CB_P35_cluster3", "HIP_E17.5_cluster1", "HIP_E17.5_cluster2", 
                             "HIP_E17.5_cluster3",  "HIP_P7_cluster1", "HIP_P7_cluster3", "HIP_P35_cluster1", 
                             "HIP_P35_cluster2", "HIP_P35_cluster3")
tbl_for_heatmap <- tbl_for_heatmap[,new_col_order_by_region]

for (i in seq(length(meta_tbl$Go_terms))){
  tbl_for_heatmap[meta_tbl$Go_terms[i],meta_tbl$Sample_Type[i]] = 
    as.numeric(meta_tbl$Log2_Fold_Enrichment[i])
}

annotation <- data.frame(Sample = colnames(tbl_for_heatmap),row.names = colnames(tbl_for_heatmap))

pdf("./final/Redo_Go_heatmap_selected_new_figure_size_by_regions.pdf", width = 22, height = 16)
pall <- pheatmap((tbl_for_heatmap),  
                 color = colorRampPalette(c("dark blue","white","yellow","red"))(41), 
                 cluster_rows = F, cluster_cols = F, 
                 show_rownames = T, show_colnames = T, 
                 fontsize = 20, fontsize_row = 35, fontsize_col = 35, 
                 legend = T,
                 legend_labels = c("Fold Enrichment"),
                 border_color = "NA", 
                 breaks = NA, 
                 cellwidth = 40, cellheight = 40)
dev.off()

# For single sample type
final_select <- readxl::read_xlsx(
  "../Meta-analysis_selected_common_Go_terms_with_fold_enrichment_MAv2.xlsx", sheet = 3) %>%
  pull(Go_terms)

single_go_terms <- readxl::read_xlsx(
  "../Meta-analysis_selected_common_Go_terms_with_fold_enrichment_MAv2.xlsx", sheet = 4) %>%
  as.data.frame() %>%
  filter(Go_terms %in% final_select)

unique_types_for_single <- unique(single_go_terms$Sample_Type)
unique_go_terms_for_single <- unique(single_go_terms$Go_terms)

tbl_for_heatmap_single <- expand_grid(Go_terms = unique_go_terms_for_single, 
                                      Sample_Type = unique_types_for_single, 
                                      Log2_Fold_Enrichment = 0) %>% 
                          as.data.frame() %>%
                          reshape(timevar = "Sample_Type",
                                  idvar = "Go_terms",
                                  direction = "wide")

colnames(tbl_for_heatmap_single) <- gsub("Log2_Fold_Enrichment.", "", colnames(tbl_for_heatmap_single))
rownames(tbl_for_heatmap_single) <- tbl_for_heatmap_single$Go_terms
tbl_for_heatmap_single$Go_terms = NULL

new_col_order_single <- c(  
  "CX_E17.5_cluster2", "CB_E17.5_cluster2","HIP_E17.5_cluster1", "HIP_E17.5_cluster3", "CB_P7_cluster1",
  "CB_P7_cluster3", "HIP_P7_cluster3", "CX_P35_cluster2", "CX_P35_cluster3", "CX_P35_cluster4", 
  "CB_P35_cluster1", "HIP_P35_cluster3"
                     )

new_col_order_single_by_region <- c(  
  "CX_E17.5_cluster2", "CX_P35_cluster2", "CX_P35_cluster3", "CX_P35_cluster4", 
  "CB_E17.5_cluster2", "CB_P7_cluster1","CB_P7_cluster3", "CB_P35_cluster1",
  "HIP_E17.5_cluster1", "HIP_E17.5_cluster3",  "HIP_P7_cluster3", "HIP_P35_cluster3"
)

tbl_for_heatmap_single <- tbl_for_heatmap_single[final_select,new_col_order_single_by_region]

for (i in seq(length(single_go_terms$Go_terms))){
  tbl_for_heatmap_single[single_go_terms$Go_terms[i],single_go_terms$Sample_Type[i]] = 
    as.numeric(single_go_terms$Log2_Fold_Enrichment[i])
}

annotation_single <- data.frame(Sample = colnames(tbl_for_heatmap_single),
                                row.names = colnames(tbl_for_heatmap_single))

pdf("./final/Go_heatmap_single_final_without_clustering_by_region_shorter_font_35.pdf", width = 18, height = 15)
pall <- pheatmap((tbl_for_heatmap_single),  
                 color = colorRampPalette(c("dark blue","white","yellow","red"))(41), 
                 cluster_rows = F, cluster_cols = F, 
                 show_rownames = T, show_colnames = T, 
                 fontsize = 20, fontsize_row = 35, fontsize_col = 35, 
                 border_color = "NA", 
                 breaks = NA, 
                 cellwidth = 30, cellheight = 30, treeheight_col = 0)
dev.off()








