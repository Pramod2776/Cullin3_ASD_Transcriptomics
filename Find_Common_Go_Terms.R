library(dplyr)
library(stringr)
library(tidyverse)
library(writexl)
library(gProfileR) # GO-Enrichment
library(pheatmap) # For Heatmap


Single_dge_bg_files <- list.files(path = "./GO_Enrich_Results_Single_Fold_Enrich",
                                  full.names = T)


elements_names <- lapply(Single_dge_bg_files, function(x){
  intermediate <- str_split(x, "/")[[1]][3] %>%
    str_split(.,"_")
  result <- paste(intermediate[[1]][[8]],
                  intermediate[[1]][[9]],
                  intermediate[[1]][[10]],
                  sep = "_")
})

extract_term_names <- lapply(Single_dge_bg_files, function(x){
  temp_col <- readxl::read_xlsx(x, sheet = 1)
  if (length(colnames(temp_col)) != 0){
    result <- temp_col %>% 
     pull(term.name)
  } else {
    result <- character()
  }
}) %>%
  set_names(elements_names)

extract_fold_enrich <- lapply(Single_dge_bg_files, function(x){
  temp_col <- readxl::read_xlsx(x, sheet = 1)
  if (length(colnames(temp_col)) != 0){
    result <- temp_col %>% pull(log2_fold_enrichment)
  } else {
    result <- character()
  }
})

selected_go_terms <- readxl::read_xlsx("./selected_common_Go_terms-old.xlsx",sheet = 1) %>%
  pull(Go_terms)

common_Go_terms <- unlist(extract_term_names,use.names = T) %>%
  as.data.frame() %>%
  set_names(.,"Go_terms") %>%
  mutate(Sample_Type = names(unlist(extract_term_names,use.names = T)) %>%
           gsub(pattern = "[0-9]+$", "",.),
         Log2_Fold_Enrichment = unlist(extract_fold_enrich) %>%
           as.numeric()) %>% 
  group_by(Go_terms) %>%
  filter(Go_terms %in% selected_go_terms) #%>%
  # summarise(Sample_Type = paste(Sample_Type, collapse=", "),
  #           Log2_Fold_Enrichment = paste(Log2_Fold_Enrichment, collapse=", ") )

#write out the common_go_term sheet with log_fold_enrichment
writexl::write_xlsx(common_Go_terms, "./common_Go_terms_with_fold_enrichment.xlsx",col_names = T)

selected_color_go_terms <- readxl::read_xlsx("./Meta_Combined_Go-terms_sorted_colored.xlsx",sheet = 1) %>%
  pull(Go_terms)

common_Go_terms <- common_Go_terms %>% filter(Go_terms %in% selected_color_go_terms)

writexl::write_xlsx(common_Go_terms, "./selected_common_Go_terms_with_fold_enrichment.xlsx",col_names = T)


test <- expand.grid("Go_terms" = common_Go_terms$Go_terms %>% unique(), 
                    "Sample_Type" = common_Go_terms$Sample_Type %>% unique()) %>%
  mutate(additional_col = 1)

common_Go_terms <- common_Go_terms %>%
  left_join(test, ., by = c("Go_terms","Sample_Type")) %>%
  select(c("Go_terms","Sample_Type","Log2_Fold_Enrichment"))

common_Go_terms[is.na(common_Go_terms)] <- 1

common_Go_terms <- common_Go_terms %>%
  reshape(.,timevar = "Sample_Type",
          idvar = "Go_terms",
          direction = "wide")

rownames(common_Go_terms) <- common_Go_terms$Go_terms; common_Go_terms$Go_terms=NULL

colnames(common_Go_terms) <- gsub(replacement = "",pattern = "Log2_Fold_Enrichment\\.",colnames(common_Go_terms))

annotation <- data.frame(Sample = colnames(common_Go_terms))
row.names(annotation) <- colnames(common_Go_terms)

common_Go_terms <- common_Go_terms[order(row.names(common_Go_terms)),]

pdf("./Go_heatmap_selected.pdf", width = 40, height = 250)
pall <- pheatmap((common_Go_terms),  
                 color = colorRampPalette(c("dark blue","white","yellow","red"))(41), 
                 cluster_rows = T, cluster_cols = F, 
                 show_rownames = T, show_colnames = T, 
                 fontsize = 15, fontsize_row = 12, fontsize_col = 10, 
                 annotation = annotation,
                 border_color = "NA", 
                 breaks = NA, 
                 cellwidth = 10, cellheight = 10)
dev.off()

writexl::write_xlsx(common_Go_terms, path = "./common_Go_terms.xlsx",col_names = T)  



