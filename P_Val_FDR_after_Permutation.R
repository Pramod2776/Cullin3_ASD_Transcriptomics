setwd("~/Desktop/2019_Summer/Iakoucheva_Lab/R_Codes/Cul3_Permutation/CUL3_MICE/data/DE_NULL_DISTRIBUTIONS")
library(tidyverse)

# Data Processing & Null Distribution Generation

deg_files <- list.files(pattern = "deg", full.names = T)
dep_files <- list.files(pattern = "dep", full.names = T)

test <- read.csv("./null_distribution_deg_PERMUTATION_1.csv")
nine_dge_name <- paste(str_split(
  str_split(unique(test$FILE), pattern = "/", simplify = T)[,2], 
  pattern = "_", 
  simplify = T)[,1],str_split(
    str_split(unique(test$FILE), pattern = "/", simplify = T)[,2], 
    pattern = "_", 
    simplify = T)[,2], 
  sep = "_")

test_1 <- read.csv("./null_distribution_dep_PERMUTATION_9.csv")
six_dpe_name <- paste(str_split(
  str_split(unique(test_1$FILE), pattern = "/", simplify = T)[,2], 
  pattern = "_", 
  simplify = T)[,1],str_split(
    str_split(unique(test_1$FILE), pattern = "/", simplify = T)[,2], 
    pattern = "_", 
    simplify = T)[,2], 
  sep = "_")


lists_names <- as.character(unique(test$GENELIST))

input_matrix <- expand_grid(nine_dge_name, lists_names)

null_dist_list <- lapply(seq(length(input_matrix$nine_dge_name)), function(x){
  temp_result <- lapply(deg_files, function(fn){
    temp_sheet <- read.csv(fn) %>%
      filter(grepl(as.character(input_matrix[x,1]), FILE)) %>%
      filter(GENELIST == as.character(input_matrix[x,2]))
  }) %>%
    bind_rows() %>%
    arrange(OVERLAP)
}) %>%
  set_names(nm = paste(input_matrix$nine_dge_name,input_matrix$lists_names, sep = "_"))



lists_names_dep <- as.character(unique(test_1$GENELIST))

input_matrix_dep <- expand_grid(six_dpe_name, lists_names_dep)

null_dist_list_dep <- lapply(seq(length(input_matrix_dep$six_dpe_name)), function(x){
  temp_result <- lapply(dep_files, function(fn){
    temp_sheet <- read.csv(fn) %>%
      filter(grepl(as.character(input_matrix_dep[x,1]), FILE)) %>%
      filter(GENELIST == as.character(input_matrix_dep[x,2]))
  }) %>%
    bind_rows() %>%
    arrange(OVERLAP)
}) %>%
  set_names(nm = paste(input_matrix_dep$six_dpe_name,input_matrix_dep$lists_names_dep, sep = "_"))


# P_val & FDR Calculation
setwd("~/Desktop/2019_Summer/Iakoucheva_Lab/R_Codes/Cul3_Permutation/Permutation_Result")

dge_lists <- list.files("./src/DGE", full.names = T)
dpe_lists <- list.files("./src/DPE", full.names = T)

P_val_deg <- lapply(seq(length(input_matrix$nine_dge_name)), function(x){
  
  temp_dge <- read.csv(dge_lists[grepl(as.character(input_matrix[x,1]), dge_lists)]) %>%
    filter(FDR <= 0.1) %>%
    pull(Gene_name) %>%
    toupper()
    
  target_list_genes <- readxl::read_xlsx("./Input_files/ASDRelevantGeneListsFromLiterature.xlsx", 
                                   sheet = as.character(input_matrix[x,2])) %>%
    pull(gene_symbol)
    
  temp_overlap <- length(intersect(temp_dge, target_list_genes))
  
  null_dist <- null_dist_list[[paste(as.character(input_matrix[x,1]), 
                                     as.character(input_matrix[x,2]),
                                     sep = "_")]] %>%
    pull(OVERLAP)
  
  p_val <- (sum(null_dist > temp_overlap)+1)/(length(null_dist)+1)
  
  data.frame("P_val" = c(p_val),
             "Overlap" = c(temp_overlap))
  
}) %>%
  bind_rows()

final_dge_permutation_result <- bind_cols(input_matrix, P_val_deg)


P_val_dep <- lapply(seq(length(input_matrix_dep$six_dpe_name)), function(x){
  
  temp_dpe <- read.csv(dpe_lists[grepl(as.character(input_matrix_dep[x,1]), dpe_lists)]) %>%
    filter(adj.P.Val <= 0.15) %>%
    pull(GeneSymbol) %>%
    toupper()
  
  target_list_genes <- readxl::read_xlsx("./Input_files/ASDRelevantGeneListsFromLiterature.xlsx", 
                                         sheet = as.character(input_matrix_dep[x,2])) %>%
    pull(gene_symbol)
  
  temp_overlap <- length(intersect(temp_dpe, target_list_genes))
  
  null_dist <- null_dist_list_dep[[paste(as.character(input_matrix_dep[x,1]), 
                                     as.character(input_matrix_dep[x,2]),
                                     sep = "_")]] %>%
    pull(OVERLAP)
  
  p_val <- (sum(null_dist > temp_overlap)+1)/(length(null_dist)+1)
  
  data.frame("P_val" = c(p_val),
             "Overlap" = c(temp_overlap))
  
}) %>%
  bind_rows()


final_dpe_permutation_result <- bind_cols(input_matrix_dep, P_val_dep)


final_dge_permutation_result <- lapply(unique(input_matrix$nine_dge_name), function(x){
  temp_result <- final_dge_permutation_result %>%
    filter(nine_dge_name == x) %>%
    mutate(FDR_BH = p.adjust(P_val, method = "BH"),
           FDR_bonferroni = p.adjust(P_val, method = "bonferroni"))
  
}) %>%
  bind_rows()

final_dpe_permutation_result <- lapply(unique(input_matrix_dep$six_dpe_name), function(x){
  temp_result <- final_dpe_permutation_result %>%
    filter(six_dpe_name == x) %>%
    mutate(FDR_BH = p.adjust(P_val, method = "BH"),
           FDR_bonferroni = p.adjust(P_val, method = "bonferroni"))
  
}) %>%
  bind_rows()

final_dge_permutation_result <- final_dge_permutation_result %>%
  rename("nine_dge_name" = "DGE_List") %>%
  mutate(Period = str_split(DGE_List, "_", simplify = T)[,1],
         Region = str_split(DGE_List, "_", simplify = T)[,2])

final_dpe_permutation_result <- final_dpe_permutation_result %>%
  rename("six_dpe_name" = "DPE_List") %>%
  mutate(Period = str_split(DPE_List, "_", simplify = T)[,1],
         Region = str_split(DPE_List, "_", simplify = T)[,2])

writexl::write_xlsx(final_dge_permutation_result, "./Output_files/CUL3_DGE_Permutation_Result_015.xlsx")
writexl::write_xlsx(final_dpe_permutation_result, "./Output_files/CUL3_DPE_Permutation_Result_015.xlsx")

final_dge_permutation_result <- final_dge_permutation_result %>%
  mutate(Period = str_split(DGE_List, "_", simplify = T)[,1],
         Region = str_split(DGE_List, "_", simplify = T)[,2])


