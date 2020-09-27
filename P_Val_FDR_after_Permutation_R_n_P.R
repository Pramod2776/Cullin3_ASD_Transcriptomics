setwd("~/Desktop/2019_Summer/Iakoucheva_Lab/R_Codes/Cul3_Permutation/CUL3_MICE_Region_Period/DE_NULL_DISTRIBUTIONS")
source("plotting.R")
library(tidyverse)
library(ggplot2)

# Data Processing & Null Distribution Generation

deg_region_file <- list.files(pattern = "deg_PERMUTATION_[0-9]*_Region", full.names = T)
deg_period_file <- list.files(pattern = "deg_PERMUTATION_[0-9]*_Period", full.names = T)
dep_region_file <- list.files(pattern = "dep_PERMUTATION_[0-9]*_Region", full.names = T)
dep_period_file <- list.files(pattern = "dep_PERMUTATION_[0-9]*_Period", full.names = T)


deg_region <- as.character(unique(read.csv(deg_region_file[[1]])$FILE))
deg_period <- as.character(unique(read.csv(deg_period_file[[1]])$FILE))
dep_region <- as.character(unique(read.csv(dep_region_file[[1]])$FILE))
dep_period <- as.character(unique(read.csv(dep_period_file[[1]])$FILE))

lists_names <- as.character(unique(read.csv(deg_region_file[[1]])$GENELIST))

# deg_regionwise
input_matrix_deg_region <- expand_grid(deg_region, lists_names)

null_dist_list_deg_region <- lapply(seq(length(input_matrix_deg_region$deg_region)), function(x){
  temp_result <- lapply(deg_region_file, function(fn){
    temp_sheet <- read.csv(fn) %>%
      filter(FILE == as.character(input_matrix_deg_region[x,1])) %>%
      filter(GENELIST == as.character(input_matrix_deg_region[x,2]))
  }) %>%
    bind_rows() %>%
    arrange(OVERLAP)
}) %>%
  set_names(nm = paste(input_matrix_deg_region$deg_region,
                       input_matrix_deg_region$lists_names, sep = "_"))


# deg_periodwise
input_matrix_deg_period <- expand_grid(deg_period, lists_names)

null_dist_list_deg_period <- lapply(seq(length(input_matrix_deg_period$deg_period)), function(x){
  temp_result <- lapply(deg_period_file, function(fn){
    temp_sheet <- read.csv(fn) %>%
      filter(FILE == as.character(input_matrix_deg_period[x,1])) %>%
      filter(GENELIST == as.character(input_matrix_deg_period[x,2]))
  }) %>%
    bind_rows() %>%
    arrange(OVERLAP)
}) %>%
  set_names(nm = paste(input_matrix_deg_period$deg_period,
                       input_matrix_deg_period$lists_names, sep = "_"))

# dep_regionwise

input_matrix_dep_region <- expand_grid(dep_region, lists_names)

null_dist_list_dep_region <- lapply(seq(length(input_matrix_dep_region$dep_region)), function(x){
  temp_result <- lapply(dep_region_file, function(fn){
    temp_sheet <- read.csv(fn) %>%
      filter(FILE == as.character(input_matrix_dep_region[x,1])) %>%
      filter(GENELIST == as.character(input_matrix_dep_region[x,2]))
  }) %>%
    bind_rows() %>%
    arrange(OVERLAP)
}) %>%
  set_names(nm = paste(input_matrix_dep_region$dep_region,
                       input_matrix_dep_region$lists_names, sep = "_"))


# dep_periodwise

input_matrix_dep_period <- expand_grid(dep_period, lists_names)

null_dist_list_dep_period <- lapply(seq(length(input_matrix_dep_period$dep_period)), function(x){
  temp_result <- lapply(dep_period_file, function(fn){
    temp_sheet <- read.csv(fn) %>%
      filter(FILE == as.character(input_matrix_dep_period[x,1])) %>%
      filter(GENELIST == as.character(input_matrix_dep_period[x,2]))
  }) %>%
    bind_rows() %>%
    arrange(OVERLAP)
}) %>%
  set_names(nm = paste(input_matrix_dep_period$dep_period,
                       input_matrix_dep_period$lists_names, sep = "_"))


# P_val & FDR Calculation
setwd("~/Desktop/2019_Summer/Iakoucheva_Lab/R_Codes/Cul3_Permutation/Permutation_Result")

deg_region_lists <- lapply(setNames(nm = deg_region), function(x){
  list.files("./src/DGE", pattern = x, full.names = T)
})

deg_period_lists <- lapply(setNames(nm = deg_period), function(x){
  list.files("./src/DGE", pattern = paste0("^",x), full.names = T)
})

dep_region_lists <- lapply(setNames(nm = dep_region), function(x){
  list.files("./src/DPE", pattern = x, full.names = T)
})

dep_period_lists <- lapply(setNames(nm = dep_period), function(x){
  list.files("./src/DPE", pattern = paste0("^",x), full.names = T)
})


#deg_regionwise
P_val_deg_region <- lapply(seq(length(input_matrix_deg_region$deg_region)), function(x){
  
  temp_dge_list <- lapply(deg_region_lists[[as.character(input_matrix_deg_region[x,1])]],
                          function(y){
                            temp_dge <- read.csv(y) %>%
                              filter(FDR <= 0.15) %>%
                              pull(Gene_name) %>%
                              toupper()
                          }) %>%
    unlist()
  
  target_list_genes <- readxl::read_xlsx("./Input_files/ASDRelevantGeneListsFromLiterature.xlsx", 
                                         sheet = as.character(input_matrix_deg_region[x,2])) %>%
    pull(gene_symbol)
  
  temp_overlap <- length(intersect(temp_dge_list, target_list_genes))
  
  null_dist <- null_dist_list_deg_region[[paste(as.character(input_matrix_deg_region[x,1]), 
                                     as.character(input_matrix_deg_region[x,2]),
                                     sep = "_")]] %>%
    pull(OVERLAP)
  
  p_val <- (sum(null_dist > temp_overlap)+1)/(length(null_dist)+1)
  
  data.frame("P_val" = c(p_val),
             "Overlap" = c(temp_overlap))
  
}) %>%
  bind_rows()

final_dge_region_permutation_result <- bind_cols(input_matrix_deg_region, P_val_deg_region)

final_dge_region_permutation_result <- lapply(unique(input_matrix_deg_region$deg_region), function(x){
  temp_result <- final_dge_region_permutation_result %>%
    filter(deg_region == x) %>%
    mutate(FDR_BH = p.adjust(P_val, method = "BH"),
           FDR_bonferroni = p.adjust(P_val, method = "bonferroni"),
           Signed_log10_FDR_BH = -log10(FDR_BH),
           Signed_log10_FDR_Bonferroni = -log10(FDR_bonferroni),
           Signed_log10_P_Val = -log10(P_val),
           Data = "DEG",
           STAR_BH = ifelse(FDR_BH <= 0.1, "*", ""),
           STAR_Bonderroni = ifelse(FDR_bonferroni <= 0.1, "*", ""),
           STAR_P_Val = ifelse(P_val <= 0.05, "*", ""),
           label_BH = paste(paste(Overlap, STAR_BH, sep = ""),paste0("(",round(FDR_BH,3),")"), sep = "\n"),
           label_Bonferroni = paste(paste(Overlap, STAR_Bonderroni, 
                                          sep = ""),paste0("(",round(FDR_bonferroni,3),")"), sep = "\n"),
           label_P_Val = paste(paste(Overlap, STAR_P_Val, 
                                     sep = ""),paste0("(",round(P_val,3),")"), sep = "\n")
    )
  
}) %>%
  bind_rows() %>%
  mutate(deg_region = factor(deg_region, levels = c("CX", "CB", "HIP")))

writexl::write_xlsx(final_dge_region_permutation_result, 
                    "./Output_files/final_dge_region_permutation_result_redo.xlsx")

#deg_periodwise
P_val_deg_period <- lapply(seq(length(input_matrix_deg_period$deg_period)), function(x){
  
  temp_dge_list <- lapply(deg_period_lists[[as.character(input_matrix_deg_period[x,1])]],
                          function(y){
                            temp_dge <- read.csv(y) %>%
                              filter(FDR <= 0.15) %>%
                              pull(Gene_name) %>%
                              toupper()
                          }) %>%
    unlist()
  
  target_list_genes <- readxl::read_xlsx("./Input_files/ASDRelevantGeneListsFromLiterature.xlsx", 
                                         sheet = as.character(input_matrix_deg_period[x,2])) %>%
    pull(gene_symbol)
  
  temp_overlap <- length(intersect(temp_dge_list, target_list_genes))
  
  null_dist <- null_dist_list_deg_period[[paste(as.character(input_matrix_deg_period[x,1]), 
                                                as.character(input_matrix_deg_period[x,2]),
                                                sep = "_")]] %>%
    pull(OVERLAP)
  
  p_val <- (sum(null_dist > temp_overlap)+1)/(length(null_dist)+1)
  
  data.frame("P_val" = c(p_val),
             "Overlap" = c(temp_overlap))
  
}) %>%
  bind_rows()

final_dge_period_permutation_result <- bind_cols(input_matrix_deg_period, P_val_deg_period)

final_dge_period_permutation_result <- lapply(unique(input_matrix_deg_period$deg_period), function(x){
  temp_result <- final_dge_period_permutation_result %>%
    filter(deg_period == x) %>%
    mutate(FDR_BH = p.adjust(P_val, method = "BH"),
           FDR_bonferroni = p.adjust(P_val, method = "bonferroni"),
           Signed_log10_FDR_BH = -log10(FDR_BH),
           Signed_log10_FDR_Bonferroni = -log10(FDR_bonferroni),
           Signed_log10_P_Val = -log10(P_val),
           Data = "DEG",
           STAR_BH = ifelse(FDR_BH <= 0.1, "*", ""),
           STAR_Bonderroni = ifelse(FDR_bonferroni <= 0.1, "*", ""),
           STAR_P_Val = ifelse(P_val <= 0.05, "*", ""),
           label_BH = paste(paste(Overlap, STAR_BH, sep = ""),paste0("(",round(FDR_BH,3),")"), sep = "\n"),
           label_Bonferroni = paste(paste(Overlap, STAR_Bonderroni, 
                                          sep = ""),paste0("(",round(FDR_bonferroni,3),")"), sep = "\n"),
           label_P_Val = paste(paste(Overlap, STAR_P_Val, 
                                     sep = ""),paste0("(",round(P_val,3),")"), sep = "\n")
           )
  
}) %>%
  bind_rows() %>%
  mutate(deg_period = factor(deg_period, levels = c("E17.5", "P7", "P35")))

writexl::write_xlsx(final_dge_period_permutation_result, 
                    "./Output_files/final_dge_period_permutation_result_redo.xlsx")

#dep_regionwise
P_val_dep_region <- lapply(seq(length(input_matrix_dep_region$dep_region)), function(x){
  
  temp_dpe_list <- lapply(dep_region_lists[[as.character(input_matrix_dep_region[x,1])]],
                          function(y){
                            temp_dpe <- read.csv(y) %>%
                              filter(adj.P.Val <= 0.1) %>%
                              pull(GeneSymbol) %>%
                              toupper()
                          }) %>%
    unlist()
  
  target_list_genes <- readxl::read_xlsx("./Input_files/ASDRelevantGeneListsFromLiterature.xlsx", 
                                         sheet = as.character(input_matrix_dep_region[x,2])) %>%
    pull(gene_symbol)
  
  temp_overlap <- length(intersect(temp_dpe_list, target_list_genes))
  
  null_dist <- null_dist_list_dep_region[[paste(as.character(input_matrix_dep_region[x,1]), 
                                                as.character(input_matrix_dep_region[x,2]),
                                                sep = "_")]] %>%
    pull(OVERLAP)
  
  p_val <- (sum(null_dist > temp_overlap)+1)/(length(null_dist)+1)
  
  data.frame("P_val" = c(p_val),
             "Overlap" = c(temp_overlap))
  
}) %>%
  bind_rows()

final_dpe_region_permutation_result <- bind_cols(input_matrix_dep_region, P_val_dep_region)

final_dpe_region_permutation_result <- lapply(unique(input_matrix_dep_region$dep_region), function(x){
  temp_result <- final_dpe_region_permutation_result %>%
    filter(dep_region == x) %>%
    mutate(FDR_BH = p.adjust(P_val, method = "BH"),
           FDR_bonferroni = p.adjust(P_val, method = "bonferroni"),
           Signed_log10_FDR_BH = -log10(FDR_BH),
           Signed_log10_FDR_Bonferroni = -log10(FDR_bonferroni),
           Signed_log10_P_Val = -log10(P_val),
           Data = "DEP",
           STAR_BH = ifelse(FDR_BH <= 0.1, "*", ""),
           STAR_Bonderroni = ifelse(FDR_bonferroni <= 0.1, "*", ""),
           STAR_P_Val = ifelse(P_val <= 0.05, "*", ""),
           label_BH = paste(paste(Overlap, STAR_BH, sep = ""),paste0("(",round(FDR_BH,3),")"), sep = "\n"),
           label_Bonferroni = paste(paste(Overlap, STAR_Bonderroni, 
                                          sep = ""),paste0("(",round(FDR_bonferroni,3),")"), sep = "\n"),
           label_P_Val = paste(paste(Overlap, STAR_P_Val, 
                                     sep = ""),paste0("(",round(P_val,3),")"), sep = "\n")
    ) 
  
}) %>%
  bind_rows() %>%
  mutate(dep_region = factor(dep_region, levels = c("CX", "CB")))

writexl::write_xlsx(final_dpe_region_permutation_result, 
                    "./Output_files/final_dpe_region_permutation_result_redo.xlsx")


#dep_periodwise
P_val_dep_period <- lapply(seq(length(input_matrix_dep_period$dep_period)), function(x){
  
  temp_dpe_list <- lapply(dep_period_lists[[as.character(input_matrix_dep_period[x,1])]],
                          function(y){
                            temp_dpe <- read.csv(y) %>%
                              filter(adj.P.Val <= 0.1) %>%
                              pull(GeneSymbol) %>%
                              toupper()
                          }) %>%
    unlist()
  
  target_list_genes <- readxl::read_xlsx("./Input_files/ASDRelevantGeneListsFromLiterature.xlsx", 
                                         sheet = as.character(input_matrix_dep_period[x,2])) %>%
    pull(gene_symbol)
  
  temp_overlap <- length(intersect(temp_dpe_list, target_list_genes))
  
  null_dist <- null_dist_list_dep_period[[paste(as.character(input_matrix_dep_period[x,1]), 
                                                as.character(input_matrix_dep_period[x,2]),
                                                sep = "_")]] %>%
    pull(OVERLAP)
  
  p_val <- (sum(null_dist > temp_overlap)+1)/(length(null_dist)+1)
  
  data.frame("P_val" = c(p_val),
             "Overlap" = c(temp_overlap))
  
}) %>%
  bind_rows()

final_dpe_period_permutation_result <- bind_cols(input_matrix_dep_period, P_val_dep_period)

final_dpe_period_permutation_result <- lapply(unique(input_matrix_dep_period$dep_period), function(x){
  temp_result <- final_dpe_period_permutation_result %>%
    filter(dep_period == x) %>%
    mutate(FDR_BH = p.adjust(P_val, method = "BH"),
           FDR_bonferroni = p.adjust(P_val, method = "bonferroni"),
           Signed_log10_FDR_BH = -log10(FDR_BH),
           Signed_log10_FDR_Bonferroni = -log10(FDR_bonferroni),
           Signed_log10_P_Val = -log10(P_val),
           Data = "DEP",
           STAR_BH = ifelse(FDR_BH <= 0.1, "*", ""),
           STAR_Bonderroni = ifelse(FDR_bonferroni <= 0.1, "*", ""),
           STAR_P_Val = ifelse(P_val <= 0.05, "*", ""),
           label_BH = paste(paste(Overlap, STAR_BH, sep = ""),paste0("(",round(FDR_BH,3),")"), sep = "\n"),
           label_Bonferroni = paste(paste(Overlap, STAR_Bonderroni, 
                                          sep = ""),paste0("(",round(FDR_bonferroni,3),")"), sep = "\n"),
           label_P_Val = paste(paste(Overlap, STAR_P_Val, 
                                     sep = ""),paste0("(",round(P_val,3),")"), sep = "\n")
           ) 
  
}) %>%
  bind_rows() %>%
  mutate(dep_period = factor(dep_period, levels = c("E17.5", "P7", "P35")))

writexl::write_xlsx(final_dpe_period_permutation_result, 
                    "./Output_files/final_dpe_period_permutation_result_redo.xlsx")

# Plotting_BH
enr_all_deg_region_BH <- ggplot() +
  facet_nested(. ~ Data + deg_region, scales = "free") +
  geom_tile(
    data = final_dge_region_permutation_result %>%
      filter(FDR_BH > 0.1),
    mapping = aes(x = "", y = lists_names),
    colour = "black", fill = "white"
  ) +
  geom_tile(
    data = final_dge_region_permutation_result %>%
      filter(FDR_BH <= 0.1),
    mapping = aes(x = "", y = lists_names, fill = Signed_log10_FDR_BH),
    colour = "black"
  ) +
  geom_text(
    data = final_dge_region_permutation_result,
    mapping = aes(x = "", y = lists_names, label = label_BH),
    size = 10
  ) +
  scale_fill_gradient(
    low = "grey90", high = "red"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "lines")
  )


enr_all_deg_period_BH <- ggplot() +
  facet_nested(. ~ Data + deg_period, scales = "free") +
  geom_tile(
    data = final_dge_period_permutation_result %>%
      filter(FDR_BH > 0.1),
    mapping = aes(x = "", y = lists_names),
    colour = "black", fill = "white"
  ) +
  geom_tile(
    data = final_dge_period_permutation_result %>%
      filter(FDR_BH <= 0.1),
    mapping = aes(x = "", y = lists_names, fill = Signed_log10_FDR_BH),
    colour = "black"
  ) +
  geom_text(
    data = final_dge_period_permutation_result,
    mapping = aes(x = "", y = lists_names, label = label_BH),
    size = 10
  ) +
  scale_fill_gradient(
    low = "grey90", high = "red"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "lines")
  )

enr_all_dep_period_BH <- ggplot() +
  facet_nested(. ~ Data + dep_period, scales = "free") +
  geom_tile(
    data = final_dpe_period_permutation_result %>%
      filter(FDR_BH > 0.1),
    mapping = aes(x = "", y = lists_names),
    colour = "black", fill = "white"
  ) +
  geom_tile(
    data = final_dpe_period_permutation_result %>%
      filter(FDR_BH <= 0.1),
    mapping = aes(x = "", y = lists_names, fill = Signed_log10_FDR_BH),
    colour = "black"
  ) +
  geom_text(
    data = final_dpe_period_permutation_result,
    mapping = aes(x = "", y = lists_names, label = label_BH),
    size = 10
  ) +
  scale_fill_gradient(
    low = "grey90", high = "red"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "lines")
  )


enr_all_dep_region_BH <- ggplot() +
  facet_nested(. ~ Data + dep_region, scales = "free") +
  geom_tile(
    data = final_dpe_region_permutation_result %>%
      filter(FDR_BH > 0.1),
    mapping = aes(x = "", y = lists_names),
    colour = "black", fill = "white"
  ) +
  geom_tile(
    data = final_dpe_region_permutation_result %>%
      filter(FDR_BH <= 0.1),
    mapping = aes(x = "", y = lists_names, fill = Signed_log10_FDR_BH),
    colour = "black"
  ) +
  geom_text(
    data = final_dpe_region_permutation_result,
    mapping = aes(x = "", y = lists_names, label = label_BH),
    size = 10
  ) +
  scale_fill_gradient(
    low = "grey90", high = "red"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "lines")
  )


#Plotting_Bonferroni
enr_all_deg_region_Bonferroni <- ggplot() +
  facet_nested(. ~ Data + deg_region, scales = "free") +
  geom_tile(
    data = final_dge_region_permutation_result %>%
      filter(FDR_bonferroni > 0.1),
    mapping = aes(x = "", y = lists_names),
    colour = "black", fill = "white"
  ) +
  geom_tile(
    data = final_dge_region_permutation_result %>%
      filter(FDR_bonferroni <= 0.1),
    mapping = aes(x = "", y = lists_names, fill = Signed_log10_FDR_BH),
    colour = "black"
  ) +
  geom_text(
    data = final_dge_region_permutation_result,
    mapping = aes(x = "", y = lists_names, label = label_Bonferroni),
    size = 10
  ) +
  scale_fill_gradient(
    low = "grey90", high = "red"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "lines")
  )


enr_all_deg_period_Bonferroni <- ggplot() +
  facet_nested(. ~ Data + deg_period, scales = "free") +
  geom_tile(
    data = final_dge_period_permutation_result %>%
      filter(FDR_bonferroni > 0.1),
    mapping = aes(x = "", y = lists_names),
    colour = "black", fill = "white"
  ) +
  geom_tile(
    data = final_dge_period_permutation_result %>%
      filter(FDR_bonferroni <= 0.1),
    mapping = aes(x = "", y = lists_names, fill = Signed_log10_FDR_Bonferroni),
    colour = "black"
  ) +
  geom_text(
    data = final_dge_period_permutation_result,
    mapping = aes(x = "", y = lists_names, label = label_Bonferroni),
    size = 10
  ) +
  scale_fill_gradient(
    low = "grey90", high = "red"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "lines")
  )

enr_all_dep_period_Bonferroni <- ggplot() +
  facet_nested(. ~ Data + dep_period, scales = "free") +
  geom_tile(
    data = final_dpe_period_permutation_result %>%
      filter(FDR_bonferroni > 0.1),
    mapping = aes(x = "", y = lists_names),
    colour = "black", fill = "white"
  ) +
  geom_tile(
    data = final_dpe_period_permutation_result %>%
      filter(FDR_bonferroni <= 0.1),
    mapping = aes(x = "", y = lists_names, fill = Signed_log10_FDR_Bonferroni),
    colour = "black"
  ) +
  geom_text(
    data = final_dpe_period_permutation_result,
    mapping = aes(x = "", y = lists_names, label = label_Bonferroni),
    size = 10
  ) +
  scale_fill_gradient(
    low = "grey90", high = "red"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "lines")
  )


enr_all_dep_region_Bonferroni <- ggplot() +
  facet_nested(. ~ Data + dep_region, scales = "free") +
  geom_tile(
    data = final_dpe_region_permutation_result %>%
      filter(FDR_bonferroni > 0.1),
    mapping = aes(x = "", y = lists_names),
    colour = "black", fill = "white"
  ) +
  geom_tile(
    data = final_dpe_region_permutation_result %>%
      filter(FDR_bonferroni <= 0.1),
    mapping = aes(x = "", y = lists_names, fill = Signed_log10_FDR_Bonferroni),
    colour = "black"
  ) +
  geom_text(
    data = final_dpe_region_permutation_result,
    mapping = aes(x = "", y = lists_names, label = label_Bonferroni),
    size = 10
  ) +
  scale_fill_gradient(
    low = "grey90", high = "red"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "lines")
  )


#Plotting P_val
enr_all_deg_region_P_Val <- ggplot() +
  facet_nested(. ~ Data + deg_region, scales = "free") +
  geom_tile(
    data = final_dge_region_permutation_result %>%
      filter(P_val > 0.05),
    mapping = aes(x = "", y = lists_names),
    colour = "black", fill = "white"
  ) +
  geom_tile(
    data = final_dge_region_permutation_result %>%
      filter(P_val <= 0.05),
    mapping = aes(x = "", y = lists_names, fill = Signed_log10_P_Val),
    colour = "black"
  ) +
  geom_text(
    data = final_dge_region_permutation_result,
    mapping = aes(x = "", y = lists_names, label = label_P_Val),
    size = 10
  ) +
  scale_fill_gradient(
    low = "grey90", high = "red"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "lines")
  )


enr_all_deg_period_P_Val <- ggplot() +
  facet_nested(. ~ Data + deg_period, scales = "free") +
  geom_tile(
    data = final_dge_period_permutation_result %>%
      filter(P_val > 0.05),
    mapping = aes(x = "", y = lists_names),
    colour = "black", fill = "white"
  ) +
  geom_tile(
    data = final_dge_period_permutation_result %>%
      filter(P_val <= 0.05),
    mapping = aes(x = "", y = lists_names, fill = Signed_log10_P_Val),
    colour = "black"
  ) +
  geom_text(
    data = final_dge_period_permutation_result,
    mapping = aes(x = "", y = lists_names, label = label_P_Val),
    size = 10
  ) +
  scale_fill_gradient(
    low = "grey90", high = "red"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "lines")
  )

enr_all_dep_period_P_Val <- ggplot() +
  facet_nested(. ~ Data + dep_period, scales = "free") +
  geom_tile(
    data = final_dpe_period_permutation_result %>%
      filter(P_val > 0.05),
    mapping = aes(x = "", y = lists_names),
    colour = "black", fill = "white"
  ) +
  geom_tile(
    data = final_dpe_period_permutation_result %>%
      filter(P_val <= 0.05),
    mapping = aes(x = "", y = lists_names, fill = Signed_log10_P_Val),
    colour = "black"
  ) +
  geom_text(
    data = final_dpe_period_permutation_result,
    mapping = aes(x = "", y = lists_names, label = label_P_Val),
    size = 10
  ) +
  scale_fill_gradient(
    low = "grey90", high = "red"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "lines")
  )


enr_all_dep_region_P_Val <- ggplot() +
  facet_nested(. ~ Data + dep_region, scales = "free") +
  geom_tile(
    data = final_dpe_region_permutation_result %>%
      filter(P_val > 0.05),
    mapping = aes(x = "", y = lists_names),
    colour = "black", fill = "white"
  ) +
  geom_tile(
    data = final_dpe_region_permutation_result %>%
      filter(P_val <= 0.05),
    mapping = aes(x = "", y = lists_names, fill = Signed_log10_P_Val),
    colour = "black"
  ) +
  geom_text(
    data = final_dpe_region_permutation_result,
    mapping = aes(x = "", y = lists_names, label = label_P_Val),
    size = 10
  ) +
  scale_fill_gradient(
    low = "grey90", high = "red"
  ) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    panel.spacing = unit(0, "lines")
  )


pdf("./Figures/DE_ListEnrichment_Permutation_FDR_BH_redo.pdf", width = 16, height = 9)
invisible(
  lapply(
    list(
      enr_all_deg_region_BH,
      enr_all_deg_period_BH,
      enr_all_dep_region_BH,
      enr_all_dep_period_BH
    ),
    function(p) {
      p <- p + 
        theme(
          text = element_text(size = 25)
        )
      print(p)
    }
  )
)
dev.off()

pdf("./Figures/DE_ListEnrichment_Permutation_FDR_Bonferroni_redo.pdf", width = 16, height = 9)
invisible(
  lapply(
    list(
      enr_all_deg_region_Bonferroni,
      enr_all_deg_period_Bonferroni,
      enr_all_dep_region_Bonferroni,
      enr_all_dep_period_Bonferroni
    ),
    function(p) {
      p <- p + 
        theme(
          text = element_text(size = 25)
        )
      print(p)
    }
  )
)
dev.off()

pdf("./Figures/DE_ListEnrichment_Permutation_P_Val_redo.pdf", width = 16, height = 9)
invisible(
  lapply(
    list(
      enr_all_deg_region_P_Val,
      enr_all_deg_period_P_Val,
      enr_all_dep_region_P_Val,
      enr_all_dep_period_P_Val
    ),
    function(p) {
      p <- p + 
        theme(
          text = element_text(size = 25)
        )
      print(p)
    }
  )
)
dev.off()



