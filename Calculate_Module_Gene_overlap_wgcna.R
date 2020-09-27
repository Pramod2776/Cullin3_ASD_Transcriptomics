library(tidyverse)

# gene_lists from literature
gene_lists <- lapply(
  setNames(nm = readxl::excel_sheets('src/ASDRelevantGeneListsFromLiterature.xlsx')[-13]),
  function(sheet) {
    readxl::read_xlsx(
      'src/ASDRelevantGeneListsFromLiterature.xlsx',
      sheet = sheet
    )[, 1] %>%
      distinct() %>%
      pull(gene_symbol)
  }
)
gene_lists <- gene_lists[!names(gene_lists) %in% c(
  "EichlerDNM_LGD_MIS_ASD_BD_SCZ", 
  "EichlerDNM_LGD_ASD_BD_SCZ",
  "VulnerableASD",
  "ASDSanders65"
)]

# Module_membership

six_p_n_r <- c("E17.5_CX", "E17.5_CB", "P7_CX", "P7_CB", "P35_CX", "P35_CB")
permutation_results <- lapply(setNames(nm = six_p_n_r), function(x){
  
  list.files("./data/WPCNA_NULL_DISTRIBUTIONS", pattern = x, full.names = T)
  
})


final_overlap_results <- lapply(setNames(nm = names(permutation_results)), function(x){
  
  temp_overlap_results <- lapply(permutation_results[[x]], function(y){
    
    temp_sheet <- read.csv(y) %>%
      lapply(., function(r) as.character(r)) %>%
      as.data.frame()
      
    Resample <- as.character(unique(temp_sheet$RESAMPLE))
    Module_number <- as.character(unique(temp_sheet$label))
    input_matrix <- expand_grid(Module_number, Resample)
    
    sec_overlap_result <- lapply(seq(length(input_matrix$Module_number)), function(z){
      
      print(z)
      temp_genes <- temp_sheet %>%
        filter(RESAMPLE == as.numeric(input_matrix[z,2])) %>%
        filter(label == as.character(input_matrix[z,1])) %>%
        pull(null_genes) %>%
        as.character()
      
      overlaps <- lapply(names(gene_lists), function(l){
        temp_overlap <- length(intersect(temp_genes, gene_lists[[l]]))
        
        temp_row <- data.frame("Module_number" = c(as.character(input_matrix[z,1])), 
                   "Resample" = c(as.character(input_matrix[z,2])),
                   "GENELIST" = c(l),
                   "Overlap" = c(temp_overlap))
      }) %>%
        bind_rows()
      
    }) %>%
      bind_rows()
    
  }) %>%
    bind_rows()
  
  write_csv(temp_overlap_results, 
            paste0("./Permutation_after_WPCNA/Module_gene_overlap/",
                  "Null_dist_Permutation_WGCNA_",
                  x,
                  ".csv"))
  
})


















