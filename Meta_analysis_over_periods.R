library(tidyverse)
library(dplyr)
library(gProfileR) # GO-Enrichment

three_periods <- c("E17.5","P7","P35")

Meta_analysis_results <- lapply(set_names(three_periods), function(x){
  
  three_regions <- c("CB","CX","HIP")
  periods_fns <- list.files(path = "./input_files", pattern = paste0("^",x), full.names = T)
  
  DGEs <- lapply(periods_fns, function(f){
    region <- str_split(f, "_")[[1]][3]
    temp_dge <- readxl::read_xlsx(f, sheet = 1) %>%
      select(Ensembl_ID...1, logFC, PValue) %>% 
      set_names(
        c("Ensembl_Gene_ID", paste0(region,"_LFC"), paste0(region,"_P_Val"))
      )
  }) %>% 
    set_names(three_regions)
  
  FC <- DGEs[["CX"]] %>% 
    full_join(DGEs[["CB"]]) %>%
    full_join(DGEs[["HIP"]])
 
  three_regions <- c("CX","CB","HIP") 
  sign_cols <- lapply(three_regions, function(r){
    mapply(
      function(i){
        sign(i)
      }, FC[[paste0(r,"_LFC")]]
    )
  }) %>%
    set_names(c("signsFC_CX","signsFC_CB","signsFC_HIP"))
  
  FC_sum <- FC[, c(1,2,4,6)] %>%
    mutate(
      signsFC_CX = sign_cols[[1]],
      signsFC_CB = sign_cols[[2]],
      signsFC_HIP = sign_cols[[3]]
    ) %>%
    as.data.frame()
  
  # data frame for extarcted ensemblIds
  FC_sum_signs <- FC_sum %>%
    mutate(
      sumsigns = apply(FC_sum[, -c(1:4)], 1, sum)
    ) %>%
    mutate(
      commonsignFC = ifelse(abs(sumsigns) == 3, sign(sumsigns), 0)
    ) %>%
    filter(abs(commonsignFC)==1)
  
  #extractions the same direction lfc ensembl Ids in all DGEs
  FC_sum_signs_ens <- FC_sum_signs %>%
    pull(Ensembl_Gene_ID)
  
  # combining the pvalues and filtering for extracted ensembl Ids
  FC_sum_signs_ens_Pvalue <- FC[, c(1,3,5,7)]%>%
    mutate(
      P_COMBINE = mapply(
        function(cx, cb, hip){
          all_p <- c(cx, cb, hip)
          if (all(is.na(all_p))) {
            return(NA)
          }
          notNA <- sum(!is.na(all_p))
          statc <- -2 * sum(log(all_p), na.rm = TRUE)
          return(1 - pchisq(statc, df = (2 * notNA)))
          
        }, CX_P_Val, CB_P_Val, HIP_P_Val
      )
    ) %>%
    mutate(FDR_COMBINE = p.adjust(P_COMBINE, method = "BH")) %>%
    filter(Ensembl_Gene_ID %in% FC_sum_signs_ens) %>%
    full_join(FC_sum_signs, by="Ensembl_Gene_ID") %>%
    left_join(
      read_tsv("Primary_assembly_GNA.txt",  col_names = FALSE) %>%
        setNames(c("Ensembl_Gene_ID", "biotype", "external_gene_name")) %>%
        mutate(
          Ensembl_Gene_ID = gsub("\\..+", "", Ensembl_Gene_ID)
        ),
      by = "Ensembl_Gene_ID"
    )%>%
    arrange(FDR_COMBINE)
  
  writexl::write_xlsx(
    FC_sum_signs_ens_Pvalue, 
    paste0("./output_files/",x,"_regions_Pvalue_combine.xlsx"), col_names = T)

  FC_sum_signs_ens_Pvalue
}) %>%
  writexl::write_xlsx(
    .,
    "./output_files/Single_periods_combined_regions_Pvalue_combine.xlsx", col_names = T)

sheet_names <- readxl::excel_sheets("./output_files/Single_periods_combined_regions_Pvalue_combine.xlsx")

GO_Enrichment_Results <- lapply(set_names(sheet_names), function(y){
  
  Background_Genes <- readxl::read_xlsx(
    "./output_files/Single_periods_combined_regions_Pvalue_combine.xlsx", sheet = y) %>%
    pull(external_gene_name)
  
  Meta_analysis_df <- readxl::read_xlsx(
    "./output_files/Single_periods_combined_regions_Pvalue_combine.xlsx", sheet = y) %>% 
    filter(FDR_COMBINE < 0.1) %>%
    filter(external_gene_name != "Cul3") %>%
    select(external_gene_name, commonsignFC)
  
  All_genes <- Meta_analysis_df %>% 
    pull(external_gene_name) %>%
    unique()
  
  Up_reg_genes <- Meta_analysis_df %>%
    filter(commonsignFC == 1) %>% 
    pull(external_gene_name) %>%
    unique()
  
  Down_reg_gene <- Meta_analysis_df %>%
    filter(commonsignFC == -1) %>% 
    pull(external_gene_name) %>%
    unique()
  
  gene_list <- list(All_genes, Up_reg_genes, Down_reg_gene) %>% 
    set_names(c("ALL","UP","DOWN"))
  
  list_name <- c()
  result_lists <- lapply(set_names(names(gene_list)),function(z){
    sec_result <- lapply(set_names(c("none", "moderate")), function(type){
      temp_result <- gprofiler(query=gene_list[[z]],
                               organism = "mmusculus",
                               correction_method = "fdr",
                               exclude_iea = T,
                               hier_filtering = type,
                               custom_bg = Background_Genes,
                               ordered_query = T,
                               src_filter=c("GO:BP", "GO:MF")) %>%
        arrange(p.value) %>%
        filter(overlap.size>1) %>%
        filter(term.size < 400)
      
      list_name <- c(list_name, paste0("GO_Enrich_",y,"_",z,"_",type))
      temp_result
    }) 
  }) 
  
})

final_outputs <- list()

for (i in names(GO_Enrichment_Results)){
  final_outputs[[i]] = readxl::read_xlsx(
    "./output_files/Single_periods_combined_regions_Pvalue_combine.xlsx", sheet = i)
  for (j in names(GO_Enrichment_Results[[i]])){
    for (k in names(GO_Enrichment_Results[[i]][[j]])){
      final_outputs[[paste(i,j,k,sep = "_")]] = GO_Enrichment_Results[[i]][[j]][[k]]
    }
  }
}
writexl::write_xlsx(final_outputs, 
                    path = "./output_files/GO_Enrich_over_periods_after_meta_analysis.xlsx",
                    col_names = T)










