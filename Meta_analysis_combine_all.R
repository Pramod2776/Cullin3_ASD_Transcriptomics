library(tidyverse)
library(dplyr)

all_fns <- list.files(path = "./input_files", full.names = T)
p_n_r <- c()
for (fn in all_fns){
  sub_fn <- str_split(str_split(fn, "/")[[1]][3],"_")
  period_region <- paste(sub_fn[[1]][1],sub_fn[[1]][2], sep = "_")
  p_n_r <- c(p_n_r, period_region)
}

all_fns <- set_names(all_fns, p_n_r)

DGEs <- lapply(set_names(p_n_r), function(f){
  print(f)
  temp_dge <- readxl::read_xlsx(all_fns[[f]], sheet = 1) %>%
    select(Ensembl_ID...1, logFC, PValue) %>% 
    set_names(
      c("Ensembl_Gene_ID", paste0(f,"_LFC"), paste0(f,"_P_Val"))
    )
})

FC <- DGEs[[1]] 
for (e in seq(2,length(DGEs))){
  FC <- full_join(FC,DGEs[[e]])
} 
  
sign_cols <- lapply(p_n_r, function(p){
  mapply(
    function(i){
      sign(i)
    }, FC[[paste0(p,"_LFC")]]
  )
}) 

FC_sum <- FC[, c(1,seq(2,18,2))] %>%
  mutate(
    signsFC_E17.5_CB = sign_cols[[1]],
    signsFC_E17.5_CX = sign_cols[[2]],
    signsFC_E17.5_HIP = sign_cols[[3]],
    signsFC_P7_CB = sign_cols[[4]],
    signsFC_P7_CX = sign_cols[[5]],
    signsFC_P7_HIP = sign_cols[[6]],
    signsFC_P35_CB = sign_cols[[7]],
    signsFC_P35_CX = sign_cols[[8]],
    signsFC_P35_HIP = sign_cols[[9]]
  ) %>%
  as.data.frame()

# data frame for extarcted ensemblIds
FC_sum_signs <- FC_sum %>%
  mutate(
    sumsigns = apply(FC_sum[, -c(1:10)], 1, sum)
  ) %>%
  mutate(
    commonsignFC = ifelse(abs(sumsigns) == 9, sign(sumsigns), 0)
  ) %>%
  filter(abs(commonsignFC)==1)

#extractions the same direction lfc ensembl Ids in all DGEs
FC_sum_signs_ens <- FC_sum_signs %>%
  pull(Ensembl_Gene_ID)

# combining the pvalues and filtering for extracted ensembl Ids
FC_sum_signs_ens_Pvalue <- FC[, c(1,seq(3,19,2))]%>%
  mutate(
    P_COMBINE = mapply(
      function(e17.5_cb, e17.5_cx, e17.5_hip,
               p7_cb, p7_cx, p7_hip,
               p35_cb, p35_cx, p35_hip){
        all_p <- c(e17.5_cb, e17.5_cx, e17.5_hip,
                   p7_cb, p7_cx, p7_hip,
                   p35_cb, p35_cx, p35_hip)
        if (all(is.na(all_p))) {
          return(NA)
        }
        notNA <- sum(!is.na(all_p))
        statc <- -2 * sum(log(all_p), na.rm = TRUE)
        return(1 - pchisq(statc, df = (2 * notNA)))
        
      }, E17.5_CB_P_Val,E17.5_CX_P_Val,E17.5_HIP_P_Val,
      P7_CB_P_Val, P7_CX_P_Val, P7_HIP_P_Val,
      P35_CB_P_Val, P35_CX_P_Val, P35_HIP_P_Val
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
  paste0("./output_files/Combine_All_Pvalue_combine.xlsx"), col_names = T)


sheet_names <- readxl::excel_sheets("./output_files/Combine_All_Pvalue_combine.xlsx")

GO_Enrichment_Results <- lapply(set_names(sheet_names), function(y){
  
  Background_Genes <- readxl::read_xlsx(
    "./output_files/Combine_All_Pvalue_combine.xlsx", sheet = y) %>%
    pull(external_gene_name)
  
  Meta_analysis_df <- readxl::read_xlsx(
    "./output_files/Combine_All_Pvalue_combine.xlsx", sheet = y) %>% 
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
  final_outputs[["Combine_All"]] = readxl::read_xlsx(
    "./output_files/Combine_All_Pvalue_combine.xlsx", sheet = i)
  for (j in names(GO_Enrichment_Results[[i]])){
    for (k in names(GO_Enrichment_Results[[i]][[j]])){
      final_outputs[[paste(j,k,sep = "_")]] = GO_Enrichment_Results[[i]][[j]][[k]]
    }
  }
}
writexl::write_xlsx(final_outputs, 
                    path = "./output_files/GO_Enrich_combine_all_after_meta_analysis.xlsx",
                    col_names = T)

# cpm