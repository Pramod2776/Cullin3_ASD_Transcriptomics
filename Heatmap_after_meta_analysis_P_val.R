library(dplyr)
library(stringr)
library(tidyverse)
library(writexl)
library(gProfileR) # GO-Enrichment
library(pheatmap) # For Heatmap

# extracting, merging and arranging data
three_regions <- c("CX","CB","HIP")
three_periods <- c("E17_5","P7","P35")


# Reading TPM data for three periods
total_exprs_list <- list(
  total_exprs_E17.5=get(load("RData_files_TPM/CUL3_E17_5.RSEM_Quant.genes.tpm.RData")),
  total_exprs_P7=get(load("RData_files_TPM/CUL3_P7.RSEM_Quant.genes.tpm.RData")),
  total_exprs_P35=get(load("RData_files_TPM/CUL3_P35.RSEM_Quant.genes.tpm.RData"))
) %>% setNames(three_periods)

# Read Metadata file for CUL3
total_meta <- read.table("CUL3_meta.txt", header=T)

total_exprs_list <- lapply(set_names(three_periods), function(x){
  total_exprs_list[[x]] <- total_exprs_list[[x]] %>%
    as.data.frame() %>%
    mutate(Ensembl_ID = rownames(.))
})

total_expr_TPM_raw <- total_exprs_list[[1]] %>%
  dplyr::select("Ensembl_ID",everything()) %>%
  left_join(total_exprs_list[[2]], by = "Ensembl_ID") %>%
  left_join(total_exprs_list[[3]], by = "Ensembl_ID")

# Extract all expression information of a certain region from CUL3 metadata file 
total_meta_arrange <- total_meta %>%
  mutate(Region = str_extract(total_meta$Sample, "([A-Z]+$)")) %>%
  arrange(Period, Region, Genotype) %>%
  dplyr::select("Sample", "Period", "Sex", "Genotype", "Region")

# TPM for a certain region
final_total_meta <- data.frame(
  "Sample" = list(),
  "Period" = list(),
  "Sex" = list(),
  "Genotype" = list(),
  "Region" = list()
)
for (p in three_periods){
  second_total_meta <- data.frame(
    "Sample" = list(),
    "Period" = list(),
    "Sex" = list(),
    "Genotype" = list(),
    "Region" = list()
  )
  for (r in three_regions){
    temp_total_meta <- data.frame(
      "Sample" = list(),
      "Period" = list(),
      "Sex" = list(),
      "Genotype" = list(),
      "Region" = list()
    )
    for (t in c("WT","HET")){
      temp_df <- total_meta_arrange[grep(p, total_meta_arrange$Sample),] 
      temp_df <- temp_df[grep(r, total_meta_arrange$Sample),] %>% filter(Genotype == t)
      temp_total_meta <- rbind(temp_total_meta, temp_df)
    }
    second_total_meta <- rbind(second_total_meta, temp_total_meta) 
  }
  final_total_meta <- rbind(final_total_meta, second_total_meta)
}

Match_Total_exprs <- total_expr_TPM_raw %>%
  dplyr::select("Ensembl_ID", as.character(final_total_meta$Sample))

input_fns <- list(
  over_regions = "./output_files/Single_regions_Combined_periods_Pvalue_combine.xlsx",
  over_periods = "./output_files/Single_periods_combined_regions_Pvalue_combine.xlsx",
  Combined = "./output_files/Combine_All_Pvalue_combine.xlsx"
)

annotation<-read_tsv("Primary_assembly_GNA.txt", col_names = FALSE) %>%
  setNames(c("ensembl_gene_id", "biotype", "external_gene_name")) %>%
  mutate(
    ensembl_gene_id = gsub("\\..+", "", ensembl_gene_id)
  ) %>%
  as.data.frame()

final_results <- lapply(set_names(names(input_fns)),function(f){
  sheets_names <- readxl::excel_sheets(input_fns[[f]])
  temp_result <- lapply(set_names(sheets_names), function(s){
    
   meta_result <- readxl::read_xlsx(input_fns[[f]], sheet = s) %>%
     filter(FDR_COMBINE < 0.1) %>%
     dplyr::select(Ensembl_Gene_ID) %>%
     unique()
   
   # Extracting unique genes
   if (s != "Sheet1"){
     
     match_total_exprs_extract <- Match_Total_exprs %>%
       filter(Ensembl_ID %in% meta_result$Ensembl_Gene_ID) %>%
       dplyr::select(grep(s,colnames(Match_Total_exprs))) 
     ids <- Match_Total_exprs %>%
       filter(Ensembl_ID %in% meta_result$Ensembl_Gene_ID) %>%
       pull(Ensembl_ID)
     match_total_exprs_extract <- match_total_exprs_extract %>%
       mutate(Ensembl_ID = ids)
     
   } else {
     match_total_exprs_extract <- Match_Total_exprs %>%
       filter(Ensembl_ID %in% meta_result$Ensembl_Gene_ID)
   }
   
   rownames(match_total_exprs_extract) <- 
     match_total_exprs_extract$Ensembl_ID ; match_total_exprs_extract$Ensembl_ID = NULL
   
   match_total_exprs_extract_GN<-match_total_exprs_extract %>%
     mutate(
       Gene_symbol = annotation[
         match(rownames(match_total_exprs_extract), 
               annotation$ensembl_gene_id
         ),
         ]$external_gene_name
     ) %>%
     distinct(Gene_symbol, .keep_all = TRUE)
   
   # Set rowname to Gene names and delete the Gene names column
   rownames(match_total_exprs_extract_GN) <- 
     match_total_exprs_extract_GN$Gene_symbol
   match_total_exprs_extract_GN$Gene_symbol=NULL
   
   write.table(match_total_exprs_extract_GN,paste0("./80_Genes_heatmap/heatmap_input","_",f,"_",s,".tsv"),
               col.names = T,
               row.names = T,
               sep = "\t"
               )
   
   # sample_size <- 0
   # if (s == "Sheet1"){
   #   sample_size <- 18
   #   repeat_times <- as.character(rep(1, sample_size)*6)
   #   column_names <- c()
   #   for (p in three_regions){
   #     for (r in three_periods){
   #       for (t in c("WT","HET")){
   #         column_names <- c(column_names, paste0(p,"_",r,"_",t))
   #       }
   #     }
   #   }
   # } else {
   #   sample_size <- 6
   #   repeat_times <- as.character(rep(1, sample_size)*6)
   #   column_names <- c()
   #   if (s %in% three_regions){
   #     for (r in three_periods){
   #       for (t in c("WT","HET")){
   #         column_names <- c(column_names, paste0(s,"_",r,"_",t))
   #       }
   #     }
   #   } else {
   #     for (p in three_regions){
   #       for (t in c("WT","HET")){
   #          column_names <- c(column_names, paste0(p,"_",s,"_",t))
   #       }
   #     }
   #   }
   # }
   # 
   # 
   # annotation_A <- data.frame(
   #   Genotype = factor(
   #     rep(column_names,
   #         repeat_times)
   #   )
   # )
   # rownames(annotation_A) <- colnames(match_total_exprs_extract_GN)
   # 
   # pdf(paste0("./final/cluster_P_Val_Comb_0.1_arrange_PC_names_with_rownames_zoomed_in_for_"
   #             ,f,"_",s,"_",".pdf")
   #      , width = 15, height = 10)
   # pall <- pheatmap((match_total_exprs_extract_GN),
   #                  color = colorRampPalette(c("blue", "white", "red"))(41),
   #                  cluster_rows = T, cluster_cols = F,
   #                  show_rownames = T, show_colnames = F,
   #                  scale = "row",
   #                  fontsize = 15, fontsize_row = 7, fontsize_col = 4,
   #                  annotation = cbind(annotation_A),
   #                  border_color = "NA",
   #                  breaks = NA,
   #                  cellwidth = 5, cellheight = 6)
   # dev.off()
   # 
   # ordered_genes<-rownames(match_total_exprs_extract_GN[pall$tree_row$order,])
   
   
   
  })
})
save(final_results, file = "./RData_Results/Heatmap_Genes.RData")


# A function to extract genes from ordered gene list based 
# on user input range
extract_cluster_genes<-function(x,y,ordered_genes){
  A<-as.numeric(match(x, ordered_genes))
  B<-as.numeric(match(y, ordered_genes))
  ordered_genes[A[1]:B[1]]
}

gene_list_for_check_for_E17_5_Down <- c()

picked_gene_ranges_for_E17_5_Down <- list(
  "E17.5" = list(
    c("Rgcc","Stk39"),
    c("Lrmp","Gm28845"),
    c("Lepr","Ehd2"),
    c("Bglap","Tymp")
  ),
  "P7" = list(
    c()
  ),
  "P35" = list(
    c("Prph","Fblim1"),
    c("Cd151","Col5a2"),
    c("2210418O10Rik","Insl3"),
    c("Fbn1","5031439G07Rik"),
    c("Bcor","Xlr3a")
  )
)

gene_list_for_check_for_E17_5_Up <- c()

picked_gene_ranges_for_E17_5_Up <- list(
  "E17.5" = list(
    c("Rgcc","Stk39"),
    c("Lrmp","Gm28845"),
    c("Lepr","Ehd2"),
    c("Bglap","Tymp")
  ),
  "P7" = list(
    c()
  ),
  "P35" = list(
    c("Prph","Fblim1"),
    c("Cd151","Col5a2"),
    c("2210418O10Rik","Insl3"),
    c("Fbn1","5031439G07Rik"),
    c("Bcor","Xlr3a")
  )
)

gene_list_for_check_for_Comb_Down <- c("Slc11a1","Vtn",
                                       "Flna","Mgst1")

picked_gene_ranges_for_Comb_Down <- list(
  "E17.5" = list(
    c("Slc11a1","Vtn")
  ),
  "P7" = list(
    c()
  ),
  "P35" = list(
    c("Flna","Mgst1")
    
  )
)

heatmap_genes <- final_results[["Combined"]][["Sheet1"]]

GO_Enrichment_outputs <- lapply(set_names(names(picked_gene_ranges_for_Comb_Down)),function(x){
  query_genes <- c()
  for (i in picked_gene_ranges_for_Comb_Down[[x]]){
    if (length(i) != 0){
      temp_genes <- extract_cluster_genes(i[[1]],i[[2]],ordered_genes = ordered_genes)
      query_genes <- c(query_genes, temp_genes)
    }
    else{
      break
    }
  }
  
  dge_fns_period <- list.files(path = "./input_files", pattern = x, full.names = T)

  dge_list_bg <- lapply(dge_fns_period, function(f){
    temp_dge <- readxl::read_xlsx(f, sheet = 1) %>%
      filter(!Gene_name=="Cul3")%>%
      pull(Gene_name)
  }) %>%
    unlist() %>%
    unique()
  
  if (length(query_genes) != 0){
    cluster_moderate<-gprofiler(query=query_genes,
                                organism = "mmusculus",
                                correction_method = "fdr",
                                exclude_iea = T,
                                hier_filtering = "none",
                                #custom_bg = rownames(match_total_exprs_extract_GN),
                                custom_bg = dge_list_bg,
                                ordered_query = T,
                                src_filter=c("GO:BP", "GO:MF")) %>%
      arrange(p.value) %>%
      filter(overlap.size>1) %>%
      filter(term.size < 400)
    
    # writexl::write_xlsx(cluster_moderate, 
    #                     path = paste0("GO_Enrich_for_DownReg/GO_enrichment_",x,"_periods_none_.xlsx"),
    #                     col_names = T)
  }
  else{
    cluster_moderate <- data.frame()
  }
  cluster_moderate
}) %>% 
  writexl::write_xlsx(., 
                      path = "GO_Enrich/GO_enrichment_all_periods_for_All_none_Down_regulated_FT.xlsx", 
                      col_names = T)






