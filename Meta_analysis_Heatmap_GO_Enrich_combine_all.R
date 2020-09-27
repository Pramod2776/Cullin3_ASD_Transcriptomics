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
  select("Ensembl_ID",everything()) %>%
  left_join(total_exprs_list[[2]], by = "Ensembl_ID") %>%
  left_join(total_exprs_list[[3]], by = "Ensembl_ID")

##### Get TPM of all genes expressed in ########
# three periods within the specific brain region
# total_expr_TPM <- c()
# for (p in three_periods){
#   total_exprs_temp <- total_exprs_list[[p]] %>% 
#     as.data.frame() %>%
#     select(contains(p)) %>%
#     mutate(Ensembl_ID = rownames(total_exprs_list[[p]]))
#   temp_df_list <- c()
#   for (x in three_regions){
#     temp_result <- total_exprs_temp %>%
#       as.data.frame() %>%
#       select(contains(x))
#     temp_df_list <- c(temp_df_list,temp_result)
#   }
#   total_expr_TPM <- c(total_expr_TPM, temp_df_list)
# } 
# total_expr_TPM <- data.frame(total_expr_TPM)
#   


# Extract all expression information of a certain region from CUL3 metadata file 
total_meta_arrange <- total_meta %>%
  mutate(Region = str_extract(total_meta$Sample, "([A-Z]+$)")) %>%
  arrange(Period, Region, Genotype) %>%
  select("Sample", "Period", "Sex", "Genotype", "Region")

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
  select("Ensembl_ID", as.character(final_total_meta$Sample))

# sanity check
cbind((colnames(Match_Total_exprs)[-1]), as.character(final_total_meta$Sample)) %>% 
  as.data.frame()%>%
  setNames(c("ID","Sample")) %>% left_join(total_meta_arrange, by = "Sample")


# All DGE lists
DGE_List_fns <- list.files(path = "./input_files", full.names = T)
DGE_list <- lapply(DGE_List_fns, function(x){
  temp_dge_list <- readxl::read_xlsx(x, sheet = 1) %>%
    filter(FDR<0.1) %>%
    #filter(Biotype=="protein_coding") %>%
    pull(Ensembl_ID...1)
}) %>% 
  unlist() %>%
  unique()

# Extracting unique genes
match_total_exprs_extract <- Match_Total_exprs %>%
  filter(Ensembl_ID %in% DGE_list)

rownames(match_total_exprs_extract) <- 
  match_total_exprs_extract$Ensembl_ID ; match_total_exprs_extract$Ensembl_ID = NULL

# Write out the output
writexl::write_xlsx(match_total_exprs_extract, 
                    "output_files/match_total_exprs_extract.xlsx", 
                    col_names = TRUE)

# Annotation - match gene names and Ensembl_IDs
annotation<-read_tsv("Primary_assembly_GNA.txt", col_names = FALSE) %>%
  setNames(c("ensembl_gene_id", "biotype", "external_gene_name")) %>%
  mutate(
    ensembl_gene_id = gsub("\\..+", "", ensembl_gene_id)
  ) %>%
  as.data.frame()

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


repeat_times <- as.character(rep(1, 18)*6)
column_names <- c()
for (p in three_periods){
  for (r in three_regions){
    for (t in c("WT","HET")){
      column_names <- c(column_names, paste0(p,"_",r,"_",t))
    }
  }
}

annotation_A <- data.frame(
  Genotype = factor(
    rep(column_names,
        repeat_times)
    )
  )
rownames(annotation_A) <- colnames(match_total_exprs_extract_GN)


# Plotting Heatmap and save the file as pdf
pdf("cluster_DGE_0.1_arrange_PC_names_with_rownames_zoomed_in_for_all.pdf", width = 20, height = 180)
pall <- pheatmap((match_total_exprs_extract_GN),  
                 color = colorRampPalette(c("blue", "white", "red"))(41), 
                 cluster_rows = T, cluster_cols = F, 
                 show_rownames = T, show_colnames = T, 
                 scale = "row", 
                 fontsize = 15, fontsize_row = 7, fontsize_col = 4, 
                 annotation = cbind(annotation_A), 
                 border_color = "NA", 
                 breaks = NA, 
                 cellwidth = 5, cellheight = 4.5)
dev.off()


ordered_genes<-rownames(match_total_exprs_extract_GN[pall$tree_row$order,])

# A function to extract genes from ordered gene list based 
# on user input range
extract_cluster_genes<-function(x,y,ordered_genes){
  A<-as.numeric(match(x, ordered_genes))
  B<-as.numeric(match(y, ordered_genes))
  ordered_genes[A[1]:B[1]]
}

# All_genes: the genes fall in the boundary (only for validity check)
gene_list_for_check_for_all <- c("Gm43490","4833418N02Rik","Gm43176","Gm8822","Rps6-ps4","Scarna2",
                         "Lsp1","Mrc1","Ap3s1-ps2","Rpl30-ps10","Hist2h2ab","Ddr2","Ppm1f","Inmt")
# All_genes: the genes fall in the boundary of each cluster
picked_gene_ranges_for_all <- list(
  "E17.5" = list(
    c("Gm43490","4833418N02Rik"),
    c("Gm43176","Gm8822"),
    c("Rps6-ps4","Scarna2")
  ),
  "P7" = list(
    c("Lsp1","Mrc1"),
    c("Ap3s1-ps2","Rpl30-ps10"),
    c("Hist2h2ab","Ddr2")
  ),
  "P35" = list(
    c("Ppm1f","Inmt")
  )
)

# All_genes: the genes fall in the boundary (only for validity check)
gene_list_for_check_for_all_down <- c("Slc25a35","Neat1","Naaa","Emcn","Scp2-ps2","Gm996",
                                      "Ecel1","Tdrd1","Chodl","Gm6311","Rps25-ps1",
                                      "Gm10076","Etv4","Lamb1","Gm6395","Rpl30-ps10")
# All_genes: the genes fall in the boundary of each cluster
picked_gene_ranges_for_all_down <- list(
  "E17.5" = list(
    c("Slc25a35","Neat1"),
    c("Naaa","Emcn"),
    c("Scp2-ps2","Gm996")
  ),
  "P7" = list(
    c()
  ),
  "P35" = list(
    c("Ecel1","Tdrd1"),
    c("Chodl","Gm6311"),
    c("Rps25-ps1","Gm10076"),
    c("Etv4","Lamb1"),
    c("Gm6395","Rpl30-ps10")
  )
)

# Protein_Coding: the genes fall in the boundary (only for validity check)
gene_list_for_check_for_protein_coding <- c("Raver1","Bcor","Xlr3a","Igf2bp3",
                                            "Acp5","Casp4","Gm14434","Tmem254c","Adra2c","Gm9780",
                                            "Sytl2","Ifit3b","Acad10","Acta2","Map6d1","Gm21985")
# Protein_Coding: the genes fall in the boundary of each cluster
picked_gene_ranges_for_protein_coding <- list(
  "E17.5" = list(
    c("Raver1","Igf2bp3")
  ),
  "P7" = list(
    c("Acp5","Casp4"),
    c("Gm14434","Tmem254c"),
    c("Adra2c","Gm9780")
  ),
  "P35" = list(
    c("Sytl2","Ifit3b"),
    c("Acad10","Acta2"),
    c("Map6d1","Gm21985")
  )
)

# Protein_Coding_DownReg
gene_list_for_check_for_protein_coding_Down <- c("Rgcc","Stk39","Lrmp","Gm28845","Lepr","Ehd2","Bglap","Tymp",
                                                 "Prph","Fblim1","Cd151","Col5a2","2210418O10Rik","Insl3",
                                                 "Fbn1","5031439G07Rik","Bcor","Xlr3a")

picked_gene_ranges_for_protein_coding_Down <- list(
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

# When run this chunk of code, change the name of the gene_ranges_list based on what you want
GO_Enrichment_outputs <- lapply(set_names(names(picked_gene_ranges_for_all)),function(x){
  query_genes <- c()
  for (i in picked_gene_ranges_for_all[[x]]){
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
      arrange(p.value)%>%
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
                    path = "GO_Enrich_Up/GO_enrichment_over_periods_for_All_none_Up_regulated_FT.xlsx", 
                    col_names = T)


#!!!!! to select by single genes 
# over regions
#  P7 downregulated exclusively 
