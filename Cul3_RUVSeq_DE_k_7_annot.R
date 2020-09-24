getwd()
setwd("/Users/pramod/Desktop/Clu3/Cul3_gender_interaction/")
library(limma)
BiocManager::install("Rcpp")
install.packages("edgeR", dependencies = TRUE)
update.packages(repos='http://cran.rstudio.com/', ask=FALSE, checkBuilt=TRUE)
library(Rcpp)
library(edgeR)
library(EDASeq)
library(RUVSeq)
library(ffpe)
library(RColorBrewer)
library(tidyverse)
library(ggplot2)


# Function Name: concise_RUVSeq
# 
# Parameters: 
#         metatable_name: The file name of the metadata sheet, which contains Samples' ID, Gender,
#                         Genotype (A, B) and Brain Regions (X,Y,Z).
#         period: The period of the input sample.
#         region: The region of the input sample.
#         RUV_input_filename: count matrix after mapping to mouse genome.
#         negctrl_filename: The file name of the file which contains negative control genes.
#         num_k: The number of differnent k's that the user want to try on. (Ex. If the input is 3, this 
#               function will run RUVSeq for k = 1,2,3)
#
concise_RUVSeq <- function(metatable_name, period, region, 
                           negctrl_filename, num_k){
  
  if(!dir.exists(paths = paste0("./Plots/RUVSeq/", paste0(period, "_", region)))){
    dir.create(path = paste0("./Plots/RUVSeq/", paste0(period, "_", region)), recursive = T)
  }
  if(!dir.exists(paths = paste0("./Results/RUVSeq/", paste0(period, "_", region)))){
    dir.create(path = paste0("./Results/RUVSeq/", paste0(period, "_", region)), recursive = T)
  }
  
  metatable <- read.table(metatable_name, header = T)
  metatable <- metatable %>%
    mutate(Region = gsub(".*[0-9]_", "", Sample)) %>%
    mutate(Period = as.character(Period)) %>%
    mutate(Period = ifelse(Period == "E17_5", yes = "E17.5", no = Period)) %>%
    filter(Period == period & Region == region) %>%
    column_to_rownames(var = "Sample")
  
  Ruv_input <- read.table(paste0("./data/Counts/", period, 
                                 "_TPM_0.1_80_", region, 
                                 "HETWT_sorted_counts.txt"), 
                          header = T, row.names = 1)
  negControls <- read.table(negctrl_filename, sep = "\t", header = TRUE, as.is = TRUE)
  
  Region_Sex_new <- metatable[order(metatable$Sex), ] %>%
    dplyr::select(Sex) 
  
  Ruv_input <- as.matrix(Ruv_input[ ,match(rownames(Region_Sex_new), colnames(Ruv_input))])
  metatable <- metatable[match(rownames(Region_Sex_new), rownames(metatable)),]
  
  neg_Con <- intersect(negControls[,2], rownames(Ruv_input))
  
  colors <- brewer.pal(9, "Set1")
  colLib <- colors[as.factor(Region_Sex_new$Sex)]
  
  uq_region <- betweenLaneNormalization(Ruv_input, which = "upper")
  pdf(paste0("./Plots/RUVSeq/", paste0(period, "_", region, "/"),
             "RUVSeq_uq_",region,"_new_font_size.pdf"),
      width = 16, height = 12)
  plotRLE(uq_region, col = colLib, outline=FALSE, las = 3, ylim = c(-.2, .2),
          ylab = "Relative log Expression", cex.axis = 0.6, cex.lab = 1.5)
  plotPCA(uq_region, col = colLib, cex=2,cex.axis = 1.5,
          cex.lab = 1.5, xlim = c(-.6, .9), ylim = c(-.7,.6) )
  dev.off()
  
  num_female <- Region_Sex_new %>%
    filter(Sex == "F") %>%
    pull(Sex) %>%
    length()
  
  num_male <- Region_Sex_new %>%
    filter(Sex == "M") %>%
    pull(Sex) %>%
    length()
  
  sex_diff <- abs(num_female - num_male)
  
  if(num_female > num_male){
    groups_region_new <- matrix(data = c(c(1:num_female),
                                         c((num_female+1):12),rep(-1, times = sex_diff)), 
                                nrow = 2, byrow = T)
  } else {
    groups_region_new <- matrix(data = c(c(1:num_female, rep(-1, times = sex_diff)),
                                         c((num_female+1):12)), 
                                nrow = 2, byrow = T)
  }
  
  S_Regions <- list("uq_region" = uq_region)
  for (i in seq(num_k)) {
    S_Region <- RUVs(uq_region, neg_Con, k=i, groups_region_new)
    S_Regions <- append(S_Regions,list(S_Region))
    pdf(paste0("./Plots/RUVSeq/", paste0(period, "_", region, "/"),"RUVSeq_",
               period,"_",region,"_k_",i,"_new_font_size.pdf"),
        width = 16, height = 12)
    plotRLE(S_Region$normalizedCounts, col = colLib,outline = FALSE,las = 3,ylim = c(-.2, .2),
            ylab = "Relative log Expression", cex.axis = 0.6, cex.lab = 1.5)
    plotPCA(S_Region$normalizedCounts, col = colLib, cex = 1.5,cex.axis = 1.5,
            cex.lab = 1.5, xlim = c(-.6, .9), ylim = c(-.7,.6) )
    dev.off()
    
  }
  S_Regions <- set_names(S_Regions, nm = c("uq_region", seq(num_k)))
  
  options <- expand.grid(names(S_Regions),c(TRUE,FALSE))
  
  for (i in seq(length(options[[1]]))){
    
    if (options[i,][[2]] & options[i,][[1]] != "uq_region"){
      design_new <- model.matrix(~Region_Sex_new$Sex + S_Regions[[as.character(options[i,][[1]])]]$W + 
                                   metatable$Genotype+ Region_Sex_new$Sex:metatable$Genotype)
    } else if(options[i,][[2]] & options[i,][[1]] == "uq_region") {
      design_new <- model.matrix(~Region_Sex_new$Sex + metatable$Genotype+Region_Sex_new$Sex:metatable$Genotype)
    } else if(options[i,][[1]] != "uq_region"){
      design_new <- model.matrix(~Region_Sex_new$Sex + S_Regions[[as.character(options[i,][[1]])]]$W)
    } else {
      design_new <- model.matrix(~Region_Sex_new$Sex)
    }
    
    
    y <- DGEList(counts=Ruv_input, group=as.factor(Region_Sex_new$Sex)) %>% 
      calcNormFactors(., method="upperquartile") %>% 
      estimateGLMCommonDisp(., design_new, verbose=TRUE) %>% 
      estimateGLMTagwiseDisp(., design_new)
    
    fit <- glmFit(y, design_new)
    lrt <- glmLRT(fit, coef=2)
    topFC <- topTags(lrt, n=Inf)$table
    annotation <- read.table("./data/Metadata/Primary_assembly_GNA.txt") %>%
      mutate("Ensembl_ID" = gsub("[//.][0-9]*", "", V1),
             "Biotype" = V2,
             "External_Gene_name" = V3) %>%
      dplyr::select(c("Ensembl_ID", "Biotype", "External_Gene_name"))
    Final_DE_Result <- left_join(topFC %>%
                                   rownames_to_column(var = "Ensembl_ID"),
                                 annotation, by = "Ensembl_ID")
    writexl::write_xlsx(Final_DE_Result, path = if_else(options[i,][[2]],
                                                        true = paste0("./Results/RUVSeq/", paste0(period, "_", region,"/"),period, "_", region, 
                                                                      "_", options[i,][[1]],"_S_", "edgeR_results_level_by_Sex.xlsx"),
                                                        false = paste0("./Results/RUVSeq/", paste0(period, "_", region,"/"),period, "_", region, 
                                                                       "_", options[i,][[1]], "_edgeR_results_level_by_Sex.xlsx")))
    
    
    pdf(paste0("./Plots/RUVSeq/", paste0(period, "_", region, "/"),
               "RUVSeq_",period,"_",region, options[i,][[1]],"p_val_FDR_new_font_size.pdf"),
        width = 12, height = 16)
    hist(topFC$PValue, main="", xlab="p-value", breaks=50, ylim=c(0, 1400), 
         cex = 1.5, cex.axis = 1.5, cex.lab = 1.5)
    hist(topFC$FDR, main="", xlab="FDR", breaks=50, ylim=c(0, 1400),
         cex = 1.5, cex.axis = 1.5, cex.lab = 1.5)
    
    plot(topFC[,1], -log10(topFC$PValue), pch=20, col="gray", cex=1.5,
         ylab="-log10(p-value)", xlab="log2(FC)", ylim=c(0, 85), 
         xlim=c(-2, 4), cex.lab=2, cex.axis=2)
    de <- rownames(topFC[topFC$FDR<=0.01,])
    points(topFC[de,1], -log10(topFC[de, "PValue"]), pch=20, col=colors[2], cex=2)
    dev.off()
    
  }
  
}



for (p in c("E17.5", "P7", "P35")){
  for(r in c("CX", "CB", "HIP")){
    
    concise_RUVSeq("./data/Metadata/CUL3_meta.txt",
                   p,
                   r,
                   "./data/Metadata/Peixoto_NegativeControls_mouse_HK.txt",
                   7)
    
  }
}

