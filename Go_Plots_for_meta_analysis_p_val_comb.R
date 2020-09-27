library(tidyverse)
library(doParallel)
registerDoParallel()
library(gProfileR)
library(UniProt.ws)
library(mygene)
options(tibble.width = Inf)


Go_terms_fns <- list.files(path = ".",pattern = "selected_Go_terms") %>%
  set_names(c("Combined_all","Over_periods","Over_regions"))

for (f in names(Go_terms_fns)){
  sheet_names <- readxl::excel_sheets(Go_terms_fns[[f]])
  
  pdf(paste0("Go_Plots_for_",f,"_after_meta_analysis_P_val_Comb.pdf"), width = 50, height = 16)
  for (s in sheet_names){
    
    go_enrich_results <- readxl::read_xlsx(Go_terms_fns[[f]],sheet = s)
    
    if (length(colnames(go_enrich_results)) != 0){
      GO_file_color_select <- go_enrich_results %>% 
        mutate(log.p.value = -log10(p.value)) %>%
        arrange(desc(log.p.value))
    } else {
      GO_file_color_select <- data.frame()
    }
    
    
    if (nrow(GO_file_color_select) != 0){
      output_plot <- ggplot(
        GO_file_color_select,
        aes(
          x = reorder(str_wrap(term.name, 40), log.p.value), y = log.p.value
        )
      ) +
        ggtitle(s)+
        geom_bar(stat = "identity", position = "dodge", width = 0.75, fill = "blue") +
        coord_flip() +
        labs(y = "-Log10(FDR)")+
        theme_bw()+
        theme(
          plot.title = element_text(size = 32, face = "bold"),
          text = element_text(size = 32),
          axis.text.y = element_text(size = 32),
          axis.title.y = element_blank()
        )
      print(output_plot)
    }
  }
  dev.off()
}
