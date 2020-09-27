setwd("/Users/pramod/Desktop/GO_plots_DGE/")
library(tidyverse)
library(doParallel)
registerDoParallel()
BiocManager::install("gProfileR")
library(gProfileR)
library(UniProt.ws)

# BiocManager::install("mygene")
library(mygene)
options(tibble.width = Inf)


dpe_1m_0.1_no_16p<-readxl::read_xlsx("P35_HIP_edgeR_results_S_HIP6_S_GN_level.xlsx", sheet = 1)
dpe_1m_0.1_no_16p<-as.data.frame(dpe_1m_0.1_no_16p)
dpe_1m_0.1_GN_up<-dpe_1m_0.1_no_16p%>%
  subset(FDR <= 0.1 & logFC > 0)

dpe_1m_0.1_GN_down<-dpe_1m_0.1_no_16p%>%
  subset(FDR <= 0.1 & logFC < 0)

GO_1m<-readxl::read_xlsx("Supplementry table-5 all-goterms_combined.xlsx", sheet = 9)

GO_1m_GN<-GO_1m %>%
  separate_rows(intersection)

GO_1m_GN_color<-GO_1m_GN %>%
  mutate (Genes = 
            ifelse (GO_1m_GN$intersection %in% toupper(dpe_1m_0.1_GN_up$Gene_name),
                    yes ="up",
                    no = ifelse(GO_1m_GN$intersection %in% toupper(dpe_1m_0.1_GN_down$Gene_name),
                                yes ="down", 
                                no="none"
                                
                    )))

GO_1m_GN_color_select<-GO_1m_GN_color%>%
  dplyr::select(p.value, term.id, domain, term.name, intersection, Genes)%>%
  mutate(log.p.value = -log10(p.value))

plotting_df1<-GO_1m_GN_color_select %>%  
  group_by(Genes, term.name, p.value, log.p.value) %>% 
  summarise(Freq = n())%>%
  group_by(term.name) %>% 
  mutate(Sum = (sum(Freq)))%>%
  mutate(percent = Freq/(sum(Freq)))%>%
  mutate(Value=log.p.value * percent) %>%
  arrange(desc(Value))
plotting_df1$Genes <- factor(plotting_df1$Genes, levels = c("up","down"))

term.name <- plotting_df1 %>% 
  group_by(term.name) %>%
  summarise(Value = sum(Value)) %>%
  arrange(Value) %>%
  pull(term.name)

#term.name <- unique(plotting_df1$term.name)

p2<-plotting_df1 %>% 
  ggplot(aes(x = term.name , y = (Value), group = Genes, fill = Genes)) +
  # NB: use "y = Freq" instead of "weight = Freq"
  geom_bar(stat = "identity", position = "stack", width = 0.75) +
  coord_flip() +
  scale_x_discrete(limits = rev(term.name))+
  scale_fill_brewer(breaks=c("up", "down"), palette = "Set1")+
  labs(y = "-Log10(FDR)")+
  theme_bw(32)+
  theme(
    text = element_text(size = 32),
    axis.text.y = element_text(size = 32),
    axis.title.y = element_blank()
  )
# guides(
#   size = guide_legend(title = "Overlapping Genes"), 
#   fill = guide_colourbar(title = "Pvalue", override.aes = (list(size = 10))),
#   colour = FALSE
# ) 
  p2
ggsave(
  filename = "Hippocampus-p35-all.pdf",
  plot = p2,
  device = "pdf", width = 20, height = 16
)
