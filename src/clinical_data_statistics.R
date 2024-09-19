######### alpha diversity of samples grouped acccoriding to all different variables
# to examine diversity we use rarefied data without removing low prevalent and low abandant spicies

#Overal strategy
# 1. create a alpha diversity tables, using 3 diversity index
# 1.1 add all factor vars to the diversity table
# 2. create a table of statistics and plates (remove outliers)
# 3. identify variables resulted in significant difference
# 4. plot significant variables



library(microbiome)
library(knitr)
library(dplyr)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(ggstatsplot)
library(gridExtra)

setwd("/Users/hhakimjavadi/Dropbox (UFL)/FMT_microbiome_smichail/")
project="FMT"
microbiome_seq_strategy="WGS"
method="Bracken"

psq_list= readRDS(file=paste("objects/",project, "_psq_list_",microbiome_seq_strategy,"_",method,".rds", sep = ""))
metadata<-read.csv(paste("data/metadata/metadata_",project,"_",microbiome_seq_strategy,"_processed",".csv", sep = ""))
vars_list<-readRDS(file = paste("objects/vars_list_",project,".rds", sep = ""))

factor_vars=c("Treatment_Timepoints_string","Placebo_Timepoints_string")
factor_vars_select=factor_vars[1]
continues_var="pucai"
data=select(metadata, Subject_ID, Sample_ID, factor_vars_select, continues_var)
data=na.omit(data)
levels<-data[[factor_vars_select]] %>% as.character() %>% as.factor() %>%levels()



plot_list<-list()
ommit_plots<-c()
statistics<-c("statistic", "p.value", "method","conf.low", "conf.high")

plot<-ggbetweenstats(
  data = data,
  x = !!factor_vars_select,
  y = pucai,
  type = "Bayes Factor",
  plot.type = "violin",
  #ggsignif.args = list(textsize = 1, tip_length = 0.01),
  pairwise.comparisons = TRUE,
  title = paste0("PUCAI Score Within Treatment/Placebo Group During the Treatment"),
  p.adjust.method = "holm",
  outlier.color = "grey",
  caption = TRUE,
  results.subtitle = TRUE
)

p<-plot + 
  theme(plot.margin = margin(1, 1, 1, 1, "cm"))

stat<-plot %>% extract_stats()
stat<-as.data.frame(stat[["subtitle_data"]][,statistics])



write.csv(stat, (paste("tables/pucai_vs_placebo_timepoint_stats.csv", sep = "")))
pdf(file = "figures/statistics/pucai_vs_treatment_timepoint_stats.pdf", height = 5.5, width = 6.5)
p
dev.off()





