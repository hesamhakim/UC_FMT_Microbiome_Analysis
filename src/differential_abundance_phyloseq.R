library(microbiome)
library(tidyverse)
library(DESeq2)
library('EnhancedVolcano')
library(gridExtra)


setwd("/Users/hhakimjavadi/Dropbox (UFL)/FMT_microbiome_smichail/")
project="FMT"
microbiome_seq_strategy="WGS"
method="Bracken"
psq_list= readRDS(file=paste("objects/",project, "_psq_list_",microbiome_seq_strategy,"_",method,".rds", sep = ""))

#var that being used to seperate sample groups
vars_list<-readRDS(file = paste("objects/vars_list_",project,".rds", sep = ""))
factor_vars_all<-vars_list$vars_track$factor_vars
factor_vars_select<-factor_vars_all[-c(27, 28)]
factor_vars_select<-factor_vars_all


res_table<-data.frame(variable="", ref_group="", compared_to="", Taxa_level="", Number_up="", Number_down="", up_species="", down_species="")
volcano_plot<-list()

tr="original"
taxa_levels<-c("phylum","order","genus", "species")
taxa_levels_select<-taxa_levels[4]

##calculate differentialy abundant species using DESeq2
for (var in factor_vars_select){
  for (tx in taxa_levels_select) {
    psq_analyzed<-psq_list[[tr]][[tx]]
    var_valus<-sample_data(psq_analyzed)[[var]] %>% as.character() %>%as.factor()
    psq_analyzed_clean<-prune_samples(!is.na(var_valus), psq_analyzed)
    ## to prevent the error due to many 0 for many speciment add 1 sudo count to all
    psq_analyzed_clean@otu_table<-psq_analyzed_clean@otu_table +1
    psq_dds = phyloseq_to_deseq2(psq_analyzed_clean, as.formula(paste("~ ",var)))
    psq_dds = DESeq(psq_dds)
    alpha = 0.05
    comp_levels<-levels(var_valus)
    comps<-combn(comp_levels, 2, simplify = FALSE)
    for (C in 1:length(comps)){
      contrast<-c(var,comps[[C]])[c(1,3,2)]
      comparison<-paste(comps[[C]][2], comps[[C]][1], sep = " vs ")
      res = DESeq2::results(psq_dds, contrast = contrast,alpha=alpha)
      res = res[!is.na(res$padj), ]
      res_sig = res[(res$padj < alpha), ]
      sig_up<-res_sig[res_sig$log2FoldChange>0,] %>% rownames()
      sig_down<-res_sig[res_sig$log2FoldChange<0,] %>% rownames()
      title<-paste(comparison, "Taxa Level: ",tx)
      res_table<-rbind(res_table, c(var, comps[[C]][1], comps[[C]][2], tx, length(sig_up), length(sig_down), 
                                    paste(sig_up, collapse=", "), paste(sig_down, collapse = ", ")))
      volcano_plot[[var]][[tx]][[comparison]]<-EnhancedVolcano(res,
                                                            lab = rownames(res),
                                                            x = 'log2FoldChange',
                                                            y = 'padj',
                                                            title = paste("Comparison: ",var,"\n", 
                                                                          "Taxa level: ",tx,"\n",
                                                                          "Case level: ",comps[[C]][2], "\n",
                                                                          "Compared to the reference level: ",comps[[C]][1], "\n",
                                                                          
                                                                          sep = ""),
                                                            titleLabSize = 13,
                                                            labSize=3,
                                                            pCutoff = 0.05,
                                                            FCcutoff = 1,
                                                            caption = NULL,
                                                            subtitle = NULL)
    }
  }
}


for (var in factor_vars_select){
  
  file_name = paste ("figures/differential_abundances/volcanoPlot_",var, ".pdf", sep = "")
  plots<-volcano_plot[[var]] %>% unlist(recursive = FALSE)
  plots_no<-length(plots)
  ncols=5
  nrows=(plots_no %/% ncols)+1
  pdf(width = 30, height =7* nrows)
  ggsave(file = file_name, arrangeGrob(grobs = plots ,ncol = ncols, 
                                       nrow =nrows ), limitsize = FALSE) 
  dev.off()
} 

write_csv(res_table, "tables/significant_up_down_differnetial_abundances_species_for_functional.csv")





# Function to split a long string into a data frame
split_string_to_df <- function(long_string, delimiter = ",") {
  # Split the string based on the provided delimiter
  string_parts <- strsplit(long_string, delimiter)[[1]]
  
  # Trim leading/trailing white spaces from each substring (optional)
  string_parts <- trimws(string_parts)
  
  # Create a one-column data frame
  df <- data.frame(substring = string_parts, stringsAsFactors = FALSE)
  
  return(df)
}


long_sp= res_table %>% filter(Taxa_level=="species") %>% 
  filter(variable=="pucai_range") %>% 
  filter(compared_to=="35_64") %>%
  filter(ref_group=="0_9") %>%
  .$down_species

df=split_string_to_df(long_sp) %>% 
  filter(!str_detect(substring, "sp.")) 

write_csv(df, "tables/species.csv")




