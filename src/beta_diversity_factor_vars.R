#library(microbiome)
#library(knitr)

library(tidyverse)
library(microbiome)
library(ggpubr)
library(phyloseq)
library(vegan)
library(gridExtra)

################# exploratory data anlaysis: beta diversity, permanova, pca #######


#cd /blue/chamala/he.hakimjavadi/projects/TCGA_metagenomics/Poore2020/
#srundev --time=420 --ntasks=1 --cpus-per-task=16 --mem=48gb
# module load R  
#setwd("/blue/chamala/he.hakimjavadi/projects/TCGA_metagenomics/Poore2020")

setwd("/Users/hhakimjavadi/Dropbox (UFL)/FMT_microbiome_smichail/")
project="FMT"
microbiome_seq_strategy="WGS"
method="Bracken"
psq_list= readRDS(file=paste("objects/",project, "_psq_list_",microbiome_seq_strategy,"_",method,".rds", sep = ""))


####PERMANOVA using psq_clr vs all factor variables

vars_list<-readRDS(file = paste("objects/vars_list_",project,".rds", sep = ""))
factor_vars_all<-vars_list$vars_track$factor_vars
continues_vars<-vars_list$vars_track$continues_vars
vars_select<-union(factor_vars_all, continues_vars)

permanova_tbl<-data.frame("variable"=1,"R2"=1, "p.val"=1,"taxa_level"=1, "abundance_data_transformation"=1,"q.val"=1)
taxa_levels<-c("phylum","order","genus", "species")
taxa_levels_select=taxa_levels[c(2:4)]
psq_transofrms<-c("original", "hlg", "cmp")
select_transformations=c("cmp", "hlg")
innit_psq<-psq_list[[select_transformations[1]]][[taxa_levels_select[1]]]
ommi_single_lvl_var<-c()
permutations=5000
dist_method<-"jsd"#bray or jsd wunifrac
for (tx in taxa_levels_select) {
  for (tr in select_transformations) {
    for (i in which(sample_variables(innit_psq) %in% vars_select)) {
      psq_analyzed<-psq_list[[tr]][[tx]]
      psq_analyzed<-subset_samples(psq_analyzed, sample_data(psq_analyzed)[,i]!="NA")
      levels.no<-sample_data(psq_analyzed)[[i]] %>% as.factor() %>% levels() %>% length()
      if (levels.no>1) {
        dist_matrix <- phyloseq::distance(psq_analyzed, method = dist_method)
        permanova<-vegan::adonis2(dist_matrix ~ sample_data(psq_analyzed)[[i]],parallel = 8, permutations = permutations)
        
        permanova_tbl<-rbind(permanova_tbl, c(sample_variables(psq_analyzed)[i], permanova[["R2"]][1], permanova[["Pr(>F)"]][1], tx, tr))
      }
      else{
        ommi_single_lvl_var<-c(ommi_single_lvl_var, sample_variables(psq_analyzed)[[i]])
      }
    }
    ### adjust for multiple comparisons
    permanova_tbl$p.val=as.numeric(permanova_tbl$p.val)
    permanova_tbl=mutate(permanova_tbl, q.val=p.val*length(vars_select))
  }
}

permanova_tbl[permanova_tbl$q.val>0.05, "q.val"]="NS" ## this will be used in the pca figs
write.csv(permanova_tbl, file = (paste("tables/PERMANOVA_res_all_vars_",project,"_",microbiome_seq_strategy,"_",method,".csv", sep = "")))

###identify significant variables

permanova_tbl_sig<-filter(permanova_tbl, q.val<=0.05)

write.csv(permanova_tbl_sig, file = (paste("tables/PERMANOVA_res_sig_vars_",project,"_",microbiome_seq_strategy,"_",method,".csv", sep = "")))
sig_vars_beta<-permanova_tbl_sig$variable %>% unique()

#pca of data colored by all variables with significant beta diversity
psq_transofrms<-c("original", "hlg", "cmp")
select_transformations<-c("hlg")
ord_method<-"PCoA" #or "PCoA"
factor_vars_select<-sig_vars_beta
pca_plots<-list()
for (tx in taxa_levels_select) {
  for (tr in select_transformations) {
    for (i in vars_select){
      psq_analyzed<-psq_list[[tr]][[tx]]
      psq_analyzed<-subset_samples(psq_analyzed, !is.na(sample_data(psq_analyzed)[,i])) #to get rid os NAs in pca
      dist <- phyloseq::distance(psq_analyzed, method=dist_method, na.rm=FALSE)
      ord <- ordinate(psq_analyzed, method = ord_method, distance = dist)
      evals <- ord$values$Eigenvalues
      qval=filter(permanova_tbl, variable==i & taxa_level==tx & abundance_data_transformation==tr)$q.val
      p1<-plot_ordination(psq_analyzed, ordination = ord, color = i) +
        labs(col = i) +
        stat_ellipse(level = 0.8)+
        ggtitle(paste(project,"   Distance Method: ", dist_method, "\n",
                      "Ordination Method: ", ord_method, "\n", 
                      "Transformation: ", tr, "\n",
                      "Taxa Level: ", tx,"\n",
                      "q.val= ", qval))+
        theme(legend.position="none")+
        coord_fixed() #use this only for NMDS ordering method, otherwise use blow:
      #coord_fixed(sqrt(evals[2] / evals[1]))
      pca_plots[[i]][[tr]][[tx]]<-p1
    }
  }
}

saveRDS(psq_list, paste("objects/PERMANOVA_pcas_",project, "_",microbiome_seq_strategy,"_",method,".rds", sep = ""))

#### plot all pcas
#### function to extract a legend
get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

plots01<-list()
plot01<-list()

for (i in factor_vars_select){
  for (tr in select_transformations) {
    file_name = paste ("figures/pca/beta_divers_",i, "_", ord_method,"_",tr, "_",
                       dist_method ,"_",microbiome_seq_strategy,"_",method, ".pdf", sep = "")
    plots_all<-ggarrange(plotlist = pca_plots[[i]][[tr]], ncol=length(taxa_levels_select), nrow=1, common.legend = TRUE, legend="bottom")
    pdf(width = 12, height = 6)
    ggsave(file = file_name,plots_all , limitsize = FALSE) 
    
    dev.off()
    
  }
  
} 

plots_all<-plot01<- pca_plots[[i]][[tr]]
plots_all_1<-plot01<- pca_plots[[i]][[tr]][[1]]
plots_all_2<-plot01<- pca_plots[[i]][[tr]][[2]]
plots_all_3<-plot01<- pca_plots[[i]][[tr]][[3]]


