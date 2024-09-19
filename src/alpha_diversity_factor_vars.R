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




############ diversity index #####################

#a comprehensive list of  indicators of diversity 
div_index_names=c('observed', 'chao1', 'diversity_inverse_simpson', 'diversity_gini_simpson', 'diversity_shannon', 'diversity_fisher', 
                  'diversity_coverage', 'evenness_camargo', 'evenness_pielou', 'evenness_simpson', 'evenness_evar', 'evenness_bulla', 
                  'dominance_dbp', 'dominance_dmn', 'dominance_absolute', 'dominance_relative', 'dominance_simpson', 'dominance_core_abundance', 
                  'dominance_gini', 'rarity_log_modulo_skewness', 'rarity_low_abundance', 'rarity_rare_abundance')
div_index_select<-c("observed","chao1","diversity_shannon","evenness_simpson",
              "dominance_simpson","rarity_rare_abundance")



vars_list<-readRDS(file = paste("objects/vars_list_",project,".rds", sep = ""))
factor_vars<-vars_list$vars_track$factor_vars
continuous_vars<-vars_list$vars_track$continues_vars
vars_select<-union(factor_vars, continuous_vars)

#diversity_tbl[,factor_vars]=1


#1. create a table of diversity indexes for all samples
#diversity index is calculated for all taxa levels
taxa_levels=c("phylum","order","genus", "species")
taxa_levels_select<-c("phylum","order","genus", "species")
transformations=c("original", "cmp", "rrfd", "hlg")
select_transformations=c("original", "rrfd")
alpha_div_tbl<-data.frame(data.frame(matrix(nrow = 0, ncol = length(vars_select)+3)))
colnames(alpha_div_tbl)<-c("smaple_id","data_transformation","taxa",vars_select)
require(ggpubr)
for (tx in taxa_levels_select) {
  for (tr in select_transformations) {
    psq_analyzed=psq_list[[tr]][[tx]]
    df<-microbiome::alpha(psq_analyzed, index = div_index_select)
    #1.1 add all variables data to the div table:
    df[,vars_select]=sample_data(psq_analyzed)[, vars_select]
    df[,"Sample_ID"]=sample_data(psq_analyzed)[, "Sample_ID"]
    df[,"Subject_ID"]=sample_data(psq_analyzed)[, "Subject_ID"]
    df[,"taxa"]<-tx
    df[, "data_transformation"]<-tr
    #df<-rownames_to_column(df, "sampe_id")
    alpha_div_tbl<-rbind(alpha_div_tbl, df)
  }

}

alpha_div_tbl_long<-pivot_longer(alpha_div_tbl, cols = all_of(div_index_select), names_to = "diversity_index_name", values_to = "diversity_value")
write_csv(alpha_div_tbl_long, paste0("tables/alpha_diversities_all_samples_", project, "_", microbiome_seq_strategy,"_",method, ".csv"))


# 2. create a table of statistics for diversity (see two method 2.1 and 2.2)

### 2.1: meothod 1: using kruskal.test function from R (Not recommended see the next method)
####### because the result of this test is not compatible wit the ggbetween test that is being used for plotting
####### 2.1.1 for each taxa/variable/index create a table by filtering out outlier and single level variable
stat_params<-c("parameter1", "parameter2", "statistic","p.value","conf.low", "conf.high")
stats<-data.frame()
ommit_var<-c()
for (tx in taxa_levels_select) {
  for (var in factor_vars) {
    for (i in div_index_select){
      data<-alpha_div_tbl %>% filter(taxa==tx)
      #data[[i]]<-remove_outliers(data[[i]], 5, na.rm = TRUE)
      data<-filter(data, (!is.na(data[[i]])))
      levels<-data[[var]] %>% as.character() %>% as.factor() %>%levels()
      if (length(levels)>1) {
        formula<-paste(i, "~", var, sep="") %>% as.formula()
        stat<-kruskal.test(formula, data = data)
        stat<-data.frame("taxa"=tx, "variable"=var, "index"=i, "p.value"=stat$p.value, 
                         "statistics"=stat$statistic, "group.no"=length(levels), "obs"=length(data[[i]]))
        stats<-rbind(stats, stat)
      } else {
        ommit_var<-c(ommit_var, var) %>% unique()
      }
    }
  }
}



# 2. create a table of statistics for diversity
### 2.2: meothod 2: using ggbetweenstats function and extract statistical data only(recommended method)
####### in this method we calculate the statistics for all variable and we plot significant only in the next step
####### statistics results and plots are compatible  because both generated using the same function "ggbetweenstats"
####### 2.2.1 for each taxa/variable/index create a table by filtering out outlier and single level variable

factor_vars_all<-vars_list$vars_track$factor_vars
plot_list<-list()
ommit_plots<-c()
statistics<-c("statistic", "p.value", "method","conf.low", "conf.high")
specifics<-c("variable", "taxa_level", "transformation","diversity_index", "var_levels", "possible_comparisons")
stats_all<-matrix(nrow = 1, ncol = length(c(statistics, specifics)))
colnames(stats_all)<-c(statistics, specifics)
for (tx in taxa_levels_select) {
  for (tr in select_transformations) {
    for (var in factor_vars_all) {
      for (i in div_index_select){
        data<-alpha_div_tbl %>% filter(taxa==tx, data_transformation==tr)
        #data[[i]]<-remove_outliers(data[[i]], 5, na.rm = TRUE)
        data<-filter(data, (!is.na(data[[i]])))
        levels<-data[[var]] %>% as.character() %>% as.factor() %>%levels()
        level_count=length(levels)
        comparison_count=level_count*(level_count-1)/2
        if (level_count>1) {      
          plot<-ggbetweenstats(
            data = data,
            x = !!var,
            y = !!i,
            type = "Bayes Factor",
            plot.type = "violin",
            #ggsignif.args = list(textsize = 1, tip_length = 0.01),
            pairwise.comparisons = TRUE,
            title = paste0("Taxa: ",tx, " ;Transformation: ",tr," ;Diversity index: ", i),
            p.adjust.method = "holm",
            outlier.color = "grey",
            caption = TRUE,
            results.subtitle = TRUE
          )
          p<-plot + 
            theme(plot.margin = margin(1, 1, 1, 1, "cm"))
          plot_list[[tr]][[var]][[tx]][[i]]<-p
          stat<-plot %>% extract_stats()
          stat<-as.data.frame(stat[["subtitle_data"]][,statistics])
          if (length(stat)>0) {
            stat[, specifics]=c(var, tx, tr,i, level_count, comparison_count)
            stats_all<-rbind(stats_all, stat)
          }

        }else{
          ommit_plots<-c(ommit_plots, var)
        }
      }
    }
  }
}


write.csv(stats_all, (paste("tables/divs_vs_VarsAll_stats_",project,"_",microbiome_seq_strategy,"_",method,".csv", sep = "")))
stats_all<- read.csv(paste("tables/divs_vs_VarsAll_stats_",project,"_",microbiome_seq_strategy,"_",method,".csv", sep = ""), row.names = 1)

###investigate the statistics table and choose the right variable to plot
###identify the real significant factors
###significant factors are the one with highest ratio of number of significant comparisons 
#divided by level of the variable (e.g. if var X has )
#miximum signicant comparison is 3
sig_vars_taxa<-data.frame()
max_var_levels=7
for (tx in taxa_levels_select) {
  for (tr in select_transformations) {
    ### to calculate which variable has more frequent significant comparison per taxa per transformation
    sig_vars<-stats_all %>% filter(var_levels<=max_var_levels,taxa_level==tx, transformation==tr) %>% 
      filter(p.value<0.05) %>%
      filter(conf.high*conf.low>0)
    sig_count<-sig_vars$variable %>% table() %>% as.data.frame()
    colnames(sig_count)<-c("variable", "sig_frequency_per_taxa")
    sig_vars<-sig_vars %>% left_join(sig_count, by="variable") %>%
      arrange(desc(sig_frequency_per_taxa)) %>% 
      select(c("taxa_level","transformation","variable", "var_levels",
               "sig_frequency_per_taxa", "possible_comparisons")) %>% 
      distinct(.keep_all=TRUE)
    sig_vars_taxa<-rbind(sig_vars_taxa, sig_vars)
  }
}

### to calculate which variable has more frequent significant comparison overal:
sig_vars_taxa$var_levels=as.numeric(sig_vars_taxa$var_levels)
sig_vars_taxa$possible_comparisons=as.numeric(sig_vars_taxa$possible_comparisons)
sig_vars_all<-sig_vars_taxa %>% mutate(sig_frq_perTaxa_nrmalizd=sig_frequency_per_taxa/possible_comparisons)
sig_vars_all<- data.table::data.table(sig_vars_all, key="variable")
sig_vars_all<-sig_vars_all[, list(sig_all=sum(sig_frq_perTaxa_nrmalizd)),by=variable]

### create plots of significant variables only
factor_vars_select<-sig_vars_all %>% filter(sig_all>0.5) %>% .$variable
plot_list<-list()
ommit_plots<-c()

statistics<-c("statistic", "p.value", "method","conf.low", "conf.high")
specifics<-c("variable", "taxa_level", "transformation","diversity_index", "var_levels", "possible_comparisons")
results<-matrix(nrow = 1, ncol = length(c(statistics, specifics)))
colnames(results)<-c(statistics, specifics)


for (tx in taxa_levels_select) {
  for (tr in select_transformations) {
    for (var in factor_vars_select) {
      for (i in div_index_select){
        data<-alpha_div_tbl %>% filter(taxa==tx, data_transformation==tr)
        #data[[i]]<-remove_outliers(data[[i]], 5, na.rm = TRUE)
        data<-filter(data, (!is.na(data[[i]])))
        levels<-data[[var]] %>% as.character() %>% as.factor() %>%levels()
        level_count=length(levels)
        comparison_count=level_count*(level_count-1)/2
        if (level_count>1) {      
          plot<-ggbetweenstats(
            data = data,
            x = !!var,
            y = !!i,
            type = "nonparametric",
            plot.type = "violin",
            #ggsignif.args = list(textsize = 1, tip_length = 0.01),
            pairwise.comparisons = TRUE,
            title = paste0("Taxa: ",tx, " ;Diversity index: ", i),
            p.adjust.method = "holm",
            outlier.color = "grey",
            caption = TRUE,
            results.subtitle = TRUE
          )
          p<-plot + 
            theme(plot.margin = margin(1, 1, 1, 1, "cm"))
          plot_list[[tr]][[var]][[tx]][[i]]<-p
          stat<-plot %>% extract_stats()
          stat<-as.data.frame(stat[["subtitle_data"]][,statistics])
          stat[, specifics]=c(var, tx, tr,i, level_count, comparison_count)
          stats_all<-rbind(stats_all, stat)
          
        }else{
          ommit_plots<-c(ommit_plots, var)
        }
      }
    }
  }
}

#plot all the plots in the list
#separate plot-sets for each variable
for (tr in select_transformations) {
  for (var in factor_vars_select){
    plots<-plot_list[[tr]][[var]]
    file_name = paste("figures/alpha_diversity/",project,"_",microbiome_seq_strategy,"_",method,"_",tr,"_",var,"_alpha_div.pdf", sep = "")
    pdf(width = 55, height =25)
    ggsave(file = file_name, arrangeGrob(grobs = unlist(plots, recursive = FALSE),ncol = length(div_index_select), 
                                         nrow =length(taxa_levels_select) ), limitsize = FALSE) 
    dev.off()
  } 
  
}





######### For the repeated measure #####
#### ############## ##############

alpha_div_tbl_long=read_csv(paste0("tables/alpha_diversities_all_samples_", project, "_", microbiome_seq_strategy,"_",method, ".csv"))

repeat_measure_vars=c("Treatment_Timepoints_week", "Placebo_Timepoints_week")
plot_list<-list()
ommit_plots<-c()

statistics<-c("p-value")
specifics<-c("variable", "levels","taxa_level", "transformation","diversity_index")
stats_all<-matrix(nrow = 1, ncol = length(c( specifics,statistics)))
colnames(stats_all)<-c(specifics,statistics)


# Visualize paired (repeated measure) comparisons
## statistics: Linear mixed effects model (LME)
library(nlme)
library(gridExtra)
for (tx in taxa_levels_select) {
  for (tr in select_transformations) {
    for (var in repeat_measure_vars) {
      for (i in div_index_select){
        data_long<-alpha_div_tbl_long %>% filter(taxa==tx, data_transformation==tr) %>% as.data.frame()
        data_long<-filter(data_long, (diversity_index_name==i))
        data_long$subject=data_long$Subject_ID %>% as.factor()## subject.id will be used by package to detect pairs
        data_long<-data_long %>% select (subject, all_of(var), diversity_value)
        data_long<-na.omit(data_long)
        data_long[, var]=as.numeric(data_long[,var])
        data_long[, var]=as.factor(data_long[, var])
        levels<-data_long[[var]] %>% as.factor() %>%levels()
        data_long<-data_long %>% dplyr::rename("repeated_variable"=all_of(var))
        formula<-paste0("diversity_value ~ ", var) %>% as.formula()
        mixed_effects_model <- lme(diversity_value ~ repeated_variable , random = ~ 1 | subject, data = data_long)
        glm_model <- glm(diversity_value ~ repeated_variable , data = data_long)
        data_long$fitted_values <- predict(mixed_effects_model)
        p_val_df=coef(summary(mixed_effects_model)) %>% as.data.frame() %>% select(`p-value`)
        p_val_df$`p-value`=p_val_df$`p-value` %>% round(digits = 2)
        p_val_df$levels=levels
        p_val_df<-p_val_df[, c("levels", "p-value")]
        rownames(p_val_df)=levels
        stats=p_val_df
        stats$taxa_level=tx
        stats$transformation=tr
        stats$variable=var
        #stats$level=levels
        stats$diversity_index=i
        stats_all=rbind(stats_all, stats)
        
        # Plotting
        plot<-ggplot(data_long, aes(x = repeated_variable, y =diversity_value, 
                                    group = subject, color = factor(subject))) +
          geom_line(show.legend = FALSE) +
          geom_point(show.legend = FALSE) +
          labs(title = paste0(" Transformation: ", tr,"; Taxa Level: ",tx, "\n",
                              " Paired comparison: ", var, "\n",
                              " Diversity index: ", i, "\n",
                              " p.values: ", toString(p_val_df$`p-value`)),
               x = var,
               y = i,
               color = "Subject") +
          annotation_custom(tableGrob(p_val_df, rows=NULL),xmin=5, xmax = 9)

        
        plot
        
          p<-plot + 
          theme(plot.margin = margin(1, 6, 1, 1, "cm"))
          plot_list[[tr]][[var]][[tx]][[i]]<-p
        }
      }
    }
  }

#plot all the plots in the list
#separate plot-sets for each variable

for (tr in select_transformations) {
  for (var in repeat_measure_vars){
    plots<-plot_list[[tr]][[var]]
    file_name = paste("figures/alpha_diversity/",project,"_",microbiome_seq_strategy,"_",method,"_",tr,"_",var,"_repeated_measure_alpha_div.pdf", sep = "")
    pdf(width = 45, height =25)
    ggsave(file = file_name, arrangeGrob(grobs = unlist(plots, recursive = FALSE),ncol = length(div_index_select), 
                                         nrow =length(taxa_levels_select) ), limitsize = FALSE) 
    dev.off()
  } 
  
}




#### scatterplot regression####
select_transformations=c("original", "cmp", "hlg")
vars_list<-readRDS(file = paste("objects/vars_list_",project,".rds", sep = ""))
continuous_vars<-vars_list$vars_track$continues_vars
plot_list<-list()
ommit_plots<-c()

statistics<-c("statistic", "p.value", "method","conf.low", "conf.high")
specifics<-c("variable", "taxa_level", "transformation","diversity_index")
results<-matrix(nrow = 1, ncol = length(c(statistics, specifics)))
colnames(results)<-c(statistics, specifics)
stats_all<-matrix(nrow = 1, ncol = length(c(statistics, specifics)))
colnames(stats_all)<-c(statistics, specifics)

for (tx in taxa_levels_select) {
  for (tr in select_transformations) {
    for (var in continuous_vars) {
      for (i in div_index_select){
        data<-alpha_div_tbl %>% filter(taxa==tx, data_transformation==tr)
        #data[[i]]<-remove_outliers(data[[i]], 5, na.rm = TRUE)
        data<-filter(data, (!is.na(data[[i]])))
          plot<-ggscatterstats(
            data = data,
            x = !!var,
            y = !!i,
            type = "parametric",
            title = paste0("Taxa: ",tx, " ;Diversity index: ", i),
            caption = TRUE,
            subtitle = TRUE
          )
          p<-plot + 
            theme(plot.margin = margin(1, 1, 1, 1, "cm"))
          plot_list[[tr]][[var]][[tx]][[i]]<-p
          stat<-plot %>% extract_stats()
          stat<-as.data.frame(stat[["subtitle_data"]][,statistics])
          stat[, specifics]=c(var, tx, tr,i)
          stats_all<-rbind(stats_all, stat)
        }
      }
    }
  }

#plot all the plots in the list
#separate plot-sets for each variable
for (tr in select_transformations) {
  for (var in continuous_vars){
    plots<-plot_list[[tr]][[var]]
    file_name = paste("figures/alpha_diversity/",project,"_",microbiome_seq_strategy,"_",method,"_",tr,"_",var,"_correlation_scatterPlot.pdf", sep = "")
    pdf(width = 40, height =25)
    ggsave(file = file_name, arrangeGrob(grobs = unlist(plots, recursive = FALSE),ncol = length(div_index_select), 
                                         nrow =length(taxa_levels_select) ), limitsize = FALSE) 
    dev.off()
  } 
  
}


