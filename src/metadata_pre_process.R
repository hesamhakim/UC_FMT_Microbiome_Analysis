library(SummarizedExperiment)
library(UCSCXenaTools)
library(dplyr)
library(tidyverse)
library(readxl)
#this script download and merges 4 metadata sources:
# metadata from xena and cBioportal, metadata & biospecimen from gdc, and metadat from Poore et.al. paper: 
setwd("/Users/hhakimjavadi/Dropbox (UFL)/FMT_microbiome_smichail/")
project="FMT"
microbiome_seq_strategy="WGS"
metadata<-read_excel("data/metadata/FMT_Metadata_08.17.24_pucai_responder_updates.xlsx")
metadata<-as.data.frame(metadata)
colnames(metadata)=gsub("-", "_", colnames(metadata))

## ad fastq file names to clinical data
file_name_sample_map <- read_delim("data/metadata/file_name_sample_map.txt")
metadata<-left_join(metadata, file_name_sample_map, by="Sample_ID")
#Create Categorical Variable from Continuous 
metadata$age <-as.numeric(metadata$age)
metadata$age_group<-metadata$age %>% 
  cut(breaks =c(0,10,15,20,30,40,50),
      labels = c("1_9", "10_14", "15_19", "20_29", "30_39", "40_50"))

metadata$pucai_range<-metadata$pucai %>% 
  cut(breaks =c(-1,10,35,65, 100),
      labels = c("0_9", "10_34", "35_64", "65+"))

metadata$pucai_base_uc_category<-metadata$pucai %>% 
  cut(breaks =c(-1,10,35,65, 100),
      labels = c("Remission", "Mild", "Moderate", "Severe"))

metadata$bmi_category<-metadata$bmi %>% 
  cut(breaks =c(0,18.59,24.99,29.99,34.99, 39.99, 44.99),
      labels = c("<18.5", "18.6_24.9", "25_29.9", "30_34.9", "35_39.9", "40_44.9"))

metadata$Treatment_Timepoints_week=as.numeric(metadata$Treatment_Timepoints_week)
metadata$Placebo_Timepoints_week=as.numeric(metadata$Placebo_Timepoints_week)

vars_levels<-list()
for (i in 1:length(colnames(metadata))) {
  name=colnames(metadata)[i]
  levels_count= length(levels(as.factor(metadata[,i])))  
  levels_names= levels(as.factor(metadata[,i]))
  class=class(metadata[,i])
  temp=list(var_name=name, levels_count=levels_count, levels_names=levels_names, class=class)
  vars_levels[[name]]<-temp
}

#separate continues and categorical variables
require(rlist)
continues_vars_levels<-vars_levels %>% subset((class=="numeric" | class=="integer")& levels_count>3)
factor_vars_levels1<-vars_levels %>% subset((class=="character" |class=="factor") & levels_count>1 & levels_count<=8)
factor_vars_levels2<-vars_levels %>% subset((class=="numeric" | class=="integer") & levels_count>1 & levels_count<=8)
factor_vars_levels<-c(factor_vars_levels1, factor_vars_levels2)
single_level_vars<-vars_levels %>% subset(levels_count==1) %>% names()
two_levels_vars<-vars_levels %>% subset(levels_count==2) %>% names()
three_levels_vars<-vars_levels %>% subset(levels_count==3) %>% names()
four_levels_vars<-vars_levels %>% subset(levels_count==4) %>% names()
five_levels_vars<-vars_levels %>% subset(levels_count==5) %>% names()
two_five_levels_vars<-vars_levels %>% subset(levels_count>=2 & levels_count<=5) %>% names()
two_four_levels_vars<-vars_levels %>% subset(levels_count>=2 & levels_count<=4) %>% names()
three_eight_levels_vars<-vars_levels %>% subset(levels_count>=3 & levels_count<=8) %>% names()

continues_vars<-names(continues_vars_levels)
factor_vars<-names(factor_vars_levels)

rm(factor_vars_levels1, factor_vars_levels2)
length(intersect(continues_vars,factor_vars))==0 #make sure there are no common elements between continues and factor vars
#if want to remove single level vars:
metadata<-dplyr::select(metadata, -c(which(colnames(metadata) %in% single_level_vars)))

#to convert all factor vars to factor class (factor_vars are still character):
#metadata[sapply(metadata, is.character)] <- lapply(metadata[sapply(metadata, is.character)], as.factor)
metadata[,factor_vars]<-lapply(metadata[,factor_vars], as.factor)

## ad fastq file names to clinical data
file_name_sample_map <- read_delim("data/metadata/file_name_sample_map.txt")


#store metadata and use it in the next step scripts
write.csv(metadata, paste("data/metadata/metadata_",project,"_",microbiome_seq_strategy,"_processed",".csv", sep = ""), row.names = FALSE)


#list all vars required for next step all variables
#create list step-wise so that a non-exist var cannot disrupt creation of the list
vars_track<-list()
vars_track[["pucai"]]=metadata$pucai
vars_track[["pucai_range"]]=metadata$pucai_range
vars_track[["pucai_base_uc_category"]]=metadata$pucai_base_uc_category
vars_track[["Group"]]=metadata$Group
vars_track[["responder"]]=metadata$responder
vars_track[["age_group"]]=metadata$age_group
vars_track[["gender"]]=metadata$gender
vars_track[["race"]]=metadata$race
vars_track[["bmi_category"]]=metadata$bmi_category
vars_track[["Open_label_follow_up"]]=vars_track[["Open_label_follow_up"]]=metadata$Open_label_follow_up
vars_track[["Treatment_Placebo_timepoint"]]=metadata$Treatment_Placebo_timepoint
vars_track[["Placebo_Timepoints"]]=metadata$Placebo_Timepoints


vars_track[["vars_all"]]=c(factor_vars, continues_vars) %>% intersect(colnames(metadata))
vars_track[["factor_vars"]]=factor_vars%>% intersect(colnames(metadata))
vars_track[["continues_vars"]]=continues_vars%>% intersect(colnames(metadata))
vars_track[["single_level_vars"]]=single_level_vars%>% intersect(colnames(metadata))
vars_track[["two_levels_vars"]]=two_levels_vars%>% intersect(colnames(metadata))
vars_track[["three_levels_vars"]]=three_levels_vars%>% intersect(colnames(metadata))
vars_track[["four_levels_vars"]]=four_levels_vars%>% intersect(colnames(metadata))
vars_track[["two_five_levels_vars"]]=two_five_levels_vars%>% intersect(colnames(metadata))
vars_track[["two_four_levels_vars"]]=two_four_levels_vars%>% intersect(colnames(metadata))
vars_track[["three_eight_levels_vars"]]=three_eight_levels_vars%>% intersect(colnames(metadata))


vars_list<-list(vars_levels=vars_levels, vars_track=vars_track)
saveRDS(vars_list, file = paste("objects/vars_list_",project,".rds", sep = ""))

#1-metadata_pre_process
#2_phyloseq_obj_innitiate
#3_filtering_transforming
#4_sample_clustering (re-run from metadata_pre... follwing new clustering)
#5_categorical_vars_dependency
#6_diversity_all_vars



