library(phyloseq)
library(microbiome)
library(DESeq2)
library(vegan)
library(ggplot2)

# this script is for various transformation of OTU abundance data

setwd("/Users/hhakimjavadi/Dropbox (UFL)/FMT_microbiome_smichail/")
project="FMT"
microbiome_seq_strategy="WGS"
method="Bracken"

psq_orig=readRDS(paste("objects/",project, "_psq_orig_",microbiome_seq_strategy,"_",method,".RDS", sep = ""))

#remove samples with total reads lower than 500 (we did not use it)
which(colSums(otu_table(psq_orig)) < 1000)
psq_high_abund <- prune_samples(colSums(otu_table(psq_orig)) > 1000, psq_orig)
### counts in origanal form (Poore 2020) contain negative values that does not allow several transformation

#power-transform abundance 

rrfd = rarefy_even_depth(psq_orig, sample.size = 1000,rngseed = 711, replace = FALSE)
##first transform original psq, then aggregate (the other way (aggregate > transform) does not work late)
psq_transformed<-list(
  original=psq_orig,
  #log = transform_sample_counts(psq_orig, function(x) log(1+x)),
  #normalizations
  hlg = microbiome::transform(psq_orig, "hellinger"),
  #idnt = microbiome::transform(psq_orig, "identity"),
  #clr = microbiome::transform(exp, "clr"),
  cmp = microbiome::transform(psq_orig, "compositional"), #same as relative
  #to rerefy, first find a right sample.size!
  rrfd = rrfd
)

# create multiple phyloseq objects at the all taxonomic levels by aggregating 
#transformed psq based on taxa levels
psq_transforms<-names(psq_transformed)
taxa_levels<-c("superkingdom","clade","phylum", "class","order","family","genus", "species") 
taxa_levels<-c("phylum","order","genus", "species") 


psq_list<-list()

#to store psq object at various taxa level into a list:
for (tr in psq_transforms) {
  for (tx in taxa_levels) {
    psq_list[[tr]][[tx]]=
      microbiome::aggregate_taxa(psq_transformed[[tr]], level =tx)
  }
}

## if tree is available (very slow)
for (tr in psq_transforms) {
  for (tx in taxa_levels) {
    psq_list[[tr]][[tx]]=
      tip_glom(psq_transformed[[tr]])
  }
}

saveRDS(psq_list, file=paste("objects/",project, "_psq_list_",microbiome_seq_strategy,"_",method,".rds", sep = ""))



### to update metadata for all tranformed psq objects
## because tranformation is lengthy process
psq_list<-readRDS(file=paste("objects/",project, "_psq_list_",microbiome_seq_strategy,"_",method,".rds", sep = ""))
metadata<-read.csv(paste("data/metadata/metadata_",project,"_",microbiome_seq_strategy,"_processed",".csv", sep = ""))
rownames(metadata)<-metadata[, "Sample_ID"]

psq_transforms<-names(psq_list)
taxa_levels<-c("phylum","order","genus", "species")

#to store psq object at various taxa level into a list:
for (i in 1:length(psq_transforms)) {
  for (j in 1:length(taxa_levels)) {
    sample_data(psq_list[[psq_transforms[i]]][[taxa_levels[j]]])=metadata
  }
}

saveRDS(psq_list, file=paste("objects/",project, "_psq_list_",microbiome_seq_strategy,"_",method,".rds", sep = ""))


