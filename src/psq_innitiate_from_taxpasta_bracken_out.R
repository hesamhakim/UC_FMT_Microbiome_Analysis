
library(phyloseq)
library(tidyverse)
setwd("/Users/hhakimjavadi/Dropbox (UFL)/FMT_microbiome_smichail/")
project="FMT"
method="Bracken"
microbiome_seq_strategy="WGS"
### This script takes the bracken output of the taxprofiler/taxpasta and genrates phyloseq object
taxpasta_bracken_abund<-read_tsv("workflow/taxaprofiler/results_all/taxpasta/kraken2_std_bracken-bracken_name_clean.txt")


########1:taxa table ##############
######## using the same bracken table we can create a taxa table to be used as taxa annotation in the psq object
taxa_ids=taxpasta_bracken_abund[c("taxonomy_id","name")]
taxa_ids$taxonomy_id=as.character(taxa_ids$taxonomy_id)
taxa_rnks=taxpasta_bracken_abund %>% 
  select("taxonomy_id","name", "rank" , "id_lineage","rank_lineage") %>% filter(rank=="species")

###to create a long table of taxa, rank, rank name and rank lineage
taxa_tbl=data_frame(taxonomy_id=character(),id_lineage=character(), rank_lineage=character())
for (i in taxa_rnks$taxonomy_id) {
  table=filter(taxa_rnks, taxonomy_id==i)
  id_lineages=table$id_lineage %>% str_split_1(pattern = ";")
  rank_lineages=table$rank_lineage %>% str_split_1(pattern = ";")
  taxa_tbl_tmp=data_frame(taxonomy_id=table$taxonomy_id,
                       id_lineage=id_lineages, rank_lineage=rank_lineages)
  taxa_tbl=rbind(taxa_tbl, taxa_tbl_tmp)
}

## add taxa names to ids
taxa_tbl<-left_join(taxa_tbl, taxa_ids, by=c("id_lineage"="taxonomy_id"))
taxa_tbl<-filter(taxa_tbl, !is.na(name))

taxa_tbl_wide=pivot_wider(taxa_tbl, id_cols =taxonomy_id, names_from = rank_lineage, values_from = name )
taxa_tbl_wide=taxa_tbl_wide %>% select(phylum, class,order,family,genus,species)

rownames(taxa_tbl_wide)=taxa_tbl_wide$species


write.csv(taxa_tbl_wide, "workflow/taxaprofiler/results_all/taxpasta/ps_taxa_table.csv", row.names = TRUE)
taxa_tbl_wide=read.csv("workflow/taxaprofiler/results_all/taxpasta/ps_taxa_table.csv", row.names = 1)
taxa_tbl_mt<-as.matrix(taxa_tbl_wide)

#### 2: abundance table ##########
abund_tbl=taxpasta_bracken_abund %>% filter(rank=="species") %>%
  select(-taxonomy_id,	-rank,	-id_lineage,	-rank_lineage) %>% rename(taxa_name=name)
## to replace file name with sample Id
abund_tbl_long=pivot_longer(abund_tbl, cols = ends_with(".gz"), values_to = "abundance", names_to = "file_name" )
file_sample_map=read_delim("data/metadata/file_name_sample_map.txt")
abund_tbl_long=left_join(abund_tbl_long, file_sample_map, by=c("file_name"="fastq_file_name"))
abund_tbl_long=abund_tbl_long %>% select(-file_name)
abund_tbl=abund_tbl_long %>% pivot_wider(names_from = Sample_ID, values_from = abundance )
abund_tbl=abund_tbl %>% column_to_rownames("taxa_name")

####3: metadata (sample) table
metadata <- read_csv("data/metadata/metadata_FMT_WGS_processed.csv")
metadata<-column_to_rownames(metadata, "Sample_ID") 

### I generated a general tree including all Bacteria data from NCBI taxonomy 
## check  the python script phylo_tree_from_ncbi_taxonomy_ete3.ipynb
tree <- ape::read.tree("/Users/hhakimjavadi/UFL Dropbox/Hesamedin Hakimjavadi/FMT_microbiome_smichail/objects/ncbi_taxonomy/ncbi_tree.nwk")


abund = otu_table(abund_tbl, taxa_are_rows = TRUE)
TAX = tax_table(taxa_tbl_mt)
sample = sample_data(metadata)
Tree=phy_tree(tree)
psq_orig <- phyloseq(abund,TAX,sample)



#plot_tree(psq_orig, label.tips = "taxa_names", ladderize = "left", color = "Phylum")

saveRDS(psq_orig, paste("objects/",project, "_psq_orig_",microbiome_seq_strategy,"_",method,".RDS", sep = ""))






ncbi_taxa_number <- taxa_ids[,"taxonomy_id"] %>% distinct() 

write_tsv(ncbi_taxa_number, "/Users/hhakimjavadi/UFL Dropbox/Hesamedin Hakimjavadi/FMT_microbiome_smichail/workflow/taxaprofiler/results_all/taxpasta/ncbi_taxa_number_all.txt")

