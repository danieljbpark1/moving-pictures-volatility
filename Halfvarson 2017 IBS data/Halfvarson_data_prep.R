## OTUs 
otus <- read.delim("./data_application_code/Halfvarson_2017_BIOM/43623/otu_table.txt")
otus[1:5,1:5]

colnames(otus) <- gsub("X1629.", "", colnames(otus))
rownames(otus) <- otus$OTU_ID
otus$OTU_ID <- NULL

otus[1:5,1:5]
dim(otus)

write.table(otus, file = "./data_application_code/Halfvarson_2017_OtuTable.txt", 
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)


## Sample info 
mapping <- read.delim("./data_application_code/Halfvarson_2017_mapping/43623_mapping_file.txt")
dim(mapping)
colnames(mapping)


## Metadata
metadata <- read.delim("./data_application_code/Halfvarson_2017_SampleMetadata.txt")
metadata = metadata %>% 
  select(sample_name, bmi, calprotectin, cd_behavior, cd_location, cd_resection, 
         collection_timestamp, description, diagnosis_full, host_subject_id, 
         host_taxid, ibd_subtype, patientnumber, perianal_disease, timepoint, 
         sex, uc_extent, year_diagnosed)
metadata$sample_name <- gsub("1629.", "", metadata$sample_name)

head(metadata)
all(colnames(otus) %in% metadata$sample_name)

write.table(metadata, file = "./data_application_code/Halfvarson_2017_Metadata.txt", 
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)



## Tree 
library(ape)
hvs.tree <- read.tree("./data_application_code/Halfvarson_2017_BIOM/60941/insertion_tree.relabelled.tre")
hvs.tree
all(rownames(otus) %in% hvs.tree$tip.label)

save(hvs.tree, file = "./data_application_code/Halfvarson_2017_Tree.Rda") 
