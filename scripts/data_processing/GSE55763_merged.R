
# Preparing GSE55763 Beta Value Datasets for Modeling 

  
#This script processes the GSE55763 methylation dataset to create a tidy dataset with samples (GSM IDs) 
#as rows, CpG probes as columns, and associated metadata including sex and age. 
#The resulting dataset is ready for downstream modeling and analysis.


#Load require packages
library(tidyverse)
library(dplyr)
library(tibble)
library(tidyr)
library(readr)


#Load metadata


metadata_GSE55763 <- read.csv("GSE55763_metadata_clean.csv", header = TRUE, stringsAsFactors = FALSE)

# Display first 3 rows



metadata_GSE55763 <- metadata_GSE55763 %>%
  filter(grepl("population study", dataset))
dim(metadata_GSE55763)



#This dataset contains the mapping between array IDs and their corresponding GSM sample IDs.


array_to_GSM_mapping <- read.csv("array_to_GSM_mapping_clean.csv", header = TRUE, stringsAsFactors = FALSE)
head(array_to_GSM_mapping)



# Keep only samples that are also present in the metadata
array_to_GSM_mapping <- array_to_GSM_mapping %>%
  filter(GSM_ID %in% metadata_GSE55763$X)

dim(array_to_GSM_mapping)


#The GSE55763_top1000_all_samples.csv dataset contains beta values for DNA methylation 
#at the 1000 most variable CpG sites across all available samples. 


betas_sub1000 <- read.csv("GSE55763_top1000_all_samples.csv", header = TRUE, stringsAsFactors = FALSE)
head(betas_sub1000 )
dim(betas_sub1000)


# Remove 'X' prefix from column names
colnames(betas_sub1000) <- sub("^X", "", colnames(betas_sub1000))

# Keep only columns that exist in the mapping
array_cols <- intersect(colnames(betas_sub1000), array_to_GSM_mapping$Array_ID)
length(array_cols)  

# Filter dataset to keep only ID_REF + array columns
betas_sub1000_filtered <- betas_sub1000 %>%
  select(ID_REF, all_of(array_cols))
```

betas_long <- betas_sub1000_filtered %>%
  pivot_longer(
    cols = -ID_REF,
    names_to = "Array_ID",
    values_to = "beta"
  ) %>%
  inner_join(array_to_GSM_mapping, by = "Array_ID")  # add GSM_ID



#Map Array_ID to GSM_ID



betas_wide <- betas_long %>%
  select(GSM_ID, ID_REF, beta) %>%
  pivot_wider(
    id_cols = GSM_ID,
    names_from = ID_REF,
    values_from = beta
  )

dim(betas_wide)  
head(betas_wide[,1:10])


# Beta values were reshaped and merged with sample metadata using GSM_ID.
# The final dataset was saved as a CSV for further analysis.


betas_wide <- betas_long %>%
  select(-Array_ID) %>%
  pivot_wider(
    names_from = ID_REF,
    values_from = beta
  )


#Merge with metadata


metadata_filtered <- metadata_GSE55763 %>%
  filter(X %in% betas_wide$GSM_ID)

merged_data <- metadata_filtered %>%
  inner_join(betas_wide, by = c("X" = "GSM_ID"))
merged_data <- merged_data %>%
  select(-tissue, -dataset)

dim(merged_data)  
head(merged_data,3)


#Save data set in CSV format


write.csv(merged_data, "GSE55763_merged1000.csv", row.names = FALSE)