#"M-values"

  
  
#In this document, all CpG Beta values were converted into M values using a logit-like transformation,
#log2(beta / (1 - beta)), to stabilize variance and improve the suitability of the data for statistical modeling.


#Load require packages
library(dplyr)
library(ggplot2)
library(kableExtra)



#Load data set
GSE55763_imputed <- read.csv("GSE55763_beta_500_imputed.csv", 
                             header = TRUE, 
                             stringsAsFactors = FALSE)

head(GSE55763_imputed,3)




# Identify CpG columns
cpg_columns <- grep("^cg", colnames(GSE55763_imputed), value = TRUE)

# Copy the dataset in object called GSE55763_M
GSE55763_M <- GSE55763_imputed

# Extract Beta values
beta_matrix <- GSE55763_M[, cpg_columns]

beta_matrix[beta_matrix == 0] <- 1e-6
beta_matrix[beta_matrix == 1] <- 1 - 1e-6

# Transform Beta to M values (logit transformation)
M_matrix <- log2(beta_matrix / (1 - beta_matrix))

# Replace CpG columns with M values
GSE55763_M[, cpg_columns] <- M_matrix



#Verify first 3 rows
head(GSE55763_M,3)



# Save file in CSV format
write.csv(GSE55763_M, "GSE55763_M_values.csv", row.names = FALSE)
