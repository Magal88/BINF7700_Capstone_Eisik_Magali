
#"Select cpgs most correlated with age"

  
  

#Load require packages
library(dplyr)    
library(readr)    
library(knitr)    
library(kableExtra)
library(knitr)



#Load data set
merged_500 <- read_csv(
  "GSE55763_beta_values_merged.csv",
  show_col_types = FALSE
)



#Metadata columns

meta_cols <- c("X", "gender", "age")

#Extract only CpG columns

cpg_matrix <- merged_500 %>% select(-all_of(meta_cols))

#Ensure numeric (for cpgs)

cpg_matrix <- cpg_matrix %>% mutate(across(everything(), as.numeric))


#Compute correlation and p-value with age for each CpG (Pearson correlation coefficient)

cors_pvals <- sapply(cpg_matrix, function(x) {
  test <- cor.test(x, merged_500$age)
  c(cor = test$estimate, pval = test$p.value)
})



#Convert to dataframe

cors_df <- as.data.frame(t(cors_pvals))

#Filter CpGs with p < 0.01

cors_sig <- cors_df[cors_df$pval < 0.01, ]



#Select top 500 CpGs by absolute correlation among significant CpGs

top500 <- rownames(cors_sig[order(abs(cors_sig$cor), decreasing = TRUE), ])[1:500]

subset_500 <- merged_500 %>% select(all_of(meta_cols), all_of(top500))


head(subset_500)


# Save the top 500 CpGs with metadata to a CSV
write_csv(subset_500, "top500_CpGs_age.csv")
