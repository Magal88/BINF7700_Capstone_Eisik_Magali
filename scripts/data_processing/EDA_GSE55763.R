#“Exploratory Data Analysis GSE55763”

#Load require packages

library(dplyr)
library(ggplot2)
library(gridExtra)
library(kableExtra)
library(knitr)
library(htmltools)
library(factoextra)
library(FactoMineR)
library(tibble)

#Load data set
GSE55763_beta <- read.csv("top500_CpGs_age.csv", header = TRUE, stringsAsFactors = FALSE)

#Missing data Analysis

#Missing data quantity
total_na <- sum(is.na(GSE55763_beta))
cat("Total number of missing values:", total_na, "\n")
## Total number of missing values: 992
#  Percentage of missing values
percent_na <- mean(is.na(GSE55763_beta)) * 100
cat(" Percentage of missing values:", round(percent_na, 2), "%\n")
##  Percentage of missing values: 0.07 %
#In this case, missing values are assumed to be Missing Completely at Random (MCAR) 
#and will be imputed using the median of each CpG.

X_cpg <- GSE55763_beta[, 4:503]  #select cpg columns
# Impute missing values using the median 
X_cpg_imputed <- as.data.frame(lapply(X_cpg, function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  return(x)
}))
GSE55763_beta <- cbind(GSE55763_beta[, 1:3], X_cpg_imputed) #combine with other variables
sum(is.na(GSE55763_beta)) #verify NAs

#Save data set in csv file

write.csv(GSE55763_beta, "GSE55763_beta_500_imputed.csv", row.names = FALSE)
#Visualizations

#Histogram for age

#Convert variables for visualization and analysis
GSE55763_beta $age <- as.numeric(GSE55763_beta $age)
GSE55763_beta $gender <- factor(GSE55763_beta $gender, levels = c("M", "F"))
ggplot(GSE55763_beta, aes(x = age)) +
  geom_histogram(binwidth = 5, fill = "violet", color = "black") +
  theme_minimal() +
  labs(title = "Age Distribution", x = "Age", y = "Samples")
#Barplot for gender

ggplot(GSE55763_beta, aes(x = gender, fill = gender)) +
  geom_bar() +
  scale_fill_manual(values = c("steelblue", "pink")) +
  theme_minimal() +
  labs(title = "Gender Distribution", x = "Gender", y = "Samples")


# Count of each gender
gender_count <- table(GSE55763_beta$gender)

# Convert to data.frame
gender_df <- as.data.frame(gender_count)
names(gender_df) <- c("Gender", "Count")

# Calculate percentage
gender_df$Percent <- round((gender_df$Count / sum(gender_df$Count)) * 100, 2)
gender_df %>%
  kbl(
    caption = tags$b("Gender distribution in GSE55763")
  ) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed"),
    full_width = FALSE
  ) %>%
  column_spec(3, bold = TRUE)

#Gender distribution in GSE55763


ggplot(GSE55763_beta, aes(x = gender, y = age, fill = gender)) +
  geom_violin(alpha = 0.6) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  scale_fill_manual(values = c("green2", "coral1")) +
  theme_minimal() +
  labs(title = "Age distribution by Gender",
       x = "Gender",
       y = "Age") +
  theme(legend.position = "none")


#Boxplot

ggplot(GSE55763_beta, aes(x = gender, y = age, fill = gender)) +
  geom_boxplot(alpha = 0.6, outlier.size = 1) +
  scale_fill_manual(values = c("steelblue", "magenta")) +
  theme_minimal() +
  labs(title = "Age distribution by Gender",
       x = "Gender",
       y = "Age") +
  theme(legend.position = "none")


#PCA

# Scale numeric features (top 50 CpGs)
numeric_data <- GSE55763_beta %>%
  select(where(is.numeric))

# Select top 50 CpG columns 
X_cpg_top50 <- numeric_data[, 1:50]
X_cpg_scaled <- scale(X_cpg_top50)

# Run PCA 
pca_res <- PCA(X_cpg_scaled, graph = FALSE)

# PCA plot for individuals (samples)
fviz_pca_ind(
  pca_res,
  geom.ind = "point",
  pointsize = 2,
  alpha.ind = 0.6,                        
  col.ind = GSE55763_beta$gender, 
  palette = c("steelblue", "tomato"),
  addEllipses = TRUE,                    
  legend.title = "Gender"
) +
  theme_minimal() +
  ggtitle("PCA of Top 50 CpGs (colored by Gender)")


fviz_pca_var(
  pca_res,
  col.var = "contrib",                     
  select.var = list(contribution = 20)     
) +
  theme_minimal() +
  ggtitle("Top 10 CpGs Contributing to PCA")


cpg_data <- GSE55763_beta %>%
  select(starts_with("cg"))

# Scale
X_cpg_top20 <- scale(cpg_data[, 1:20])

# PCA
pca_res <- PCA(X_cpg_top20, graph = FALSE)

# Contribuciones de variables
var_contribution <- get_pca_var(pca_res)$contrib
var_contrib_df <- as.data.frame(var_contribution) %>%
  rownames_to_column(var = "CpG") %>%
  select(CpG, Dim.1, Dim.2) %>%        # Solo PC1 y PC2
  rename(PC1_Contribution = Dim.1,
         PC2_Contribution = Dim.2) %>%
  mutate(Total_Contribution = PC1_Contribution + PC2_Contribution)

# Top 20
top20_cpgs <- var_contrib_df %>%
  arrange(desc(Total_Contribution)) %>%
  slice(1:20)

# Display table
top20_cpgs %>%
  kbl(caption = "Top 20 CpGs contributing to PCA") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE) %>%
  column_spec(2:4)

g1<- ggplot(GSE55763_beta, aes(x = gender, y = age, fill = gender)) +
  geom_boxplot(alpha = 0.6, outlier.size = 1) +
  scale_fill_manual(values = c("steelblue", "magenta")) +
  theme_minimal() +
  labs(title = "Age distribution by Gender",
       x = "Gender",
       y = "Age") +
  theme(legend.position = "none")
ggsave(filename = "boxplot.png", plot = g1, width = 6, height = 4, dpi = 300)
g2<- ggplot(GSE55763_beta, aes(x = gender, y = age, fill = gender)) +
  geom_violin(alpha = 0.6) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  scale_fill_manual(values = c("green2", "coral1")) +
  theme_minimal() +
  labs(title = "Age distribution by Gender",
       x = "Gender",
       y = "Age") +
  theme(legend.position = "none")
ggsave(filename = "violin_plot.png", plot = g2, width = 6, height = 4, dpi = 300)
g3<- ggplot(GSE55763_beta, aes(x = gender, fill = gender)) +
  geom_bar() +
  scale_fill_manual(values = c("steelblue", "pink")) +
  theme_minimal() +
  labs(title = "Gender Distribution", x = "Gender", y = "Samples")
ggsave(filename = "barplot.png", plot = g3, width = 6, height = 4, dpi = 300)
g4<- ggplot(GSE55763_beta, aes(x = age)) +
  geom_histogram(binwidth = 5, fill = "violet", color = "black") +
  theme_minimal() +
  labs(title = "Age Distribution", x = "Age", y = "Samples")
ggsave(filename = "histogram.png", plot = g4, width = 6, height = 4, dpi = 300)