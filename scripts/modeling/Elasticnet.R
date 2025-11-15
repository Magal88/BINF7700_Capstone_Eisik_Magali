 "Elastic Net"

  
#In this notebook, an Elastic Net regression model is implemented to build an epigenetic clock for age prediction. 
 #The model includes gender as a covariate to account for sex-specific differences in DNA methylation.

#Load require packages
library(dplyr)    
library(readr) 
library(tibble)
library(ggplot2)  
library(caret)
library(knitr)    
library(glmnet)
library(kableExtra)
library(knitr)
library(gridExtra)



#Load data set
subset_500 <- read.csv("GSE55763_M_values.csv", header = TRUE, stringsAsFactors = FALSE)


## Prepare matrix for Elastic Net model


set.seed(500)  # for reproducibility

# Metadata columns
meta_cols <- c("X", "gender", "age")

# CpG matrix (select cpgs)
X_cpg <- dplyr::select(subset_500, -all_of(meta_cols))
# Scale CpGs: mean 0, SD 1
X_scaled <- scale(X_cpg)

# Incorporate the variable Gender (factor) as covariate in the model
gender <- as.factor(subset_500$gender)
# Encode gender as dummy variable
X_gender <- model.matrix(~ gender - 1)  

# Predictors
X_pred <- cbind(X_scaled, X_gender)

# Response variable (outcome=age)
y_age <- subset_500$age



## Splitting Data into Training and Test Sets


# Split Data: Training (80%) and Test (20%)
train_index <- createDataPartition(y_age, p = 0.8, list = FALSE) 
X_train <- X_pred[train_index, ] 
X_test <- X_pred[-train_index, ] 
y_train <- y_age[train_index] 
y_test <- y_age[-train_index]


## Cross validation


# 5-fold cross-validation (training data)
cv_enet <- cv.glmnet(X_train, y_train, alpha = 0.5, nfolds = 5)
# Select optimal lambda
best_lambda <- cv_enet$lambda.min
cat("Optimal lambda:", best_lambda, "\n")
# Train final model with best lambda
enet_model <- glmnet(X_train, y_train, alpha = 0.5, lambda = best_lambda)


# Cross-validation results-plot
plot(cv_enet)
abline(v = log(best_lambda), col = "red", lty = 2)


## Model predictions

# Predictions on training set
y_train_pred <- predict(enet_model, newx = X_train, s = best_lambda)
# Predictions on test set
y_test_pred <- predict(enet_model, newx = X_test, s = best_lambda)

## Metrics

# Calculate training metrics
train_mae <- mean(abs(y_train - y_train_pred))
train_rmse <- sqrt(mean((y_train - y_train_pred)^2))
train_r2 <- cor(y_train, y_train_pred)^2
# Calculate test metrics
test_mae <- mean(abs(y_test - y_test_pred))
test_rmse <- sqrt(mean((y_test - y_test_pred)^2))
test_r2 <- cor(y_test, y_test_pred)^2
# Display training results
cat("=== Training Set ===\n")
cat("MAE:", round(train_mae, 2), "years\n")
cat("RMSE:", round(train_rmse, 2), "years\n")
cat("R²:", round(train_r2, 2), "\n\n")
# Display test results (FINAL METRICS)
cat("=== Test Set (FINAL METRICS) ===\n")
cat("MAE:", round(test_mae, 2), "years\n")
cat("RMSE:", round(test_rmse, 2), "years\n")
cat("R²:", round(test_r2, 2), "\n")


# Create comparison table
metrics_df <- data.frame(
  Metric = c("MAE (years)", "RMSE (years)", "R²"),
  Training = c(train_mae, train_rmse, train_r2),
  Test = c(test_mae, test_rmse, test_r2)
)
# Display formatted table
kable(metrics_df, digits = 3, 
      caption = "Model Performance: Training vs Test") %>%
  kable_styling(full_width = FALSE, 
                bootstrap_options = c("striped", "hover", "condensed"))



## Data visualization


# Data frame for plot
plot_data_combined <- data.frame(
  Actual = c(y_train, y_test),
  Predicted = c(as.vector(y_train_pred), as.vector(y_test_pred))
)
# Metrics text for test set
metrics_text <- sprintf("\nR² = %.2f\nMAE = %.2f \nRMSE = %.2f", 
                        test_r2, test_mae, test_rmse)
# Scatter plot
ggplot(plot_data_combined, aes(x = Actual, y = Predicted)) +
  geom_point(alpha = 0.5, color = "steelblue", size = 2) +
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed", color = "gray30", linewidth = 0.8) +
  annotate("label",
           x = min(plot_data_combined$Actual), 
           y = max(plot_data_combined$Predicted),
           label = metrics_text, 
           hjust = 0, vjust = 1,
           size = 3.5, color = "black", 
           fill = "white", alpha = 0.8, fontface = "bold") +
  labs(title = "Elastic Net Epigenetic Clock: Predicted vs Actual Age",
       x = "Actual Age (years)", 
       y = "Predicted Age (years)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))



#The R² of 0.82 indicates that only about 82% of the variance in chronological age is explained by the model.
#The MAE of 3.50 years shows that, on average, 
#the predicted age deviates from the actual age by approximately 3 years, 
#and the RMSE of 4.46 years reflects the overall magnitude of prediction errors.


# Calculate residuals
plot_data <- data.frame(
  Actual = c(y_train, y_test),
  Predicted = c(y_train_pred, y_test_pred),
  Set = c(rep("Train", length(y_train)), rep("Test", length(y_test)))
)


plot_data$Residuals <- plot_data$Predicted - plot_data$Actual
# Create residual vs actual age plot
p1 <- ggplot(plot_data, aes(x = Actual, y = Residuals, color = Set)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  labs(title = "Residuals vs Actual Age",
       x = "Actual Age (years)", 
       y = "Residuals (years)") +
  theme_minimal() +
  theme(legend.position = "top")
# Create residual distribution plot
p2 <- ggplot(plot_data, aes(x = Residuals, fill = Set)) +
  geom_histogram(alpha = 0.6, bins = 30, position = "identity") +
  labs(title = "Residual Distribution",
       x = "Residuals (years)", 
       y = "Count") +
  theme_minimal() +
  theme(legend.position = "top")
# Display plots
grid.arrange(p1, p2, ncol = 2)




## Coefficients


# Extract coefficients
coef_mat <- coef(enet_model, s = best_lambda)

# Convert to data frame and select top 20
coef_df <- as.data.frame(as.matrix(coef_mat))
colnames(coef_df) <- "Coefficient"  
coef_df <- coef_df %>%
  rownames_to_column("Variable") %>%
  filter(Variable != "(Intercept)") %>%
  arrange(desc(abs(Coefficient))) %>%
  slice(1:20)
# Display top coefficients table
kable(coef_df, digits = 4, 
      caption = "Top 20 Features by Coefficient Magnitude") %>%
  kable_styling(full_width = FALSE, 
                bootstrap_options = c("striped", "hover", "condensed"))

# Extract coefficients
coef_mat <- coef(cv_enet, s = "lambda.min")

# Convert to data frame
coef_df <- as.data.frame(as.matrix(coef_mat)) %>%
  rownames_to_column("CpG") %>%
  rename(Weight = !!colnames(.)[2]) %>%  # Rename the second column dynamically
  filter(CpG != "(Intercept)") %>%
  arrange(desc(abs(Weight))) %>%
  slice(1:20)

# Display table 
kable(coef_df, digits = 4, caption = "Top 20 Elastic Net Model Coefficients") %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover"))





## Top variables Contributing to Age Prediction


top_cpgs <- head(coef_df, 20)  # top 20 variables

ggplot(top_cpgs, aes(x = reorder(CpG, Weight), y = Weight)) +
  geom_bar(stat = "identity", fill = "coral1") +
  coord_flip() +
  labs(title = "Top variables Contributing to Age Prediction ",
       x = "Variables",
       y = "Coefficient (Weight)") +
  theme_minimal()





g1<- ggplot(top_cpgs, aes(x = reorder(CpG, Weight), y = Weight)) +
  geom_bar(stat = "identity", fill = "coral1") +
  coord_flip() +
  labs(title = "Top variables Contributing to Age Prediction ",
       x = "Variables",
       y = "Coefficient (Weight)") +
  theme_minimal()
#save plot in png format
ggsave("variables_contribution_age.png", plot = g1, width = 12, height = 10, dpi = 150)



# Extract coefficients
coef_mat <- coef(enet_model, s = best_lambda)
# Convert to data frame and select top 20
coef_df <- as.data.frame(as.matrix(coef_mat)) %>%
  rownames_to_column("Variable") %>%
  setNames(c("Variable", "Coefficient")) %>%
  filter(Variable != "(Intercept)") %>%
  arrange(desc(abs(Coefficient))) %>%
  slice(1:20)
# Display top coefficients table
kable(coef_df, digits = 4, 
      caption = "Top 20 Features by Coefficient Magnitude") %>%
  kable_styling(full_width = FALSE, 
                bootstrap_options = c("striped", "hover", "condensed"))
#CpG-Specific Analysis
# Select only CpGs (exclude gender variables)
top_cpgs <- coef_df %>%
  filter(grepl("^cg", Variable)) %>%
  arrange(desc(abs(Coefficient)))
# Classify as hypomethylated or hypermethylated
top_cpgs <- top_cpgs %>%
  mutate(Methylation_Change = ifelse(Coefficient > 0, "Hypermethylated", "Hypomethylated"))
# Display number of CpGs in top features
cat("Number of CpGs in top 20 features:", nrow(top_cpgs), "\n")
#Top CpGs Visualization (Colored by Methylation)
# Plot with methylation direction
p_cpgs_colored <- ggplot(top_cpgs, aes(x = reorder(Variable, Coefficient), 
                                       y = Coefficient, 
                                       fill = Methylation_Change)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Hypermethylated" = "firebrick2", 
                               "Hypomethylated" = "steelblue")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = "Top CpGs for Age Prediction",
    x = "CpG Sites",
    y = "Elastic Net Coefficient (Weight)",
    fill = "Methylation Change"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 9)
  )

print(p_cpgs_colored)
#Top CpGs Visualization (Single Color)
# Plot with single color
p_cpgs_simple <- ggplot(top_cpgs, aes(x = reorder(Variable, Coefficient), 
                                      y = Coefficient)) +
  geom_bar(stat = "identity", fill = "coral1") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = "Top CpGs Contributing to Age Prediction",
    x = "CpG Sites",
    y = "Elastic Net Coefficient (Weight)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 9)
  )

print(p_cpgs_simple)

#Save Plots
# Create metrics text for combined plot
metrics_text <- sprintf("Test Metrics:\nR² = %.3f\nMAE = %.2f years\nRMSE = %.2f years", 
                        test_r2, test_mae, test_rmse)
# Plot 1: Predicted vs Actual (all points combined)
p1 <- ggplot(plot_data_combined, aes(x = Actual, y = Predicted)) +
  geom_point(alpha = 0.5, color = "steelblue", size = 2) +
  geom_smooth(method = "lm", color = "magenta", se = FALSE) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray30") +
  annotate("label",
           x = min(plot_data_combined$Actual), 
           y = max(plot_data_combined$Predicted),
           label = metrics_text, 
           hjust = 0, vjust = 1,
           size = 3.5, color = "black", 
           fill = "white", alpha = 0.8, fontface = "bold") +
  labs(
    title = "Chronological vs Biological Age",
    x = "Chronological Age (years)",    
    y = "Biological Age (years)"        
  ) +
  theme_minimal(base_size = 13) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
# Plot 2: Predicted vs Actual (Training vs Test)
p2 <- ggplot(plot_data, aes(x = Actual, y = Predicted, color = Set)) +
  geom_point(alpha = 0.6, size = 2.5) +
  geom_abline(slope = 1, intercept = 0, 
              linetype = "dashed", color = "gray30", linewidth = 0.8) +
  annotate("label",
           x = min(plot_data$Actual), 
           y = max(plot_data$Predicted),
           label = metrics_text, 
           hjust = 0, vjust = 1,
           size = 3.5, color = "black", 
           fill = "white", alpha = 0.8, fontface = "bold") +
  labs(title = "Elastic Net Epigenetic Clock: Predicted vs Actual Age",
       x = "Actual Age (years)", 
       y = "Predicted Age (years)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "top") +
  scale_color_manual(values = c("Training" = "#3498db", "Test" = "#e74c3c"))
# Plot 3: Top CpGs (colored by methylation)
p3 <- ggplot(top_cpgs, aes(x = reorder(Variable, Coefficient), 
                           y = Coefficient, 
                           fill = Methylation_Change)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_manual(values = c("Hypermethylated" = "firebrick2", 
                               "Hypomethylated" = "steelblue")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = "Top CpGs for Age Prediction",
    x = "CpG Sites",
    y = "Elastic Net Coefficient (Weight)",
    fill = "Methylation Change"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 9)
  )
# Plot 4: Top CpGs (simple version)
p4 <- ggplot(top_cpgs, aes(x = reorder(Variable, Coefficient), 
                           y = Coefficient)) +
  geom_bar(stat = "identity", fill = "coral1") +
  coord_flip() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = "Top CpGs Contributing to Age Prediction",
    x = "CpG Sites",
    y = "Elastic Net Coefficient (Weight)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 9)
  )
# Save all plots
ggsave("chronological_vs_biological_age.png", plot = p1, 
       width = 10, height = 8, dpi = 300)

ggsave("predicted_vs_actual_train_test.png", plot = p2, 
       width = 10, height = 8, dpi = 300)

ggsave("top_cpgs_methylation_colored.png", plot = p3, 
       width = 12, height = 10, dpi = 300)

ggsave("top_cpgs_simple.png", plot = p4, 
       width = 12, height = 10, dpi = 300)

all_coef_full <- as.data.frame(as.matrix(coef(cv_enet, s = "lambda.min"))) %>%
  rownames_to_column("Variable") %>%
  rename(Coefficient = 2) %>%   #
  filter(Variable != "(Intercept)")   #remove intercept

# Filter gender coefficient
gender_coef <- all_coef_full %>%
  filter(grepl("gender", Variable, ignore.case = TRUE))


# Table with coefficients (gender)
cat("\nGender Covariate Coefficients (All):\n")
kable(gender_coef, digits = 4, caption = "Gender Covariate Coefficients") %>%
  kable_styling(full_width = FALSE, position = "center")




## Select top 20cpgs 

# Extract coefficients from cv.glmnet at lambda.min
coef_mat <- coef(cv_enet, s = "lambda.min")


top20_cpgs <- as.data.frame(as.matrix(coef_mat)) %>%
  rownames_to_column("CpG") %>%
  rename(Weight = !!colnames(.)[2]) %>%   
  filter(CpG != "(Intercept)" & !grepl("^gender", CpG)) %>%  # exclude gender
  arrange(desc(abs(Weight))) %>%
  slice(1:20) %>%
  mutate(Methylation_Change = ifelse(Weight > 0,
                                     "Hypermethylated",
                                     "Hypomethylated"))

# Save to CSV
write.csv(dplyr::select(top20_cpgs, CpG, Weight, Methylation_Change),
          "Top20_CpGs.csv",
          row.names = FALSE)

# Extract coefficients from the trained model (not cv)
coef_mat <- coef(enet_model, s = best_lambda)

top20_cpgs <- as.data.frame(as.matrix(coef_mat)) %>%
  rownames_to_column("CpG") %>%
  setNames(c("CpG", "Weight")) %>%
  filter(CpG != "(Intercept)" & !grepl("^gender", CpG)) %>%  # Exclude intercept and gender
  arrange(desc(abs(Weight))) %>%
  slice(1:20) %>%
  mutate(Methylation_Change = ifelse(Weight > 0,
                                     "Hypermethylated",
                                     "Hypomethylated"))
# Display the top 20 CpGs
kable(top20_cpgs, digits = 4, 
      caption = "Top 20 CpGs Selected by Elastic Net") %>%
  kable_styling(full_width = FALSE, 
                bootstrap_options = c("striped", "hover"))
# Save to CSV
write.csv(dplyr::select(top20_cpgs, CpG, Weight, Methylation_Change),
          "Top20_CpGs_ElasticNet.csv",
          row.names = FALSE)







The predicted ages from the Elastic Net model (Pred_ENet) were saved for each sample. These predictions are used for statistical comparisons between models, including paired tests, correlation analysis, and Bland–Altman plots to assess consistency and potential bias.

# Predictions from the Elastic Net model
pred_enet <- predict(cv_enet, newx = X_pred, s = "lambda.min")

# Create a data.frame that
results_enet <- data.frame(
  SampleID = subset_500$X,
  Observed_Age = y_age,
  Pred_ENet = as.numeric(pred_enet)
)

# Save to CSV
write.csv(results_enet, "ENet_Predicted_Ages.csv", row.names = FALSE)




