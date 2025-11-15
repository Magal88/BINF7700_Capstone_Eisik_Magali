 "Bland-Altman Plot"


#Load required packages
library(knitr)
library(kableExtra)
library(ggplot2)
library(blandr)



#Save metrics of each models in a data frame
metrics_models <- data.frame(
  Model = c("Elastic Net", "Random Forest"),
  R2 = c(0.82, 0.60),
  MAE = c(3.50, 5.03),
  RMSE = c(4.46, 6.19)
)
# Display table
metrics_models %>%
  kbl(
    caption = "<b style='text-align:center;'>Performance Metrics for Machine Learning Models</b>",
    align = "c",
    escape = FALSE  
  ) %>%
  kable_styling(
    full_width = FALSE,
    position = "center",
    bootstrap_options = c("striped", "hover")
  )

#Load datasets with predicted age
rf_pred <- read.csv("RF_Predicted_Ages.csv", header = TRUE)
enet_pred <- read.csv("ENet_Predicted_Ages.csv", header = TRUE)

#Verify dimensions and first 3 rows
dim(rf_pred)
dim(enet_pred)
head(rf_pred,3)
head(enet_pred,3)

#Save in data frame
pred_comparison <- data.frame(
  Pred_RF = rf_pred$Pred_RF,
  Pred_ENet = enet_pred$Pred_ENet
)





p <- blandr.draw(
  pred_comparison$Pred_RF,
  pred_comparison$Pred_ENet,
  method1name = "Random Forest",
  method2name = "Elastic Net",
  plotTitle = "Bland-Altman Plot: RF vs ENet Age Predictions"
)

p + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size=14, face="bold")) +
  geom_point(color="steelblue", size=2) +
  geom_hline(yintercept=0, linetype="dashed", color="coral1") 
```

#Statistics
ba_stats <- blandr.statistics(
  pred_comparison$Pred_RF,
  pred_comparison$Pred_ENet
)
print(ba_stats)


