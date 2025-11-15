# "Cohort summary GSE55763"

#Load required packages

library(dplyr)
library(knitr)
library(kableExtra)



#Load data set
beta_data <- read.csv("GSE55763_beta_500_imputed.csv", header = TRUE, stringsAsFactors = FALSE)



age_summary <- beta_data %>%
  summarise(
    n = n(),
    Mean = round(mean(age, na.rm = TRUE), 1),
    SD = round(sd(age, na.rm = TRUE), 1),
    Median = median(age, na.rm = TRUE),
    Min = min(age, na.rm = TRUE),
    Q1 = quantile(age, 0.25, na.rm = TRUE),
    Q3 = quantile(age, 0.75, na.rm = TRUE),
    Max = max(age, na.rm = TRUE)
  )

kable(age_summary, caption = "<b>Age summary</b>", format = "html") %>%
  kable_styling(
    full_width = FALSE, 
    position = "center", 
    bootstrap_options = c("striped", "hover", "condensed")
  )

gender_summary <- beta_data %>%
  group_by(gender) %>%
  summarise(count = n()) %>%
  mutate(percent = round(count / sum(count) * 100, 1))

kable(gender_summary, caption = "<b>Gender summary</b>", format = "html") %>%
  kable_styling(
    full_width = FALSE,
    position = "center",
    bootstrap_options = c("striped", "hover", "condensed")
  )

age_by_gender <- beta_data %>%
  group_by(gender) %>%
  summarise(
    n = n(),
    Mean = round(mean(age, na.rm = TRUE), 1),
    SD = round(sd(age, na.rm = TRUE), 1),
    Median = round(median(age, na.rm = TRUE), 1),
    Min = round(min(age, na.rm = TRUE), 1),
    Q1 = round(quantile(age, 0.25, na.rm = TRUE), 1),
    Q3 = round(quantile(age, 0.75, na.rm = TRUE), 1),
    Max = round(max(age, na.rm = TRUE), 1)
  )
kable(age_by_gender, caption = "<b>Age Summary by Gender</b>", format = "html") %>%
  kable_styling(
    full_width = FALSE,
    position = "center",
    bootstrap_options = c("striped", "hover", "condensed")
  )


