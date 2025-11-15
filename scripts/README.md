The scripts/ directory is organized into four subfolders, each representing a specific stage of the data analysis pipeline.

1) **data_processing:**
Cleaning and restructuring of raw GEO and beta-value data, including extraction, formatting, and reading in chunks for initial quality processing.
Packages: GEOparse, glob, numpy, pandas (Python)
2) **data_preparation:**
Preparation of cleaned data for modeling.
Includes feature selection, metadata mapping, M-value conversion, and dataset assembly.
Packages: dplyr, tidyverse, tibble,tidyr,readr, kableExtra, knitr, ggplot2 (R)
3) **modeling**
Training and evaluation of machine-learning models.
Includes Elastic Net, Random Forest, SHAP interpretation, Bland–Altman plots for model comparison, tuning, and performance assessment.
Packages: glmnet, caret, ggplot2, gridExtra, dplyr, tidyverse, tibble,tidyr,readr, kableExtra, knitr, blandr (R)
Packages: numpy, pandas, matplotlib, shap, sklearn (Python)
4) **downstream_analysis:**
Annotated the top 20 CpGs to genes, separated hyper- and hypomethylated sites, and performed KEGG pathway enrichment with visualization via dotplots, cnetplots, and barplots.
Packages: dplyr, tidyr, ggplot2, kableExtra, knitr, IlluminaHumanMethylation450kanno.ilmn12.hg19, clusterProfiler, org.Hs.eg.db, enrichplot, readr, patchwork.
