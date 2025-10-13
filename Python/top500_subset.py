import pandas as pd

df = pd.read_csv("GSE55763_normalized_betas.txt.gz", sep="\t", nrows=10000)  # first 10k rows
# compute variance across samples (skip ID column)
variances = df.iloc[:, 1:].var(axis=1)
top500 = df.iloc[variances.nlargest(500).index, :101]
top500.to_csv("GSE55763_top500_var.csv", index=False)
