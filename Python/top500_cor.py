import pandas as pd
import numpy as np
import glob
from scipy.stats import pearsonr

# --------------------------
# Load metadata
# --------------------------

metadata = pd.read_csv("GSE55763_metadata_clean.csv")
metadata_dict = dict(zip(metadata['GSM'], metadata['Age']))

# --------------------------
# Correlation
# --------------------------
def compute_correlations(chunk_file, metadata_dict):
    # Read chunk (CpGs x samples)
    chunk = pd.read_csv(chunk_file, index_col=0)


    common_samples = [s for s in chunk.columns if s in metadata_dict]
    chunk = chunk[common_samples]


    age_vector = np.array([metadata_dict[s] for s in common_samples])

    # Compute Pearson correlation for each CpG
    cor_list = []
    for cpg, values in chunk.iterrows():
        r, _ = pearsonr(values.values, age_vector)
        cor_list.append((cpg, r))

    return pd.DataFrame(cor_list, columns=['CpG', 'Correlation'])

# --------------------------
#  Read all chunks
# --------------------------
all_cor = []

chunk_files = sorted(glob.glob("beta_chunk_*.csv.gz"))  # gzip CSVs
for chunk_file in chunk_files:
    print("Processing:", chunk_file)
    cor_df = compute_correlations(chunk_file, metadata_dict)
    all_cor.append(cor_df)

# Concatenate all results
all_cor_df = pd.concat(all_cor, ignore_index=True)

# --------------------------
#  Select top 500 CpGs
# --------------------------
all_cor_df['absCorrelation'] = all_cor_df['Correlation'].abs()
top500 = all_cor_df.sort_values('absCorrelation', ascending=False).head(500)

# Save results
top500.to_csv("top500_cpgs.csv", index=False)