import pandas as pd

# Load mapping files
array_to_gsm = pd.read_csv("GSE55763_family.csv")      # columns: Array_ID, GSM_ID
metadata = pd.read_csv("GSE55763_metadata_clean.csv")  # columns: GSM_ID, age

# Load a beta-value chunk
beta_chunk = pd.read_csv("beta_chunk_0.csv")  # columns: Array_IDs

# 1️⃣ Merge Array_ID → GSM_ID → age
mapping = pd.merge(array_to_gsm, metadata[['GSM_ID','age']], on='GSM_ID', how='left')

# 2️⃣ Reorder age to match beta columns
age_ordered = mapping.set_index('Array_ID').loc[beta_chunk.columns].age

# 3️⃣ Now you have a vector of ages corresponding to the beta columns
print(age_ordered.head())

# 4️⃣ Optional: rename beta columns to GSM_ID or leave as Array_ID
beta_chunk.columns = beta_chunk.columns  # keep Array_IDs or beta_chunk.columns = mapping['GSM_ID']

