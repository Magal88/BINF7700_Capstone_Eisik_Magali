#!/usr/bin/env python3
# extract_metadata_age.py
# Extract metadata (age, sex, tissue, dataset) from GSE55763 and save as CSV

import GEOparse
import pandas as pd

# Download and parse GSE55763
gse = GEOparse.get_GEO(geo="GSE55763", destdir="./data", silent=True, how="full")

# Parse all characteristics into a dictionary
metadata = {}
for gsm_name, gsm in gse.gsms.items():
    metadata[gsm_name] = {}
    for char in gsm.metadata.get('characteristics_ch1', []):
        parts = char.split(':', 1)  #
        if len(parts) == 2:
            key, value = parts
            key = key.strip().lower()
            value = value.strip()
            metadata[gsm_name][key] = value

# Convert the dictionary to a DataFrame
df_metadata = pd.DataFrame.from_dict(metadata, orient='index')

# Convert age to numeric
if 'age' in df_metadata.columns:
    df_metadata['age'] = pd.to_numeric(df_metadata['age'], errors='coerce')

# Keep only samples with age and gender (if columns exist)
if 'age' in df_metadata.columns and 'gender' in df_metadata.columns:
    df_filtered = df_metadata.dropna(subset=['age', 'gender'])
else:
    df_filtered = df_metadata

# Save the cleaned metadata to CSV
df_filtered.to_csv('./data/GSE55763_metadata_clean.csv', index=True)

# Preview first rows
print(df_filtered.head())
      
