import GEOparse
import pandas as pd
import os

# Path to the .soft.gz file located inside the data folder
soft_path = "data/GSE55763_family.soft.gz"

# Load the GEO file locally
gse = GEOparse.get_GEO(filepath=soft_path)

# Extract metadata from each GSM sample
metadata = []
for gsm_name, gsm in gse.gsms.items():
    sample_info = gsm.metadata
    metadata.append({
        "GSM_ID": gsm_name,
        "title": sample_info.get("title", [""])[0],
        "source_name": sample_info.get("source_name_ch1", [""])[0],
        "age": sample_info.get("characteristics_ch1", [""])[0] if len(sample_info.get("characteristics_ch1", [])) > 0 else "",
        "sex": sample_info.get("characteristics_ch1", [""])[1] if len(sample_info.get("characteristics_ch1", [])) > 1 else "",
    })

# Convert list of dictionaries into a DataFrame
meta_df = pd.DataFrame(metadata)

# Clean text in metadata fields
meta_df["age"] = meta_df["age"].str.replace("age: ", "", regex=False)
meta_df["sex"] = meta_df["sex"].str.replace("sex: ", "", regex=False)

# Save cleaned metadata as CSV inside the data folder
output_path = os.path.join("data", "GSE55763_metadata_clean.csv")
meta_df.to_csv(output_path, index=False)
(base) [eisik.m@explorer-02 capstone]$ cat soft_to_csv.py
import GEOparse
import pandas as pd

# Load the soft file (already downloaded/unzipped)
gse = GEOparse.get_GEO(filepath="GSE55763_family.soft", silent=True)

# Extract metadata
metadata = {}
for gsm_name, gsm in gse.gsms.items():
    for key, value in gsm.metadata.items():
        if key == 'characteristics_ch1':
            for tmp in value:
                splitUp = [i.strip() for i in tmp.split(':')]
                if len(splitUp) == 2:
                    if splitUp[0] not in metadata:
                        metadata[splitUp[0]] = {}
                    metadata[splitUp[0]][gsm_name] = splitUp[1]

df_metadata = pd.DataFrame(metadata).transpose()  # probes -> rows
df_metadata.index.name = 'GSM_ID'

# Save to CSV
df_metadata.to_csv("GSE55763_family.csv")