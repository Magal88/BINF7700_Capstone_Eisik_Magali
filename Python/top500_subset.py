import pandas as pd
import glob

# ---------------------------
# Configuration
# ---------------------------
chunks_path = "beta_chunk_*.csv.gz"       # All chunks
output_file = "GSE55763_top500_all_samples.csv"  # Name of the final CSV

# ---------------------------
# Read and concatenate all chunks
# ---------------------------
chunk_files = sorted(glob.glob(chunks_path))
print(f"Found {len(chunk_files)} chunks.")

all_chunks = []
for f in chunk_files:
    print(f"Reading {f}")
    chunk = pd.read_csv(f)
    all_chunks.append(chunk)

all_data = pd.concat(all_chunks, ignore_index=True)
print("All chunks concatenated.")

# ---------------------------
# Calculate variance per row (CpG)
# ---------------------------
# Assumes the first column is CpG_ID
variances = all_data.iloc[:, 1:].var(axis=1)

# ---------------------------
# Select the top 500 most variable CpGs
# ---------------------------
top500 = all_data.loc[variances.nlargest(500).index, :]

# ---------------------------
# Save the final CSV
# ---------------------------
top500.to_csv(output_file, index=False)

