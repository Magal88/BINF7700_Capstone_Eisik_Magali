import pandas as pd
import glob

# ---------------------------
# Configuration
# ---------------------------
chunks_path = "beta_chunk_*.csv.gz"       # All chunk files
output_file = "GSE55763_top1000_all_samples.csv"  # Final CSV output

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
# Select top 1000 most variable CpGs
# ---------------------------
top1000 = all_data.loc[variances.nlargest(1000).index, :]

# ---------------------------
# Save final CSV
# ---------------------------
top1000.to_csv(output_file, index=False)
print(f"Done! Top 1000 CpGs saved to {output_file}")
