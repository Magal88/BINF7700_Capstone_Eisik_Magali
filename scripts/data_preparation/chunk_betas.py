import pandas as pd

# Path to your file
file_path = "GSE55763_normalized_betas.txt.gz"

# Specify chunk size (number of rows per chunk)
chunksize = 100  

# Create an iterator
chunks = pd.read_csv(file_path, sep="\t", chunksize=chunksize)

for i, chunk in enumerate(chunks):
    # For example, save each chunk to a separate CSV
    chunk.to_csv(f"beta_chunk_{i}.csv.gz", index=False, compression='gzip')
