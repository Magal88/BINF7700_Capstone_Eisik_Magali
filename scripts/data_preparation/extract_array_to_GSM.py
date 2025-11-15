# extract_array_to_gsm.py
import pandas as pd

soft_file = "GSE55763_family.soft"
output_csv = "GSE55763_family.csv"

array_ids = []
gsm_ids = []

with open(soft_file, "r") as f:
    in_table = False
    for line in f:
        line = line.strip()
        # Start of the sample table
        if line.startswith("!Sample_table_begin"):
            in_table = True
            continue
        # End of the sample table
        if line.startswith("!Sample_table_end"):
            in_table = False
            continue
        if in_table:
            # Skip comment lines
            if line.startswith("!"):
                continue
            # Split the line by tab or whitespace
            parts = line.split("\t")
            if len(parts) >= 2:
                array_ids.append(parts[0])
                gsm_ids.append(parts[1])

# Create DataFrame
df = pd.DataFrame({
    "Array_ID": array_ids,
    "GSM_ID": gsm_ids
})

# Save to CSV
df.to_csv(output_csv, index=False)
print(f"Mapping saved to {output_csv}, {len(df)} entries.")
