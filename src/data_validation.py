import pandas as pd
import os

# File paths
metadata_path = "data/raw/GSE40791_final_metadata.xlsx"
expression_path = "data/raw/GSE40791_gene_expression_matrix_GSM.csv"

# Check if files exist
print("Checking files...\n")
print("Metadata file exists:", os.path.exists(metadata_path))
print("Expression matrix exists:", os.path.exists(expression_path))

# Load metadata
metadata = pd.read_excel(metadata_path)

# Load gene expression matrix
expression = pd.read_csv(expression_path)

# Display metadata information
print("\n--- METADATA ---")
print("Shape:", metadata.shape)
print("\nColumns:")
print(metadata.columns.tolist())
print("\nFirst 5 rows:")
print(metadata.head())

# Display expression matrix information
print("\n--- GENE EXPRESSION MATRIX ---")
print("Shape:", expression.shape)
print("\nFirst 5 rows:")
print(expression.head())
# Extract sample IDs from expression matrix
expression_sample_ids = set(expression.columns[1:])

# Display metadata sample IDs
print("\n--- SAMPLE ID CHECK ---")
print("Number of expression samples:", len(expression_sample_ids))

# Show metadata columns to identify the sample ID column
print("\nMetadata columns:")
print(metadata.columns.tolist())
# Compare sample IDs between metadata and expression matrix

metadata_sample_ids = set(metadata["gsm_id"].astype(str))

matching_samples = expression_sample_ids.intersection(metadata_sample_ids)
missing_in_metadata = expression_sample_ids - metadata_sample_ids
missing_in_expression = metadata_sample_ids - expression_sample_ids

print("\n--- SAMPLE ID MATCHING RESULTS ---")

print("Expression samples:", len(expression_sample_ids))
print("Metadata samples:", len(metadata_sample_ids))
print("Matching samples:", len(matching_samples))

print("\nSamples in expression but NOT in metadata:")
print(list(missing_in_metadata)[:10])

print("\nSamples in metadata but NOT in expression:")
print(list(missing_in_expression)[:10])