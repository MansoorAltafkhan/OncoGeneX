import pandas as pd
import os

metadata_path = "data/raw/GSE40791_final_metadata.xlsx"
expression_path = "data/raw/GSE40791_gene_expression_matrix_GSM.csv"

print("Checking files...\n")
print("Metadata file exists:", os.path.exists(metadata_path))
print("Expression matrix exists:", os.path.exists(expression_path))

metadata = pd.read_excel(metadata_path)

expression = pd.read_csv(expression_path)

print("\n--- METADATA ---")
print("Shape:", metadata.shape)
print("\nColumns:")
print(metadata.columns.tolist())
print("\nFirst 5 rows:")
print(metadata.head())

print("\n--- GENE EXPRESSION MATRIX ---")
print("Shape:", expression.shape)
print("\nFirst 5 rows:")
print(expression.head())

expression_sample_ids = set(expression.columns[1:])

print("\n--- SAMPLE ID CHECK ---")
print("Number of expression samples:", len(expression_sample_ids))

print("\nMetadata columns:")
print(metadata.columns.tolist())

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
