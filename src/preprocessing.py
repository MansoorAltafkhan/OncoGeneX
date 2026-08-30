import pandas as pd
import numpy as np

# Load the gene expression data
expression = pd.read_csv(
    "data/raw/GSE40791_gene_expression_matrix_GSM.csv",
    index_col=0
)

# Load metadata
metadata = pd.read_excel(
    "data/raw/GSE40791_final_metadata.xlsx"
)

print("--- DATA PREPROCESSING ---")

# Transpose expression data
# Rows = samples, Columns = genes
expression = expression.T

print("Expression shape after transpose:", expression.shape)

# Make sure expression samples follow metadata order
expression = expression.loc[metadata["gsm_id"]]

print("Aligned expression shape:", expression.shape)

# Check for missing values
print("\nMissing values:")
print(expression.isnull().sum().sum())

# Create processed folder if it doesn't exist
import os
os.makedirs("data/processed", exist_ok=True)

# Save processed data
expression.to_csv("data/processed/processed_expression.csv")

print("\nPreprocessing completed successfully!")