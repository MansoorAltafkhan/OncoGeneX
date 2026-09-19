import pandas as pd
import numpy as np

expression = pd.read_csv(
    "data/raw/GSE40791_gene_expression_matrix_GSM.csv",
    index_col=0
)

metadata = pd.read_excel(
    "data/raw/GSE40791_final_metadata.xlsx"
)

print("--- DATA PREPROCESSING ---")

expression = expression.T

print("Expression shape after transpose:", expression.shape)

expression = expression.loc[metadata["gsm_id"]]

print("Aligned expression shape:", expression.shape)

print("\nMissing values:")
print(expression.isnull().sum().sum())

import os
os.makedirs("data/processed", exist_ok=True)

expression.to_csv("data/processed/processed_expression.csv")

print("\nPreprocessing completed successfully!")
