import os
import pandas as pd

metadata_path = "data/raw/GSE40791_final_metadata.xlsx"
expression_path = "data/raw/GSE40791_gene_expression_matrix_GSM.csv"

print("Checking files...\n")

if not os.path.exists(metadata_path):
    raise FileNotFoundError(f"Metadata file not found: {metadata_path}")

if not os.path.exists(expression_path):
    raise FileNotFoundError(f"Expression file not found: {expression_path}")

metadata = pd.read_excel(metadata_path)
expression = pd.read_csv(expression_path)

print("--- METADATA ---")
print("Shape:", metadata.shape)
print("Columns:", metadata.columns.tolist())
print(metadata.head())

print("\n--- GENE EXPRESSION MATRIX ---")
print("Shape:", expression.shape)
print(expression.head())

required_metadata_columns = ["gsm_id", "group"]

missing_columns = [
    column for column in required_metadata_columns
    if column not in metadata.columns
]

if missing_columns:
    raise ValueError(
        f"Metadata is missing required columns: {missing_columns}"
    )

metadata["gsm_id"] = metadata["gsm_id"].astype(str).str.strip()
expression.columns = expression.columns.astype(str).str.strip()

if metadata["gsm_id"].isin(["", "nan", "None"]).any():
    raise ValueError("Validation failed: metadata contains empty sample IDs.")

duplicate_metadata_ids = metadata[
    metadata["gsm_id"].duplicated(keep=False)
]["gsm_id"].unique()

if len(duplicate_metadata_ids) > 0:
    raise ValueError(
        f"Validation failed: duplicate metadata sample IDs found: "
        f"{list(duplicate_metadata_ids)}"
    )

if expression.shape[1] < 2:
    raise ValueError(
        "Validation failed: expression matrix must contain "
        "a gene/probe column and sample columns."
    )

expression_sample_ids = set(expression.columns[1:])
metadata_sample_ids = set(metadata["gsm_id"])

if "" in expression_sample_ids:
    raise ValueError(
        "Validation failed: expression matrix contains an empty sample ID."
    )

matching_samples = expression_sample_ids.intersection(metadata_sample_ids)
missing_in_metadata = expression_sample_ids - metadata_sample_ids
missing_in_expression = metadata_sample_ids - expression_sample_ids

print("\n--- SAMPLE ID MATCHING RESULTS ---")
print("Expression samples:", len(expression_sample_ids))
print("Metadata samples:", len(metadata_sample_ids))
print("Matching samples:", len(matching_samples))

if missing_in_metadata:
    print("\nSamples in expression but NOT in metadata:")
    print(list(missing_in_metadata)[:10])

if missing_in_expression:
    print("\nSamples in metadata but NOT in expression:")
    print(list(missing_in_expression)[:10])

if missing_in_metadata or missing_in_expression:
    raise ValueError(
        "Validation failed: expression and metadata sample IDs do not match."
    )

expression_values = expression.iloc[:, 1:].apply(
    pd.to_numeric,
    errors="coerce"
)

if expression_values.isna().any().any():
    missing_values = int(expression_values.isna().sum().sum())
    raise ValueError(
        f"Validation failed: {missing_values} non-numeric or missing "
        "expression values detected."
    )

gene_ids = expression.iloc[:, 0].astype(str).str.strip()

if gene_ids.isin(["", "nan", "None"]).any():
    raise ValueError(
        "Validation failed: empty gene/probe IDs detected."
    )

duplicate_gene_ids = gene_ids[gene_ids.duplicated()].unique()

if len(duplicate_gene_ids) > 0:
    print(
        f"\nWarning: {len(duplicate_gene_ids)} duplicate gene/probe IDs found."
    )

groups = metadata["group"].astype(str).str.strip().str.lower()

if groups.isin(["", "nan", "none"]).any():
    raise ValueError(
        "Validation failed: missing group labels detected."
    )

print("\n--- GROUP INFORMATION ---")
print(metadata["group"].value_counts())

print("\n--- DATA VALIDATION COMPLETE ---")
print("✓ Files exist")
print("✓ Required metadata columns present")
print("✓ Sample IDs are valid")
print("✓ No duplicate metadata sample IDs")
print("✓ Expression and metadata samples match")
print("✓ Expression values are numeric and complete")
print("✓ Gene/probe IDs are present")
print("✓ Group labels are present")
print("\nAll critical validation checks passed.")