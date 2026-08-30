import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# ============================================================
# 1. CREATE OUTPUT DIRECTORY
# ============================================================

os.makedirs("results/visualization", exist_ok=True)


# ============================================================
# 2. LOAD PROCESSED EXPRESSION DATA
# ============================================================

expression = pd.read_csv(
    "data/processed/processed_expression.csv",
    index_col=0
)

# Convert column names to strings
expression.columns = expression.columns.astype(str).str.strip()


# ============================================================
# 3. LOAD UPREGULATED AND DOWNREGULATED GENES
# ============================================================

up_df = pd.read_csv(
    "results/upregulated_genes.csv"
)

down_df = pd.read_csv(
    "results/downregulated_genes.csv"
)


print("Expression data shape:", expression.shape)

print("\nUpregulated columns:")
print(up_df.columns.tolist())

print("\nDownregulated columns:")
print(down_df.columns.tolist())


# ============================================================
# 4. EXTRACT GENE IDs
# ============================================================

def extract_gene_data(df):

    # The 'Gene' column contains IDs/probe IDs
    if "Gene" not in df.columns:

        raise ValueError(
            "'Gene' column not found. "
            f"Available columns: {df.columns.tolist()}"
        )

    # Get gene/probe IDs
    gene_ids = (
        df["Gene"]
        .dropna()
        .astype(str)
        .str.strip()
        .drop_duplicates()
        .tolist()
    )

    return gene_ids


# ============================================================
# 5. SELECT TOP DIFFERENTIALLY EXPRESSED GENES
# ============================================================

# Top 15 upregulated genes
top_up_ids = extract_gene_data(
    up_df.head(15)
)

# Top 15 downregulated genes
top_down_ids = extract_gene_data(
    down_df.head(15)
)

# Combine
top_gene_ids = (
    top_up_ids +
    top_down_ids
)


print("\nTop genes selected:", len(top_gene_ids))


# ============================================================
# 6. FIND GENE IDs IN EXPRESSION MATRIX
# ============================================================

available_genes = [

    gene

    for gene in top_gene_ids

    if gene in expression.columns

]


print(
    "Genes found in expression matrix:",
    len(available_genes)
)

print(
    "Genes not found:",
    len(top_gene_ids) - len(available_genes)
)


if len(available_genes) == 0:

    print("\nFirst 10 expression columns:")
    print(expression.columns[:10].tolist())

    print("\nFirst 10 DEG Gene IDs:")
    print(top_gene_ids[:10])

    raise ValueError(
        "None of the selected Gene IDs were found "
        "in the expression matrix."
    )


# ============================================================
# 7. CREATE GENE LABELS
# ============================================================

# Create mapping from Gene ID to Gene Symbol
gene_labels = {}

combined_df = pd.concat(
    [up_df, down_df]
)

for _, row in combined_df.iterrows():

    gene_id = str(row["Gene"]).strip()

    if "Gene Symbol" in combined_df.columns:

        gene_symbol = str(
            row["Gene Symbol"]
        ).strip()

        if pd.notna(row["Gene Symbol"]) and gene_symbol != "nan":

            gene_labels[gene_id] = (
                gene_symbol +
                f" ({gene_id})"
            )

        else:

            gene_labels[gene_id] = gene_id

    else:

        gene_labels[gene_id] = gene_id


# ============================================================
# 8. CREATE HEATMAP DATA
# ============================================================

# Select available genes
heatmap_data = expression[
    available_genes
].T


# Replace probe IDs with Gene Symbols for display
heatmap_data.index = [

    gene_labels.get(
        gene,
        gene
    )

    for gene in available_genes

]


# ============================================================
# 9. Z-SCORE NORMALIZATION
# ============================================================

heatmap_zscore = heatmap_data.apply(

    lambda x: (
        x - x.mean()
    ) / x.std(),

    axis=1

)


# Replace missing/infinite values
heatmap_zscore = heatmap_zscore.replace(

    [np.inf, -np.inf],

    np.nan

)

heatmap_zscore = heatmap_zscore.fillna(0)


# ============================================================
# 10. CREATE HEATMAP
# ============================================================

plt.figure(
    figsize=(14, 10)
)


plt.imshow(

    heatmap_zscore,

    aspect="auto",

    interpolation="nearest"

)


plt.colorbar(
    label="Z-score"
)


plt.yticks(

    range(
        len(heatmap_zscore.index)
    ),

    heatmap_zscore.index,

    fontsize=8

)


# Hide sample labels to avoid overcrowding
plt.xticks([])


plt.xlabel(
    "Samples"
)


plt.ylabel(
    "Genes"
)


plt.title(
    "Heatmap of Top Upregulated and Downregulated Genes"
)


plt.tight_layout()


# ============================================================
# 11. SAVE HEATMAP
# ============================================================

output_file = (
    "results/visualization/"
    "top_DEG_heatmap.png"
)


plt.savefig(

    output_file,

    dpi=300,

    bbox_inches="tight"

)


plt.close()


# ============================================================
# 12. COMPLETE
# ============================================================

print("\n==========================================")

print(
    "DEG HEATMAP ANALYSIS COMPLETE"
)

print("==========================================")

print(
    f"\nHeatmap saved to: "
    f"{output_file}"
)