from pathlib import Path
import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

EXPRESSION_FILE = Path("data/processed/processed_expression.csv")
METADATA_FILE = Path("data/raw/GSE40791_final_metadata.xlsx")
ANNOTATION_FILE = Path("data/raw/GPL570-55999.txt")

RESULTS_DIR = Path("results")
PLOTS_DIR = RESULTS_DIR / "plots"

FDR_THRESHOLD = 0.05
LOG2FC_THRESHOLD = 1.0

HEATMAP_TOP_N = 20
PCA_TOP_VARIABLE_GENES = 2000

def ensure_required_columns(df, columns, name):
    """Raise a clear error if required columns are missing."""
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise ValueError(
            f"{name} is missing required columns: {missing}\n"
            f"Available columns: {df.columns.tolist()}"
        )

def clean_gene_symbol(value):
    """
    Clean GPL570 gene symbols.
    GEO annotations can contain multiple symbols separated by ' /// '.
    """
    if pd.isna(value):
        return ""

    value = str(value).strip()

    if value in {"", "---", "nan", "NA"}:
        return ""

    return value.split(" /// ")[0].strip()

def zscore_rows(df):
    """
    Row-wise Z-score for a gene x sample heatmap.
    Uses ddof=0 and safely handles constant rows.
    """
    mean = df.mean(axis=1)
    std = df.std(axis=1, ddof=0).replace(0, np.nan)

    scaled = df.sub(mean, axis=0).div(std, axis=0)
    return scaled.fillna(0)

def make_unique_labels(labels):
    """
    Make duplicated labels unique while retaining the readable gene symbol.
    Example: CA4, CA4_2
    """
    counts = {}
    unique = []

    for label in labels:
        label = str(label)
        counts[label] = counts.get(label, 0) + 1

        if counts[label] == 1:
            unique.append(label)
        else:
            unique.append(f"{label}_{counts[label]}")

    return unique

print("\n================================================")
print("ONCOGENEX - DIFFERENTIAL EXPRESSION PIPELINE")
print("================================================")

print("\n--- LOADING DATA ---")

expression = pd.read_csv(EXPRESSION_FILE, index_col=0)
metadata = pd.read_excel(METADATA_FILE)

ensure_required_columns(metadata, ["gsm_id", "group"], "Metadata")

expression.index = expression.index.astype(str).str.strip()
metadata["gsm_id"] = metadata["gsm_id"].astype(str).str.strip()

expression = expression.apply(pd.to_numeric, errors="coerce")

metadata["group"] = (
    metadata["group"]
    .astype(str)
    .str.strip()
    .str.lower()
)

print("Expression shape:", expression.shape)
print("Metadata shape:", metadata.shape)

print("\n--- ALIGNING SAMPLES ---")

metadata = metadata.drop_duplicates(subset="gsm_id", keep="first")
metadata = metadata.set_index("gsm_id")

common_ids = expression.index.intersection(metadata.index)

if len(common_ids) == 0:
    raise ValueError(
        "No matching sample IDs were found between expression data and metadata."
    )

expression = expression.loc[common_ids].copy()
metadata = metadata.loc[common_ids].copy()

if not expression.index.equals(metadata.index):
    raise RuntimeError("Sample alignment failed.")

print("Matching samples:", len(common_ids))
print("Expression samples after alignment:", expression.shape[0])
print("Metadata samples after alignment:", metadata.shape[0])
print("\nGroup distribution:")
print(metadata["group"].value_counts())

print("\n--- SEPARATING GROUPS ---")

required_groups = {"cancer", "normal"}
observed_groups = set(metadata["group"].dropna().unique())

if not required_groups.issubset(observed_groups):
    raise ValueError(
        f"Both 'cancer' and 'normal' groups are required.\n"
        f"Observed groups: {sorted(observed_groups)}"
    )

cancer_mask = metadata["group"].eq("cancer")
normal_mask = metadata["group"].eq("normal")

cancer_samples = expression.loc[cancer_mask].copy()
normal_samples = expression.loc[normal_mask].copy()

if len(cancer_samples) < 2 or len(normal_samples) < 2:
    raise ValueError(
        "At least two samples are required in each group for differential testing."
    )

print("Cancer samples:", len(cancer_samples))
print("Normal samples:", len(normal_samples))

print("\n--- FILTERING GENES ---")

nonempty_genes = expression.columns[
    expression.notna().any(axis=0)
]

expression = expression.loc[:, nonempty_genes]

gene_variance = expression.var(axis=0, ddof=1)
valid_genes = gene_variance[gene_variance > 1e-10].index

print("Total genes after removing empty genes:", expression.shape[1])
print("Genes retained for statistical testing:", len(valid_genes))

if len(valid_genes) == 0:
    raise ValueError("No genes remain after variance filtering.")

print("\n--- PERFORMING WELCH'S T-TEST ---")

cancer_valid = cancer_samples.loc[:, valid_genes]
normal_valid = normal_samples.loc[:, valid_genes]

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    t_statistics, p_values = ttest_ind(
        cancer_valid,
        normal_valid,
        axis=0,
        equal_var=False,
        nan_policy="omit"
    )

t_statistics = pd.Series(t_statistics, index=valid_genes, dtype=float)
p_values = pd.Series(p_values, index=valid_genes, dtype=float)

valid_pvalue_mask = np.isfinite(p_values)
print("Valid p-values:", int(valid_pvalue_mask.sum()))
print("Invalid p-values:", int((~valid_pvalue_mask).sum()))

print("\n--- PERFORMING FDR CORRECTION ---")

adjusted_p_values = pd.Series(
    np.nan,
    index=expression.columns,
    dtype=float
)

valid_p_values = p_values.loc[valid_pvalue_mask]

if len(valid_p_values) > 0:
    _, adjusted_valid, _, _ = multipletests(
        valid_p_values.values,
        method="fdr_bh"
    )
    adjusted_p_values.loc[valid_p_values.index] = adjusted_valid

print("\n--- CREATING RESULTS TABLE ---")

mean_cancer = cancer_samples.mean(axis=0)
mean_normal = normal_samples.mean(axis=0)

log2_fold_change = mean_cancer - mean_normal

all_t_statistics = pd.Series(
    np.nan,
    index=expression.columns,
    dtype=float
)
all_t_statistics.loc[t_statistics.index] = t_statistics

all_p_values = pd.Series(
    np.nan,
    index=expression.columns,
    dtype=float
)
all_p_values.loc[p_values.index] = p_values

de_results = pd.DataFrame({
    "Gene": expression.columns,
    "Mean_Cancer": mean_cancer.reindex(expression.columns).values,
    "Mean_Normal": mean_normal.reindex(expression.columns).values,
    "Log2FC": log2_fold_change.reindex(expression.columns).values,
    "T_Statistic": all_t_statistics.reindex(expression.columns).values,
    "P_Value": all_p_values.reindex(expression.columns).values,
    "Adjusted_P_Value": adjusted_p_values.reindex(expression.columns).values,
})

de_results["Significant"] = (
    (de_results["Adjusted_P_Value"] < FDR_THRESHOLD)
    & (de_results["Log2FC"].abs() >= LOG2FC_THRESHOLD)
)

de_results["Regulation"] = np.select(
    [
        de_results["Significant"] & (de_results["Log2FC"] > 0),
        de_results["Significant"] & (de_results["Log2FC"] < 0),
    ],
    [
        "Upregulated",
        "Downregulated",
    ],
    default="Not significant"
)

significant_genes = (
    de_results.loc[de_results["Significant"]]
    .sort_values(
        ["Adjusted_P_Value", "Log2FC"],
        ascending=[True, False]
    )
    .copy()
)

print("\n--- MAPPING PROBES TO GENE SYMBOLS ---")

annotation = pd.read_csv(
    ANNOTATION_FILE,
    sep="\t",
    comment="#",
    low_memory=False
)

ensure_required_columns(
    annotation,
    ["ID", "Gene Symbol", "Gene Title"],
    "GPL570 annotation file"
)

gene_annotation = annotation[
    ["ID", "Gene Symbol", "Gene Title"]
].copy()

gene_annotation = gene_annotation.rename(columns={"ID": "Gene"})
gene_annotation["Gene"] = gene_annotation["Gene"].astype(str).str.strip()
gene_annotation["Gene Symbol"] = gene_annotation["Gene Symbol"].apply(
    clean_gene_symbol
)

gene_annotation = gene_annotation.drop_duplicates(
    subset="Gene",
    keep="first"
)

de_results_annotated = de_results.merge(
    gene_annotation,
    on="Gene",
    how="left"
)

de_results_annotated["Gene Symbol"] = (
    de_results_annotated["Gene Symbol"]
    .fillna("")
    .astype(str)
    .str.strip()
)

de_results_annotated["Gene Title"] = (
    de_results_annotated["Gene Title"]
    .fillna("")
    .astype(str)
    .str.strip()
)

de_results_annotated["Display_Name"] = np.where(
    de_results_annotated["Gene Symbol"].ne(""),
    de_results_annotated["Gene Symbol"],
    de_results_annotated["Gene"]
)

significant_genes_annotated = (
    de_results_annotated.loc[
        de_results_annotated["Significant"]
    ]
    .copy()
)

print("Total significant probes:", len(significant_genes_annotated))
print(
    "Mapped significant probes:",
    int(significant_genes_annotated["Gene Symbol"].ne("").sum())
)

print("\n--- CREATING UNIQUE GENE LISTS ---")

mapped_significant = significant_genes_annotated.loc[
    significant_genes_annotated["Gene Symbol"].ne("")
].copy()

mapped_significant = (
    mapped_significant
    .sort_values(
        ["Adjusted_P_Value", "Log2FC"],
        ascending=[True, False]
    )
    .drop_duplicates(subset="Gene Symbol", keep="first")
)

upregulated_genes = (
    mapped_significant.loc[
        mapped_significant["Regulation"] == "Upregulated"
    ]
    .sort_values(
        ["Log2FC", "Adjusted_P_Value"],
        ascending=[False, True]
    )
    .copy()
)

downregulated_genes = (
    mapped_significant.loc[
        mapped_significant["Regulation"] == "Downregulated"
    ]
    .sort_values(
        ["Log2FC", "Adjusted_P_Value"],
        ascending=[True, True]
    )
    .copy()
)

top_20_upregulated = upregulated_genes.head(20)
top_20_downregulated = downregulated_genes.head(20)

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

de_results_annotated.to_csv(
    RESULTS_DIR / "all_differential_expression_results_annotated.csv",
    index=False
)

significant_genes_annotated.to_csv(
    RESULTS_DIR / "significant_genes_with_symbols.csv",
    index=False
)

upregulated_genes.to_csv(
    RESULTS_DIR / "upregulated_genes.csv",
    index=False
)

downregulated_genes.to_csv(
    RESULTS_DIR / "downregulated_genes.csv",
    index=False
)

top_20_upregulated.to_csv(
    RESULTS_DIR / "top_20_upregulated_genes.csv",
    index=False
)

top_20_downregulated.to_csv(
    RESULTS_DIR / "top_20_downregulated_genes.csv",
    index=False
)

print("\n--- CREATING VOLCANO PLOT ---")

volcano_data = de_results_annotated.dropna(
    subset=["Log2FC", "Adjusted_P_Value"]
).copy()

volcano_data["Plot_Adjusted_P"] = (
    volcano_data["Adjusted_P_Value"].clip(lower=1e-300)
)
volcano_data["Minus_Log10_P"] = -np.log10(
    volcano_data["Plot_Adjusted_P"]
)

plt.figure(figsize=(10, 8))

not_sig = volcano_data.loc[
    volcano_data["Regulation"] == "Not significant"
]
down = volcano_data.loc[
    volcano_data["Regulation"] == "Downregulated"
]
up = volcano_data.loc[
    volcano_data["Regulation"] == "Upregulated"
]

plt.scatter(
    not_sig["Log2FC"],
    not_sig["Minus_Log10_P"],
    alpha=0.25,
    s=14,
    label="Not significant"
)

plt.scatter(
    down["Log2FC"],
    down["Minus_Log10_P"],
    alpha=0.65,
    s=18,
    label="Downregulated"
)

plt.scatter(
    up["Log2FC"],
    up["Minus_Log10_P"],
    alpha=0.65,
    s=18,
    label="Upregulated"
)

plt.axvline(-LOG2FC_THRESHOLD, linestyle="--")
plt.axvline(LOG2FC_THRESHOLD, linestyle="--")
plt.axhline(-np.log10(FDR_THRESHOLD), linestyle="--")

plt.xlabel("Log2 Fold Change")
plt.ylabel("-Log10 Adjusted P-Value")
plt.title("Volcano Plot of All Tested Genes")
plt.legend()
plt.tight_layout()
plt.savefig(
    PLOTS_DIR / "volcano_plot.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

print("Volcano plot saved.")

print("\n--- CREATING MA PLOT ---")

ma_data = de_results_annotated.dropna(
    subset=["Mean_Cancer", "Mean_Normal", "Log2FC"]
).copy()

ma_data["Average_Expression"] = (
    ma_data["Mean_Cancer"] + ma_data["Mean_Normal"]
) / 2

plt.figure(figsize=(10, 7))

ma_not_sig = ma_data.loc[
    ma_data["Regulation"] == "Not significant"
]
ma_down = ma_data.loc[
    ma_data["Regulation"] == "Downregulated"
]
ma_up = ma_data.loc[
    ma_data["Regulation"] == "Upregulated"
]

plt.scatter(
    ma_not_sig["Average_Expression"],
    ma_not_sig["Log2FC"],
    alpha=0.20,
    s=10,
    label="Not significant"
)

plt.scatter(
    ma_down["Average_Expression"],
    ma_down["Log2FC"],
    alpha=0.65,
    s=16,
    label="Downregulated"
)

plt.scatter(
    ma_up["Average_Expression"],
    ma_up["Log2FC"],
    alpha=0.65,
    s=16,
    label="Upregulated"
)

plt.axhline(0, linestyle="--")
plt.axhline(-LOG2FC_THRESHOLD, linestyle="--", alpha=0.6)
plt.axhline(LOG2FC_THRESHOLD, linestyle="--", alpha=0.6)

plt.xlabel("Average Log2 Expression")
plt.ylabel("Log2 Fold Change")
plt.title("MA Plot of All Tested Genes")
plt.legend()
plt.tight_layout()
plt.savefig(
    PLOTS_DIR / "ma_plot.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

print("MA plot saved.")

print("\n--- CREATING TOP 10 UPREGULATED GENE PLOT ---")

plot_up = upregulated_genes.head(10).sort_values(
    "Log2FC",
    ascending=True
)

plt.figure(figsize=(10, 7))
plt.barh(
    plot_up["Gene Symbol"],
    plot_up["Log2FC"]
)
plt.xlabel("Log2 Fold Change")
plt.ylabel("Gene Symbol")
plt.title("Top 10 Upregulated Genes in Lung Adenocarcinoma")
plt.tight_layout()
plt.savefig(
    PLOTS_DIR / "top_10_upregulated_genes.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

print("Top upregulated gene plot saved.")

print("\n--- CREATING TOP 10 DOWNREGULATED GENE PLOT ---")

plot_down = downregulated_genes.head(10).sort_values(
    "Log2FC",
    ascending=False
)

plt.figure(figsize=(10, 7))
plt.barh(
    plot_down["Gene Symbol"],
    plot_down["Log2FC"]
)
plt.xlabel("Log2 Fold Change")
plt.ylabel("Gene Symbol")
plt.title("Top 10 Downregulated Genes in Lung Adenocarcinoma")
plt.tight_layout()
plt.savefig(
    PLOTS_DIR / "top_10_downregulated_genes.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

print("Top downregulated gene plot saved.")

print("\n--- CREATING HEATMAP ---")

heatmap_candidates = (
    mapped_significant
    .sort_values(
        ["Adjusted_P_Value", "Log2FC"],
        ascending=[True, False]
    )
    .head(HEATMAP_TOP_N)
    .copy()
)

top_probe_ids = [
    probe
    for probe in heatmap_candidates["Gene"]
    if probe in expression.columns
]

if len(top_probe_ids) < 2:
    print("Heatmap skipped: fewer than 2 selected probes were found.")
else:

    heatmap_matrix = expression.loc[:, top_probe_ids].T.copy()

    probe_to_symbol = (
        heatmap_candidates
        .set_index("Gene")["Gene Symbol"]
        .to_dict()
    )

    heatmap_matrix.index = make_unique_labels([
        probe_to_symbol.get(probe, probe)
        for probe in heatmap_matrix.index
    ])

    heatmap_scaled = zscore_rows(heatmap_matrix)

    sample_groups = metadata.loc[
        heatmap_scaled.columns,
        "group"
    ]

    group_palette = {
        "cancer": "red",
        "normal": "blue"
    }

    col_colors = sample_groups.map(group_palette)

    g = sns.clustermap(
        heatmap_scaled,
        cmap="coolwarm",
        center=0,
        col_colors=col_colors,
        xticklabels=False,
        yticklabels=True,
        figsize=(16, 10),
        method="average",
        metric="euclidean"
    )

    g.fig.suptitle(
        "Clustered Heatmap of Top 20 Differentially Expressed Genes",
        y=1.02
    )

    g.savefig(
        PLOTS_DIR / "top_20_deg_heatmap.png",
        dpi=300,
        bbox_inches="tight"
    )
    plt.close(g.fig)

    print("Heatmap saved.")

print("\n--- CREATING PCA PLOT ---")

pca_variance = expression.var(axis=0, ddof=1)
pca_genes = (
    pca_variance
    .sort_values(ascending=False)
    .head(min(PCA_TOP_VARIABLE_GENES, expression.shape[1]))
    .index
)

pca_input = expression.loc[:, pca_genes].copy()

pca_input = pca_input.apply(
    lambda col: col.fillna(col.median()),
    axis=0
)

pca_input = pca_input.dropna(axis=1, how="all")

X_scaled = StandardScaler().fit_transform(pca_input)

pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

pca_df = pd.DataFrame(
    X_pca,
    index=expression.index,
    columns=["PC1", "PC2"]
)

pca_df["Group"] = metadata.loc[
    pca_df.index,
    "group"
]

plt.figure(figsize=(10, 7))

for group in ["normal", "cancer"]:
    subset = pca_df.loc[pca_df["Group"] == group]

    if len(subset) == 0:
        continue

    plt.scatter(
        subset["PC1"],
        subset["PC2"],
        label=group,
        alpha=0.8,
        s=60
    )

plt.xlabel(
    f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}% variance)"
)
plt.ylabel(
    f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}% variance)"
)
plt.title("PCA of Lung Adenocarcinoma and Normal Samples")
plt.legend(title="Group")
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(
    PLOTS_DIR / "pca_plot.png",
    dpi=300,
    bbox_inches="tight"
)
plt.close()

print("PCA plot saved.")

print("\n================================================")
print("ANALYSIS COMPLETED SUCCESSFULLY")
print("================================================")

print("\nSummary:")
print("Total tested genes:", int(de_results["P_Value"].notna().sum()))
print("Significant genes:", len(significant_genes_annotated))
print("Unique upregulated genes:", len(upregulated_genes))
print("Unique downregulated genes:", len(downregulated_genes))

print("\nResult files:")
print(RESULTS_DIR / "all_differential_expression_results_annotated.csv")
print(RESULTS_DIR / "significant_genes_with_symbols.csv")
print(RESULTS_DIR / "upregulated_genes.csv")
print(RESULTS_DIR / "downregulated_genes.csv")
print(RESULTS_DIR / "top_20_upregulated_genes.csv")
print(RESULTS_DIR / "top_20_downregulated_genes.csv")

print("\nPlot files:")
print(PLOTS_DIR / "volcano_plot.png")
print(PLOTS_DIR / "ma_plot.png")
print(PLOTS_DIR / "top_10_upregulated_genes.png")
print(PLOTS_DIR / "top_10_downregulated_genes.png")
print(PLOTS_DIR / "top_20_deg_heatmap.png")
print(PLOTS_DIR / "pca_plot.png")
