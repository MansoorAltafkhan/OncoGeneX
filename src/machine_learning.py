import pandas as pd
import numpy as np
import os
os.makedirs("results", exist_ok=True)

from sklearn.metrics import roc_auc_score, roc_curve
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix,
    roc_auc_score
)

print("=== MACHINE LEARNING ANALYSIS ===")

# Load processed expression data
expression = pd.read_csv(
    "data/processed/processed_expression.csv",
    index_col=0
)

# Load metadata
metadata = pd.read_excel(
    "data/raw/GSE40791_final_metadata.xlsx"
)

# Set GSM ID as index
metadata = metadata.set_index("gsm_id")

# Align metadata with expression samples
metadata = metadata.loc[expression.index]
print("\nExpression shape:", expression.shape)
print("Metadata shape:", metadata.shape)

# Create target variable
y = metadata["group"]

# Features
X = expression

print("\nClass distribution:")
print(y.value_counts())

print("\nFeatures shape:", X.shape)
import pandas as pd

from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    classification_report,
    confusion_matrix,
    roc_auc_score
)

print("\n=== MACHINE LEARNING ANALYSIS ===")

# -------------------------------
# X = gene expression features
# y = cancer/normal labels
# -------------------------------

# Make sure these variables correspond to your aligned datasets
X = expression
y = metadata["group"]

print("\nFeature matrix shape:", X.shape)
print("Target distribution:")
print(y.value_counts())


# -------------------------------
# Train-test split
# -------------------------------

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.20,
    random_state=42,
    stratify=y
)

print("\nTraining samples:", X_train.shape[0])
print("Testing samples:", X_test.shape[0])


# -------------------------------
# Machine Learning Pipeline
# -------------------------------

model = Pipeline([
    ("feature_selection", SelectKBest(score_func=f_classif, k=1000)),
    ("scaler", StandardScaler()),
    ("classifier", LogisticRegression(
        max_iter=5000,
        random_state=42
    ))
])


# -------------------------------
# Train model
# -------------------------------

print("\nTraining model...")
model.fit(X_train, y_train)


# -------------------------------
# Predictions
# -------------------------------

y_pred = model.predict(X_test)
# ----------------------------------------
# ROC-AUC ANALYSIS
# ----------------------------------------

# Get probability predictions
y_prob = model.predict_proba(X_test)

print("\nClasses:", model.classes_)

# Explicitly select cancer probability
cancer_index = list(model.classes_).index("cancer")
y_prob_cancer = y_prob[:, cancer_index]

# Convert labels: cancer = 1, normal = 0
y_test_binary = (y_test == "cancer").astype(int)

# Calculate ROC-AUC
roc_auc = roc_auc_score(y_test_binary, y_prob_cancer)

print("\n=== ROC-AUC RESULTS ===")
print("Test ROC-AUC:", roc_auc)

print("\n=== MODEL RESULTS ===")

print("\nAccuracy:")
print(accuracy_score(y_test, y_pred))

print("\nClassification Report:")
print(classification_report(y_test, y_pred))

print("\nConfusion Matrix:")
print(confusion_matrix(y_test, y_pred))
# ----------------------------------------
# Cross-validation
# ----------------------------------------

cv = StratifiedKFold(
    n_splits=5,
    shuffle=True,
    random_state=42
)

cv_scores = cross_val_score(
    model,
    X,
    y,
    cv=cv,
    scoring="accuracy"
)

print("\n=== CROSS-VALIDATION RESULTS ===")
print("CV Scores:", cv_scores)
print("Mean CV Accuracy:", cv_scores.mean())
print("Standard Deviation:", cv_scores.std())
# ----------------------------------------
# CROSS-VALIDATED ROC-AUC
# ----------------------------------------

cv_auc_scores = cross_val_score(
    model,
    X,
    y,
    cv=cv,
    scoring="roc_auc"
)

print("\n=== CROSS-VALIDATED ROC-AUC RESULTS ===")
print("CV ROC-AUC Scores:", cv_auc_scores)
print("Mean CV ROC-AUC:", cv_auc_scores.mean())
print("Standard Deviation:", cv_auc_scores.std())
# ==========================================
# TOP 100 SELECTED GENE ANALYSIS
# ==========================================

print("\n=== TOP 100 SELECTED GENES ===")

# Fit SelectKBest on the complete dataset
selector = SelectKBest(score_func=f_classif, k=100)
selector.fit(X, y)

# Get selected gene names
selected_features = X.columns[selector.get_support()]

# Get F-scores
feature_scores = selector.scores_[selector.get_support()]

# Create a DataFrame
top_genes = pd.DataFrame({
    'Gene': selected_features,
    'F_Score': feature_scores
})

# Sort genes by importance
top_genes = top_genes.sort_values(
    by='F_Score',
    ascending=False
)

print(top_genes)
top_genes.to_csv(
    "results/top_100_selected_genes.csv",
    index=False
)

print("\nTop 100 genes saved successfully!")
# ==========================================
# CONVERT PROBE IDs TO GENE SYMBOLS
# ==========================================

annotation = pd.read_csv(
    "data/raw/GPL570-55999.txt",
    sep="\t",
    comment="#",
    low_memory=False
)

print(annotation.columns)
# ==========================================
# MAP PROBE IDs TO GENE SYMBOLS
# ==========================================

print("\n=== CONVERTING PROBE IDs TO GENE SYMBOLS ===")

# Keep only required columns
gene_annotation = annotation[
    ["ID", "Gene Symbol", "Gene Title"]
].copy()

# Rename ID column to match the Gene column
gene_annotation = gene_annotation.rename(
    columns={"ID": "Gene"}
)

# Merge selected genes with annotation data
top_genes_annotated = top_genes.merge(
    gene_annotation,
    on="Gene",
    how="left"
)

# Display results
print("\nTop genes with gene symbols:")
print(top_genes_annotated.head(20))

# Save annotated results
top_genes_annotated.to_csv(
    "results/top_100_genes_with_symbols.csv",
    index=False
)

print("\nGene symbol mapping completed successfully!")
print("File saved: results/top_100_genes_with_symbols.csv")
# ==========================================
# UPDATE GENE NAMES FOR PLOTS
# ==========================================

# Merge top genes with annotation data
top_genes_with_symbols = top_genes.merge(
    annotation[["ID", "Gene Symbol", "Gene Title"]],
    left_on="Gene",
    right_on="ID",
    how="left"
)

# Clean gene symbols
top_genes_with_symbols["Gene Symbol"] = (
    top_genes_with_symbols["Gene Symbol"]
    .fillna(top_genes_with_symbols["Gene"])
    .astype(str)
    .str.split(" /// ")
    .str[0]
)

# Use gene symbols as display names
top_genes_with_symbols["Display_Name"] = (
    top_genes_with_symbols["Gene Symbol"]
)

# Save updated file
top_genes_with_symbols.to_csv(
    "results/top_100_genes_with_symbols.csv",
    index=False
)

print("\nGene symbols added successfully!")
print(top_genes_with_symbols[[
    "Gene",
    "Gene Symbol",
    "Gene Title",
    "F_Score"
]].head(20))

# ==========================================
# CLUSTERED HEATMAP WITH GENE SYMBOLS
# ==========================================

import seaborn as sns
import matplotlib.pyplot as plt

print("\n=== CREATING CLUSTERED HEATMAP ===")

# Get top 20 genes with annotation information
top_20_genes = top_genes_with_symbols.head(20).copy()

# Probe IDs are needed to extract expression data
top_20_probe_ids = top_20_genes["Gene"].tolist()

# Gene symbols will be used as labels in the heatmap
top_20_display_names = top_20_genes["Display_Name"].tolist()

# Extract expression data using probe IDs
cluster_data = X[top_20_probe_ids].copy()

# Rename probe ID columns to gene symbols
cluster_data.columns = top_20_display_names

# Standardize each gene
cluster_data = cluster_data.apply(
    lambda x: (x - x.mean()) / x.std(),
    axis=0
)

# Create colours for cancer and normal samples
class_colors = y.map({
    "cancer": "red",
    "normal": "blue"
})

# Create clustered heatmap
g = sns.clustermap(
    cluster_data.T,
    cmap="coolwarm",
    center=0,
    col_colors=class_colors,
    xticklabels=False,
    yticklabels=True,
    figsize=(16, 10),
    method="average",
    metric="euclidean"
)

# Add title
g.fig.suptitle(
    "Clustered Heatmap of Top 20 Selected Genes",
    y=1.02
)

# Save figure
g.savefig(
    "results/clustered_top_20_gene_heatmap.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("\nClustered heatmap with gene symbols saved successfully!")
# ============================================
# PCA VISUALIZATION
# ============================================

from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

print("\n=== PCA VISUALIZATION ===")

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

# Standardize expression data
X_scaled = StandardScaler().fit_transform(X)

# Perform PCA
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X_scaled)

# Create dataframe for plotting
pca_df = pd.DataFrame(
    X_pca,
    columns=["PC1", "PC2"]
)

pca_df["Group"] = y.values

# Create scatter plot
plt.figure(figsize=(10, 7))

for group in pca_df["Group"].unique():
    
    subset = pca_df[pca_df["Group"] == group]
    
    plt.scatter(
        subset["PC1"],
        subset["PC2"],
        label=group,
        alpha=0.8,
        s=60
    )

# Add explained variance to axis labels
plt.xlabel(
    f"PC1 ({pca.explained_variance_ratio_[0] * 100:.2f}% variance)"
)

plt.ylabel(
    f"PC2 ({pca.explained_variance_ratio_[1] * 100:.2f}% variance)"
)

plt.title("PCA of Cancer and Normal Samples")

plt.legend(title="Group")

plt.grid(True, alpha=0.3)

# Save figure
plt.savefig(
    "results/pca_visualization.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("PCA visualization saved successfully!")
# ============================================
# TOP 20 FEATURE IMPORTANCE PLOT
# ============================================

import matplotlib.pyplot as plt

print("\n=== TOP 20 FEATURE IMPORTANCE ===")

# Select top 20 features with gene symbols
top_20_importance = top_genes_with_symbols.head(20).copy()

# Sort so the most important feature appears at the top
top_20_importance = top_20_importance.sort_values(
    by="F_Score",
    ascending=True
)

# Create horizontal bar plot
plt.figure(figsize=(10, 8))

plt.barh(
    top_20_importance["Gene Symbol"],
    top_20_importance["F_Score"]
)

plt.xlabel("F-Score")
plt.ylabel("Gene Symbol")
plt.title("Top 20 Most Important Selected Genes")

plt.tight_layout()

# Save figure
plt.savefig(
    "results/top_20_feature_importance.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Top 20 feature importance plot saved successfully!")
# ============================================
# ROC CURVE VISUALIZATION
# ============================================

import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc

print("\n=== ROC CURVE VISUALIZATION ===")

# Calculate ROC curve
fpr, tpr, thresholds = roc_curve(
    y_test_binary,
    y_prob_cancer
)

# Calculate AUC
roc_auc = auc(fpr, tpr)

# Create plot
plt.figure(figsize=(8, 6))

plt.plot(
    fpr,
    tpr,
    label=f"ROC Curve (AUC = {roc_auc:.3f})"
)

# Diagonal reference line
plt.plot(
    [0, 1],
    [0, 1],
    linestyle="--",
    label="Random Classifier"
)

plt.xlabel("False Positive Rate")
plt.ylabel("True Positive Rate")
plt.title("ROC Curve for Cancer Classification")
plt.legend(loc="lower right")
plt.grid(True, alpha=0.3)

plt.tight_layout()

# Save figure
plt.savefig(
    "results/roc_curve.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("ROC curve saved successfully!")
# ============================================
# CONFUSION MATRIX HEATMAP
# ============================================

import seaborn as sns
import matplotlib.pyplot as plt

print("\n=== CONFUSION MATRIX VISUALIZATION ===")

# Calculate confusion matrix
cm = confusion_matrix(
    y_test,
    y_pred,
    labels=["cancer", "normal"]
)

# Create plot
plt.figure(figsize=(7, 6))

sns.heatmap(
    cm,
    annot=True,
    fmt="d",
    cmap="Blues",
    xticklabels=["Cancer", "Normal"],
    yticklabels=["Cancer", "Normal"]
)

plt.xlabel("Predicted Label")
plt.ylabel("True Label")
plt.title("Confusion Matrix for Cancer Classification")

plt.tight_layout()

# Save figure
plt.savefig(
    "results/confusion_matrix.png",
    dpi=300,
    bbox_inches="tight"
)

plt.show()

print("Confusion matrix heatmap saved successfully!")