import os
import math
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt

os.makedirs("results/enrichment", exist_ok=True)

up_file = "results/upregulated_genes.csv"
down_file = "results/downregulated_genes.csv"

up_df = pd.read_csv(up_file)
down_df = pd.read_csv(down_file)

print("\nUpregulated file columns:")
print(up_df.columns.tolist())

print("\nDownregulated file columns:")
print(down_df.columns.tolist())

def extract_gene_symbols(df):

    possible_columns = [
        "Gene Symbol",
        "gene_symbol",
        "GeneSymbol",
        "SYMBOL",
        "symbol"
    ]

    gene_column = None

    for col in possible_columns:
        if col in df.columns:
            gene_column = col
            break

    if gene_column is None:
        raise ValueError(
            f"Gene symbol column not found. "
            f"Available columns: {df.columns.tolist()}"
        )

    genes = (
        df[gene_column]
        .dropna()
        .astype(str)
        .str.strip()
    )

    genes = genes[genes != ""]
    genes = genes.drop_duplicates().tolist()

    return genes

up_genes = extract_gene_symbols(up_df)
down_genes = extract_gene_symbols(down_df)

print("\nNumber of upregulated genes:", len(up_genes))
print("Number of downregulated genes:", len(down_genes))

gene_sets = [
    "GO_Biological_Process_2023",
    "KEGG_2021_Human"
]

def run_enrichment(genes, group_name):

    print(f"\nRunning enrichment analysis for {group_name} genes...")

    enr = gp.enrichr(
        gene_list=genes,
        gene_sets=gene_sets,
        organism="human",
        outdir=None
    )

    results = enr.results

    if results is None or results.empty:
        print(
            f"No enrichment results found for {group_name}"
        )
        return None

    output_file = (
        f"results/enrichment/"
        f"{group_name}_enrichment.csv"
    )

    results.to_csv(
        output_file,
        index=False
    )

    print(f"Results saved to: {output_file}")

    return results

up_results = run_enrichment(
    up_genes,
    "upregulated"
)

down_results = run_enrichment(
    down_genes,
    "downregulated"
)

def save_separate_results(results, group_name):

    if results is None or results.empty:
        return None, None

    go_results = results[
        results["Gene_set"].str.contains(
            "GO_Biological_Process",
            case=False,
            na=False
        )
    ].copy()

    kegg_results = results[
        results["Gene_set"].str.contains(
            "KEGG",
            case=False,
            na=False
        )
    ].copy()

    go_file = (
        f"results/enrichment/"
        f"GO_{group_name}.csv"
    )

    go_results.to_csv(
        go_file,
        index=False
    )

    kegg_file = (
        f"results/enrichment/"
        f"KEGG_{group_name}.csv"
    )

    kegg_results.to_csv(
        kegg_file,
        index=False
    )

    print(f"Saved: {go_file}")
    print(f"Saved: {kegg_file}")

    return go_results, kegg_results

go_up, kegg_up = save_separate_results(
    up_results,
    "upregulated"
)

go_down, kegg_down = save_separate_results(
    down_results,
    "downregulated"
)

def create_barplot(results, analysis_type, group_name):

    if results is None or results.empty:
        print(
            f"No results available for "
            f"{analysis_type} {group_name} plot."
        )
        return

    top_results = results.sort_values(
        "Adjusted P-value",
        ascending=True
    ).head(10).copy()

    top_results["Adjusted P-value"] = (
        top_results["Adjusted P-value"]
        .clip(lower=1e-300)
    )

    top_results["minus_log10_p"] = (
        -top_results["Adjusted P-value"].apply(
            math.log10
        )
    )

    top_results = top_results.sort_values(
        "minus_log10_p",
        ascending=True
    )

    plt.figure(figsize=(12, 7))

    plt.barh(
        top_results["Term"],
        top_results["minus_log10_p"]
    )

    plt.xlabel(
        "-Log10 Adjusted P-value",
        fontsize=12
    )

    plt.ylabel(
        "Pathway / Biological Process",
        fontsize=12
    )

    plt.title(
        f"Top 10 {analysis_type} Enriched Terms: "
        f"{group_name.capitalize()} Genes",
        fontsize=15
    )

    plt.tight_layout()

    analysis_filename = analysis_type.lower().replace(" ", "_")

    output_file = (
        f"results/enrichment/"
        f"{analysis_filename}_{group_name}_barplot.png"
    )

    plt.savefig(
        output_file,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()

    print(f"Plot saved: {output_file}")

create_barplot(
    go_up,
    "GO Biological Process",
    "upregulated"
)

create_barplot(
    go_down,
    "GO Biological Process",
    "downregulated"
)

create_barplot(
    kegg_up,
    "KEGG Pathway",
    "upregulated"
)

create_barplot(
    kegg_down,
    "KEGG Pathway",
    "downregulated"
)

print("\n=================================================")
print("PATHWAY ENRICHMENT ANALYSIS COMPLETE")
print("=================================================")

print("\nResults saved in:")
print("results/enrichment/")

print("\nGenerated plots:")
print("1. GO Biological Process - Upregulated")
print("2. GO Biological Process - Downregulated")
print("3. KEGG Pathway - Upregulated")
print("4. KEGG Pathway - Downregulated")
