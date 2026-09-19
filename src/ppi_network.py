import os
import pandas as pd
import requests
import networkx as nx
import matplotlib.pyplot as plt
import time

os.makedirs("results/ppi", exist_ok=True)

up_file = "results/upregulated_genes.csv"
down_file = "results/downregulated_genes.csv"

up_df = pd.read_csv(up_file)
down_df = pd.read_csv(down_file)

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
        .drop_duplicates()
        .tolist()
    )

    return genes

up_genes = extract_gene_symbols(up_df)
down_genes = extract_gene_symbols(down_df)

print("Number of upregulated genes:", len(up_genes))
print("Number of downregulated genes:", len(down_genes))

def get_string_network(genes, group_name):

    print(f"\nFetching PPI network for {group_name} genes...")

    url = "https://string-db.org/api/tsv/network"

    batch_size = 400

    all_networks = []

    total_batches = (
        len(genes) + batch_size - 1
    ) // batch_size

    for i in range(0, len(genes), batch_size):

        batch_number = (
            i // batch_size
        ) + 1

        gene_batch = genes[
            i:i + batch_size
        ]

        print(
            f"Processing batch "
            f"{batch_number}/{total_batches} "
            f"({len(gene_batch)} genes)..."
        )

        identifiers = "\r".join(
            gene_batch
        )

        params = {
            "identifiers": identifiers,
            "species": 9606,
            "required_score": 700
        }

        try:

            response = requests.post(
                url,
                data=params,
                timeout=60
            )

            if response.status_code != 200:

                print(
                    f"Batch {batch_number} failed: "
                    f"HTTP {response.status_code}"
                )

                print(
                    response.text[:500]
                )

                continue

            if not response.text.strip():

                print(
                    f"No interactions found "
                    f"for batch {batch_number}"
                )

                continue

            from io import StringIO

            batch_df = pd.read_csv(
                StringIO(response.text),
                sep="\t"
            )

            if not batch_df.empty:

                all_networks.append(
                    batch_df
                )

                print(
                    f"Interactions found in "
                    f"batch {batch_number}: "
                    f"{len(batch_df)}"
                )

        except Exception as e:

            print(
                f"Error in batch "
                f"{batch_number}: {e}"
            )

        time.sleep(1)

    if len(all_networks) == 0:

        print(
            f"No PPI network data found "
            f"for {group_name}"
        )

        return None

    network_df = pd.concat(
        all_networks,
        ignore_index=True
    )

    if (
        "preferredName_A" in network_df.columns
        and
        "preferredName_B" in network_df.columns
    ):

        network_df = network_df.drop_duplicates(
            subset=[
                "preferredName_A",
                "preferredName_B"
            ]
        )

    output_file = (
        f"results/ppi/"
        f"{group_name}_ppi_network.csv"
    )

    network_df.to_csv(
        output_file,
        index=False
    )

    print(
        f"\nTotal PPI interactions found: "
        f"{len(network_df)}"
    )

    print(
        f"Network saved to: "
        f"{output_file}"
    )

    return network_df

def create_network(network_df, group_name):

    if (
        network_df is None
        or
        network_df.empty
    ):

        print(
            f"No PPI network data available "
            f"for {group_name}"
        )

        return None

    G = nx.Graph()

    for _, row in network_df.iterrows():

        protein1 = row[
            "preferredName_A"
        ]

        protein2 = row[
            "preferredName_B"
        ]

        score = row["score"]

        G.add_edge(

            protein1,
            protein2,

            weight=score

        )

    print(
        f"\n{group_name.capitalize()} "
        f"PPI Network:"
    )

    print(
        "Nodes:",
        G.number_of_nodes()
    )

    print(
        "Edges:",
        G.number_of_edges()
    )

    return G

def identify_hub_genes(G, group_name):

    if G is None:

        return None

    print(
        f"\nIdentifying hub genes for "
        f"{group_name} genes..."
    )

    degree_centrality = (
        nx.degree_centrality(G)
    )

    betweenness_centrality = (
        nx.betweenness_centrality(G)
    )

    closeness_centrality = (
        nx.closeness_centrality(G)
    )

    hub_data = []

    for gene in G.nodes():

        hub_data.append({

            "Gene":
                gene,

            "Degree":
                G.degree(gene),

            "Degree_Centrality":
                degree_centrality[gene],

            "Betweenness_Centrality":
                betweenness_centrality[gene],

            "Closeness_Centrality":
                closeness_centrality[gene]

        })

    hub_df = pd.DataFrame(
        hub_data
    )

    hub_df = hub_df.sort_values(

        "Degree",

        ascending=False

    )

    output_file = (

        f"results/ppi/"
        f"{group_name}_hub_genes.csv"

    )

    hub_df.to_csv(

        output_file,

        index=False

    )

    print(

        f"Hub genes saved to: "
        f"{output_file}"

    )

    print(

        f"\nTop 10 hub genes "
        f"({group_name}):"

    )

    print(

        hub_df[
            [
                "Gene",
                "Degree",
                "Degree_Centrality"
            ]
        ].head(10)

    )

    return hub_df

def plot_ppi_network(
    G,
    hub_df,
    group_name
):

    if (
        G is None
        or
        hub_df is None
    ):

        return

    top_genes = (

        hub_df["Gene"]

        .head(20)

        .tolist()

    )

    subgraph = G.subgraph(
        top_genes
    )

    if (
        subgraph.number_of_nodes()
        == 0
    ):

        print(
            "No network available to plot."
        )

        return

    plt.figure(
        figsize=(12, 10)
    )

    pos = nx.spring_layout(

        subgraph,

        seed=42

    )

    node_sizes = [

        subgraph.degree(node) * 100

        for node in subgraph.nodes()

    ]

    nx.draw_networkx_nodes(

        subgraph,

        pos,

        node_size=node_sizes,

        alpha=0.8

    )

    nx.draw_networkx_edges(

        subgraph,

        pos,

        alpha=0.5

    )

    nx.draw_networkx_labels(

        subgraph,

        pos,

        font_size=8

    )

    plt.title(

        f"Top Hub Genes PPI Network: "
        f"{group_name.capitalize()}"

    )

    plt.axis("off")

    plt.tight_layout()

    output_file = (

        f"results/ppi/"
        f"{group_name}_ppi_network.png"

    )

    plt.savefig(

        output_file,

        dpi=300,

        bbox_inches="tight"

    )

    plt.close()

    print(

        f"PPI network plot saved: "
        f"{output_file}"

    )

up_network_df = get_string_network(

    up_genes,

    "upregulated"

)

up_graph = create_network(

    up_network_df,

    "upregulated"

)

up_hub_df = identify_hub_genes(

    up_graph,

    "upregulated"

)

plot_ppi_network(

    up_graph,

    up_hub_df,

    "upregulated"

)

down_network_df = get_string_network(

    down_genes,

    "downregulated"

)

down_graph = create_network(

    down_network_df,

    "downregulated"

)

down_hub_df = identify_hub_genes(

    down_graph,

    "downregulated"

)

plot_ppi_network(

    down_graph,

    down_hub_df,

    "downregulated"

)

print("\n=================================================")

print(
    "PPI NETWORK ANALYSIS COMPLETE"
)

print(
    "================================================="
)

print(
    "\nCheck the folder:"
)

print(
    "results/ppi/"
)
