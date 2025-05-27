---
jupytext:
  formats: md:myst
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.11.5
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Parametrized co-occurrence networks

:::{note} Last update 👈
:class: dropdown
David Palecek, May 27, 2025
:::

Complex networks provide metrics and to identify kestone species bridges/bottleneck species within the community. EMO-BON data are unified in methods but highly diverse in respect to the metadata (sampling station, temperature, season etc.), which can be interpreted as control/treatment. To refrase a particular question in respect to seasonality, is there structural difference in taxonomy buildup in the communities between spring, summer, autumn and winter.

We will construct association matrices between taxa over the series of sampling events grouped by the factor (season in this case).

## Recipe

1. Load data
2. Remove high taxa (non-identified sequences)
3. Pivot table from sampling events to taxonomies.
4. Remove low abundance taxa
5. Rarefy, or normalize
6. Remove replicates
7. Split to groups on chosen `factor`
8. Calculate associations (Bray-curits dissimilarity, Spearman's correlation, etc.)
9. Build network per group
10. Identify positive and negative associations
11. Downstream analysis

## Loading data and metadata

```{code-cell}
:label: import
:class: dropdown
import os
from typing import Dict
import numpy as np
import pandas as pd

from skbio.stats import subsample_counts
from skbio.diversity import beta_diversity
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
import networkx as nx

from momics.taxonomy import (
    pivot_taxonomic_data,
    separate_taxonomy)
from skbio.diversity import beta_diversity

from mgo.udal import UDAL

# All low level functions are imported from the momics package
from momics.loader import load_parquets_udal
from momics.metadata import get_metadata_udal, enhance_metadata

import seaborn as sns
import matplotlib.pyplot as plt
```

All the methods (eventually will be refactored to the marine-momics-methods) are for now defined here

```{code-cell}
:label: methods
:class: dropdown

def get_data():
    return load_parquets_udal()


# Load and merge metadata
def get_full_metadata():
    return get_metadata_udal()


def get_valid_samples():
    url = "https://raw.githubusercontent.com/emo-bon/momics-demos/refs/heads/main/data/shipment_b1b2_181.csv"
    df_valid = pd.read_csv(
        url, header=0, index_col=0,
    )
    return df_valid


def remove_high_taxa(df: pd.DataFrame, tax_level: str = 'phylum') -> pd.DataFrame:
    """
    Remove high level taxa from the dataframe.

    Args:
        df (pd.DataFrame): DataFrame containing taxonomic data.
        tax_level (str): The taxonomic level to filter by (e.g., 'phylum', 'class', 'order', etc.).

    Returns:
        pd.DataFrame: DataFrame with rows where the specified taxonomic level is not None.
    """
    if tax_level not in df.columns:
        raise ValueError(f"Taxonomic level '{tax_level}' not found in DataFrame.")
    
    # Filter out rows where the taxonomic level is None or NaN
    return df[~df[tax_level].isna()].copy()


def pivot_taxonomic_data(df: pd.DataFrame, normalize: str = None, rarefy_depth: int = None) -> pd.DataFrame:
    """
    Prepares the taxonomic data (LSU and SSU tables) for analysis. Apart from
    pivoting.
    
    Normalization of the pivot is optional. Methods include:
    - None: no normalization.
    - 'tss_sqrt': Total Sum Scaling and Square Root Transformation.
    - 'rarefy': rarefaction to a specified depth, if None, min of sample sums is used.

    TODO: refactor scaling to a new method and offer different options.

    Args:
        df (pd.DataFrame): The input DataFrame containing taxonomic information.
        normalize (str, optional): Normalization method.
            Options: None, 'tss_sqrt', 'rarefy'. Defaults to None.
        rarefy_depth (int, optional): Depth for rarefaction. If None, uses min sample sum.
            Defaults to None.

    Returns:
        pd.DataFrame: A pivot table with taxonomic data.
    """
    # Select relevant columns
    df['taxonomic_concat'] = (
        df['ncbi_tax_id'].astype(str) + 
        ';sk_' + df['superkingdom'].fillna('') + 
        ';k_' + df['kingdom'].fillna('') + 
        ';p_' + df['phylum'].fillna('') + 
        ';c_' + df['class'].fillna('') + 
        ';o_' + df['order'].fillna('') + 
        ';f_' + df['family'].fillna('') + 
        ';g_' + df['genus'].fillna('') + 
        ';s_' + df['species'].fillna('')
    )
    pivot_table = df.pivot_table(
        index=['ncbi_tax_id','taxonomic_concat'], 
        columns='ref_code', 
        values='abundance',
    ).fillna(0).astype(int)
    pivot_table = pivot_table.reset_index()
    # change inex name
    pivot_table.columns.name = None

    # normalize values
    if normalize == 'tss_sqrt':
        # Total Sum Scaling and Square Root Transformation
        pivot_table.iloc[:, 2:] = pivot_table.iloc[:, 2:].apply(lambda x: x / x.sum())
        pivot_table.iloc[:, 2:] = pivot_table.iloc[:, 2:].apply(lambda x: np.sqrt(x))
    elif normalize == 'rarefy':
        # rarefy
        pivot_table.iloc[:, 2:] = rarefy_table(pivot_table.iloc[:, 2:], depth=rarefy_depth)
    else:
        # no normalization
        pass
    return pivot_table


def prevalence_cutoff(df: pd.DataFrame, percent: float = 10, skip_columns: int = 2) -> pd.DataFrame:
    """
    Apply a prevalence cutoff to the DataFrame, removing features that do not
    appear in at least a certain percentage of samples.
    This is useful for filtering out low-prevalence features that may not be
    biologically relevant.

    Args:
        df (pd.DataFrame): The input DataFrame containing feature abundances.
        percent (float): The prevalence threshold as a percentage.
        skip_columns (int): The number of columns to skip (e.g., taxonomic information).

    Returns:
        pd.DataFrame: A filtered DataFrame with low-prevalence features removed.
    """
    # Calculate the number of samples
    num_samples = df.shape[1] - skip_columns
    # Calculate the prevalence threshold
    threshold = (percent / 100) * num_samples
    # Filter features based on prevalence
    filtered = df.loc[df.iloc[:, skip_columns:].gt(0).sum(axis=1) >= threshold]
    # Reset the index
    filtered = filtered.reset_index(drop=True)
    return filtered


def rarefy_table(df: pd.DataFrame, depth: int = None, axis: int = 1) -> pd.DataFrame:
    """
    Rarefy an abundance table to a given depth. If depth is None, uses the
    minimum sample sum across all samples.
    This function is a wrapper around the skbio.stats.subsample_counts function.
    
    Args:
        df: pd.DataFrame (rows: features, columns: samples)
        depth: int or None, rarefaction depth. If None, uses min sample sum.
        axis: int, 1 for samples in columns, 0 for samples in rows.

    Returns:
        pd.DataFrame: A rarefied DataFrame.
    """
    if axis == 1:
        sample_sums = df.sum(axis=0)
    else:
        sample_sums = df.sum(axis=1)
    if depth is None:
        depth = sample_sums.min()
    rarefied = {}
    print("Minimum rarefaction depth:", depth)
    for sample in df.columns if axis == 1 else df.index:
        counts = df[sample].values if axis == 1 else df.loc[sample].values
        if sample_sums[sample] < depth:
            # Not enough counts, fill with NaN or zeros
            rarefied[sample] = np.full_like(counts, np.nan)
        else:
            rarefied_counts = subsample_counts(counts.astype(int), int(depth))
            rarefied[sample] = rarefied_counts
    if axis == 1:
        return pd.DataFrame(rarefied, index=df.index)
    else:
        return pd.DataFrame(rarefied, columns=df.columns)


def split_metadata(metadata: pd.DataFrame, factor: str) -> Dict[str, list]:
    """
    Splits the metadata ref codes to dictionary of key being the factor value and
    value is a list of the ref codes.

    Args:
        metadata (pd.DataFrame): The input DataFrame containing metadata.
        factor (str): The column name to split the metadata by.
    Returns:
        Dict[str, list]: A dictionary with keys as unique values of the factor and
                         values as lists of ref codes.
    """
    if factor not in metadata.columns:
        raise ValueError(f"Factor '{factor}' not found in DataFrame columns.")
    # check if column is categorical
    if not isinstance(metadata[factor].dtype, pd.CategoricalDtype):
        raise ValueError(f"Column '{factor}' is not categorical (object dtype).")
    
    # for each unique value in the factor column, create a new table and append to a dictionary
    grouped_data = {}
    for value in metadata[factor].unique():
        # filter the dataframe
        filtered_df = metadata[metadata[factor] == value]
        # get the ref codes
        ref_codes = filtered_df['ref_code'].tolist()
        # add to the dictionary
        grouped_data[value] = ref_codes
    
    return grouped_data


def split_taxonomic_data_pivoted(taxonomy: pd.DataFrame, groups: Dict[str, list]) -> Dict[str, pd.DataFrame]:
    """
    Splits the taxonomic data into dictionary of DataFrames for each group.
    The split is based on the column names which need to match between the taxonomy DataFrame
    and the groups lists. The DataFrame should have a 'ncbi_tax_id' and 'taxonomic_concat' which 
    will serve as index of the resulting DataFrames.
    
    Args:
        taxonomy (pd.DataFrame): The input DataFrame containing taxonomic information.
        groups (Dict[str, list]): A dictionary where keys are unique values of a factor
            which correspond to the groups to split by.

    Returns:
        Dict[str, pd.DataFrame]: A dictionary with keys as unique values of the factor and
                                 values as DataFrames with separate columns for each taxonomic rank.
    """
    if not isinstance(groups, dict):
        raise ValueError("Groups must be a dictionary.")

    # for each unique value in the factor column, create a new table and append to a dictionary
    grouped_data = {}
    for value, ref_codes in groups.items():
        # filter the dataframe
        filtered_df = taxonomy[['ncbi_tax_id', 'taxonomic_concat'] + ref_codes]
        # remove rows with all zeros and print how many rows were removed
        len_before = len(filtered_df)
        filtered_df = filtered_df[filtered_df[ref_codes].sum(axis=1) != 0]
        len_after = len(filtered_df)
        print(f"Removed {len_before - len_after} rows with all zeros for {value}.")
        # check if the dataframe is empty
        if filtered_df.empty:
            print(f"Warning: No data for {value} in the taxonomic data.")
            continue
        # cprint how many rows were removed

        # add to the dictionary
        grouped_data[value] = filtered_df
    
    return grouped_data


def compute_bray_curtis(df: pd.DataFrame, skip_cols: int = 2, direction: str = 'samples') -> pd.DataFrame:
    """
    Compute Bray-Curtis dissimilarity and return as a pandas DataFrame.
    This function computes the Bray-Curtis dissimilarity for samples in the DataFrame.

    Args:
        df (pd.DataFrame): The input DataFrame containing sample counts.
        skip_cols (int): Number of columns to skip (e.g., taxonomic information).
        direction (str): Direction of the dissimilarity calculation, 'samples' or 'taxa'.
    
    Returns:
        pd.DataFrame: A DataFrame containing the Bray-Curtis dissimilarity matrix.
    """
    if direction not in ['samples', 'taxa']:
        raise ValueError("Direction must be either 'samples' or 'taxa'.")

    if direction == 'samples':
        # Use the sample IDs as the index
        ids = df.columns[skip_cols:].astype(str).tolist()
        result = beta_diversity(
            metric='braycurtis',
            counts=df.iloc[:, skip_cols:].T,
            ids=ids
        )
    elif direction == 'taxa':
        ids = df['ncbi_tax_id'].astype(str).tolist()
        result = beta_diversity(
            metric='braycurtis',
            counts=df.iloc[:, skip_cols:],
            ids=ids
        )

    
    bray_curtis_df = pd.DataFrame(
        result.data,
        index=ids,
        columns=ids
    )
    return bray_curtis_df


def fdr_pvals(p_spearman_df: pd.DataFrame) -> pd.DataFrame:
    """
    Apply FDR correction to the p-values DataFrame using Benjamini/Hochberg (non-negative)
    method. This function extracts the upper triangle of the p-values DataFrame.
    
    Args:
        p_spearman_df (pd.DataFrame): DataFrame containing p-values.
    
    Returns:
        pd.DataFrame: DataFrame with FDR corrected p-values.
    """
    # Extract upper triangle p-values
    pval_array = p_spearman_df.where(np.triu(np.ones(p_spearman_df.shape), k=1).astype(bool)).stack().values

    # Apply FDR correction
    _rejected, pvals_corrected, _, _ = multipletests(pval_array, alpha=0.05, method='fdr_bh')

    # Map corrected p-values back to a DataFrame
    pvals_fdr = p_spearman_df.copy()
    pvals_fdr.values[np.triu_indices_from(p_spearman_df, k=1)] = pvals_corrected
    pvals_fdr.values[np.tril_indices_from(p_spearman_df, k=0)] = np.nan  # Optional: keep only upper triangle
    return pvals_fdr


def interaction_to_graph(df, pos_cutoff=0.8, neg_cutoff=-0.6):
    """
    Create a network from the correlation matrix.
    Args:
        df (pd.DataFrame): The input DataFrame containing correlation values.
        pos_cutoff (float): Positive correlation cutoff.
        neg_cutoff (float): Negative correlation cutoff.
    Returns:
        nodes (list): List of node indices.
        edges_pos (list): List of positive edges.
        edges_neg (list): List of negative edges.
    """
    nodes, edges_pos, edges_neg = [], [], []
    count_pos, count_neg = 0, 0
    cols = df.columns.tolist()
    for i in range(df.shape[0]):
        nodes.append(cols[i])
        for j in range(i+1, df.shape[1]):
            if df.iloc[i, j] > pos_cutoff:
                edges_pos.append((cols[i], cols[j]))
                count_pos += 1
            # print(f"Sample {i} and Sample {j} have a high correlation of {df.iloc[i, j]}")
            elif df.iloc[i, j] < neg_cutoff:
                edges_neg.append((cols[i], cols[j]))
                count_neg += 1
                # print(f"Sample {i} and Sample {j} have a high negative correlation of {df.iloc[i, j]}")
    print(f"Number of positive edges: {count_pos}")
    print(f"Number of negative edges: {count_neg}")
    return nodes, edges_pos, edges_neg


def interaction_to_graph_with_pvals(df, pvals_df, pos_cutoff=0.8, neg_cutoff=-0.6, p_val_cutoff=0.05):
    """
    Create a network from the correlation matrix and p-values.
    Args:
        df (pd.DataFrame): The input DataFrame containing correlation values.
        pvals_df (pd.DataFrame): The DataFrame containing p-values.
        pos_cutoff (float): Positive correlation cutoff.
        neg_cutoff (float): Negative correlation cutoff.
    Returns:
        nodes (list): List of node indices.
        edges_pos (list): List of positive edges with p-values.
        edges_neg (list): List of negative edges with p-values.
    """
    nodes, edges_pos, edges_neg = [], [], []
    count_pos, count_neg = 0, 0
    cols = df.columns.tolist()
    for i in range(df.shape[0]):
        nodes.append(cols[i])
        for j in range(i+1, df.shape[1]):
            if df.iloc[i, j] > pos_cutoff and pvals_df.iloc[i, j] < p_val_cutoff:
                edges_pos.append((cols[i], cols[j]))
                count_pos += 1
            elif df.iloc[i, j] < neg_cutoff and pvals_df.iloc[i, j] < p_val_cutoff:
                edges_neg.append((cols[i], cols[j]))
                count_neg += 1
    print(f"Number of positive edges: {count_pos}")
    print(f"Number of negative edges: {count_neg}")
    return nodes, edges_pos, edges_neg
```

```{code-cell}
# Load metadata
full_metadata = get_full_metadata()

# filter the metadata only for valid 181 samples
valid_samples = get_valid_samples()
full_metadata = enhance_metadata(full_metadata, valid_samples)

mgf_parquet_dfs = get_data()
ssu = mgf_parquet_dfs['ssu']
ssu.head()
```

```{code-cell}
assert full_metadata.shape[0] == 181

print(full_metadata.shape)
```

Batch 1+2 contain 181 valid samples, which should be the first value printed above. Now we convert `object` columns to categorical variables

```{code-cell}
for col in full_metadata.columns:
    # print(full_metadata[col].isinstance(pd.CategoricalDtype()))
    # check if object dtype
    if full_metadata[col].dtype == 'object':
        # convert to categorical
        full_metadata[col] = full_metadata[col].astype('category')

# example how to crosscheck the abve operation
isinstance(full_metadata['season'].dtype, pd.CategoricalDtype)
```

## 2. Remove high taxa

```{code-cell}
print("Original DataFrame shape:", ssu.shape)
ssu = remove_high_taxa(ssu, tax_level='phylum')
print("Filtered DataFrame shape:", ssu.shape)
```

## 3. Pivot tables

```{code-cell}
ssu_pivot = pivot_taxonomic_data(ssu, normalize=None, rarefy_depth=None)
ssu_pivot.head()
```

## 4. Remove low abundance taxa

```{code-cell}
ssu_filt = prevalence_cutoff(ssu_pivot, percent=10, skip_columns=2)
ssu_filt.head()
```

## 5. Rarefy or normalize

```{code-cell}
ssu_rarefied = ssu_filt.copy()
ssu_rarefied.iloc[:, 2:] = rarefy_table(ssu_filt.iloc[:, 2:], depth=None, axis=1)
ssu_rarefied.head()
```

## 6. Remove duplicates

Anything which has count > 1 is replicate of some sort.

```{code-cell}
full_metadata['replicate_info'].value_counts()
```

Remove them here

```{code-cell}
filtered_metadata = full_metadata.drop_duplicates(subset='replicate_info', keep='first')
filtered_metadata.shape
```

## 7. Split to groups on chosen `factor`

Split metadata to groups

```{code-cell}
FACTOR = 'season'
groups = split_metadata(
    filtered_metadata,
    FACTOR
)

# remove groups which have less than 2 members (bad for your statistics :)
for groups_key in list(groups.keys()):
    print(f"{groups_key}: {len(groups[groups_key])} samples")
    if len(groups[groups_key]) < 3:
        del groups[groups_key]
        print(f"Warning: {groups_key} has less than 3 samples, therefore removed.")
```

Split pivoted taxonomy data to groups:

```{code-cell}
# pivot and filter the ssu
split_taxonomy = split_taxonomic_data_pivoted(
    ssu_rarefied,
    groups
)

for key, value in split_taxonomy.items():
    print(f"{key}: {value.shape[0]} rows, {value.shape[1]} columns")
# split_taxonomy.keys(),

print(ssu_rarefied.shape)
split_taxonomy['Autumn'].iloc[:, 2:].head()

```

## 8. Calculate associations (Bray-curits dissimilarity and Spearman's)

Bray-curtis dissimilarity:

```{code-cell}
bray_taxa = {}
for factor, df in split_taxonomy.items():
    bray_taxa[factor] = compute_bray_curtis(df, direction='taxa')
```

Spearman correlation:

```{code-cell}
spearman_taxa = {}
# Compute Spearman correlation
for factor, df in split_taxonomy.items():
    corr, p_spearman = spearmanr(df.iloc[:, 2:].T)
    assert corr.shape == p_spearman.shape, "Spearman correlation and p-values must have the same shape."
    corr_df = pd.DataFrame(
        corr,
        index=df['ncbi_tax_id'],
        columns=df['ncbi_tax_id']
    )
    p_spearman_df = pd.DataFrame(
        p_spearman,
        index=df['ncbi_tax_id'],
        columns=df['ncbi_tax_id']
    )
    d = {
        'correlation': corr_df,
        'p_vals': p_spearman_df
    }
    spearman_taxa[factor] = d
    assert spearman_taxa[factor]['correlation'].shape == spearman_taxa[factor]['p_vals'].shape, "Spearman correlation and p-values must have the same shape."


spearman_taxa['Summer']['correlation'].head()
```

## 9. False discovery rate (FDR) correction

```{code-cell}
for factor, d in spearman_taxa.items():
    pvals_fdr = fdr_pvals(d['p_vals'])
    spearman_taxa[factor]['p_vals_fdr'] = pvals_fdr
```

Let's crosscheck the shapes of one particular set of DFs for one factor value

```{code-cell}
season = 'Summer'
spearman_taxa[season]['correlation'].shape, spearman_taxa[season]['p_vals'].shape, spearman_taxa[season]['p_vals_fdr'].shape
```

We will plot the FDR curve a bit later, patiance please.

## 10. Build network per group

Bray-curtis

```{code-cell}
# histogram of the correlation values for setting graph cutoffs
plt.figure(figsize=(10, 5))
for factor, df in bray_taxa.items():
    plt.hist(df.values.flatten(), bins=50, alpha=0.5, label=factor, log=True)

plt.title("Histogram of Correlation Values")
plt.xlabel("Correlation")
plt.ylabel("Frequency")
plt.legend()
plt.show()
```

Spearman correlation

```{code-cell}
# histogram of the correlation values for setting graph cutoffs
plt.figure(figsize=(10, 5))
for factor, dict_df in spearman_taxa.items():
    plt.hist(dict_df['correlation'].values.flatten(), bins=50, alpha=0.5, label=factor, log=True)
plt.title("Histogram of Correlation Values")
plt.xlabel("Correlation")
plt.ylabel("Frequency")
plt.legend()
plt.show()

plt.figure(figsize=(10, 5))
for factor, dict_df in spearman_taxa.items():
    plt.scatter(
        dict_df['p_vals'].values.flatten()[::10],  # downsample for better visibility and speed
        dict_df['p_vals_fdr'].values.flatten()[::10],
        alpha=0.5,
        label=factor,
    )
plt.axvline(0.05, color='gray', linestyle='--', label='Raw p=0.05')
plt.axhline(0.05, color='black', linestyle='--', label='FDR p=0.05')
plt.xlabel('Raw p-value')
plt.ylabel('FDR-corrected p-value')
plt.title(f'Comparison of Raw and FDR-corrected p-values for {factor}')
plt.legend()
plt.show()
```

## 11. Identify positive and negative associations (Example)

```{code-cell}
for factor, df in spearman_taxa.items():
    print(f"factor: {factor}")
    nodes, edges_pos, edges_neg = interaction_to_graph_with_pvals(
        dict_df['correlation'],
        dict_df['p_vals_fdr'],
        pos_cutoff=0.7,
        neg_cutoff=-0.5,
        p_val_cutoff=0.05
    )
```

## 12. Downstream analysis

We are going to look at:

- `degree centrality`, which measures the number of direct connections a node (taxon) has. High degree centrality suggests a taxon is highly connected and may play a central role in the community. Therefore represents a so called `keystone species`.
- `betweenness`, which represents how often a node appears on the shortest paths between other nodes. Taxa with high/low betweenness may act as bridges or bottlenecks in the network, respectively.

```{code-cell}
network_results = {}
for factor, dict_df in spearman_taxa.items():
    print(f"Factor: {factor}")
    nodes, edges_pos, edges_neg = interaction_to_graph_with_pvals(dict_df['correlation'], dict_df['p_vals_fdr'], pos_cutoff=0.7, neg_cutoff=-0.5, p_val_cutoff=0.05)

    G = nx.Graph(
        mode = factor,
    )

    G.add_nodes_from(nodes)
    G.add_edges_from(edges_pos, color='green')
    G.add_edges_from(edges_neg, color='red')

    network_results[factor] = {
        "graph": G,
        "nodes": nodes,
        "edges_pos": edges_pos,
        "edges_neg": edges_neg
    }

    colors = nx.get_edge_attributes(G, 'color')

    degree_centrality = nx.degree_centrality(G)

    network_results[factor]['degree_centrality'] = sorted(degree_centrality.items(),
                                                          key=lambda x: x[1],
                                                          reverse=True)[:10]
    
    betweenness = nx.betweenness_centrality(G)

    network_results[factor]['top_betweenness'] = sorted(betweenness.items(),
                                                    key=lambda x: x[1],
                                                    reverse=True)[:10]
    network_results[factor]['bottom_betweenness'] = sorted(betweenness.items(),
                                                           key=lambda x: x[1])[:10]
    network_results[factor]['total_nodes'] = G.number_of_nodes()
    network_results[factor]['total_edges'] = G.number_of_edges()
```

Compile dataframes from the dictionary, not ideal yet.

```{code-cell}
DF = pd.DataFrame(columns=[FACTOR, 'centrality', 'top_betweenness', 'bottom_betweenness', 'total_nodes', 'total_edges'])
for factor, dict_results in network_results.items():
    DF = pd.concat([DF, pd.DataFrame([{
        FACTOR: factor,
        'centrality': dict_results['degree_centrality'],
        'top_betweenness': dict_results['top_betweenness'],
        'bottom_betweenness': dict_results['bottom_betweenness'],
        'total_nodes': dict_results['total_nodes'],
        'total_edges': dict_results['total_edges']
    }])], ignore_index=True)
DF.head()


df_centrality = pd.DataFrame(columns=[FACTOR, 'node', 'centrality'])

for factor, dict_results in network_results.items():
    for node, centrality in dict_results['degree_centrality']:
        df_centrality = pd.concat([df_centrality, pd.DataFrame({
            FACTOR: [factor],
            'node': [node],
            'centrality': [centrality]
        })], ignore_index=True)

# Set hierarchical (MultiIndex) with 'factor' above 'node'
# df_centrality.set_index([FACTOR, 'node'], inplace=True)
df_centrality
```

The node number identifiers correspond to the NCBI tax_id. Is any node shared in the high centrality table?

```{code-cell}
sum(df_centrality['node'].value_counts() > 1)
```

Whole table

```{code-cell}
DF.head()
```

### Jaccard similarity of edge sets

This allows to compare networks pair-wise. Since we implemented coloring of the positive and negative associated edges, we can perfor this 3 times (for total, positive, and negative). Three seasons pair-wise give three nuber as well.

First we calculate all the necessary sets of edges:

```{code-cell}
# this is a pair-wise comparison of the networks
edges_summer_pos = set(network_results['Summer']['edges_pos'])
edges_summer_neg = set(network_results['Summer']['edges_neg'])
edges_summer_total = set(network_results['Summer']['graph'].edges())


edges_autumn_pos = set(network_results['Autumn']['edges_pos'])
edges_autumn_neg = set(network_results['Autumn']['edges_neg'])
edges_autumn_total = set(network_results['Autumn']['graph'].edges())

edges_winter_pos = set(network_results['Winter']['edges_pos'])
edges_winter_neg = set(network_results['Winter']['edges_neg'])
edges_winter_total = set(network_results['Winter']['graph'].edges())
```

And now Jaccard distances

```{code-cell}
print('Analyse edges identified with POSITIVE correlation:')
print('###########################################')
intersection = edges_summer_pos & edges_autumn_pos
union = edges_summer_pos | edges_autumn_pos
jaccard_similarity = len(intersection) / len(union)
print(f"Jaccard similarity between Summer and Autumn: {jaccard_similarity:.4f}")

intersection = edges_summer_pos & edges_winter_pos
union = edges_summer_pos | edges_winter_pos
jaccard_similarity = len(intersection) / len(union)
print(f"Jaccard similarity between Summer and Winter: {jaccard_similarity:.4f}")

intersection = edges_autumn_pos & edges_winter_pos
union = edges_autumn_pos | edges_winter_pos
jaccard_similarity = len(intersection) / len(union)
print(f"Jaccard similarity between Autumn and Winter: {jaccard_similarity:.4f}\n")


print('Analyse all edges:')
print('###########################################')
intersection = edges_summer_total & edges_autumn_total
union = edges_summer_total | edges_autumn_total
jaccard_similarity = len(intersection) / len(union)
print(f"Jaccard similarity between Summer and Autumn: {jaccard_similarity:.4f}")

intersection = edges_summer_total & edges_winter_total
union = edges_summer_total | edges_winter_total
jaccard_similarity = len(intersection) / len(union)
print(f"Jaccard similarity between Summer and Winter: {jaccard_similarity:.4f}")

intersection = edges_autumn_total & edges_winter_total
union = edges_autumn_total | edges_winter_total
jaccard_similarity = len(intersection) / len(union)
print(f"Jaccard similarity between Autumn and Winter: {jaccard_similarity:.4f}\n")


print('Analyse edges identified with NEGATIVE correlation:')
print('###########################################')
intersection = edges_summer_neg & edges_autumn_neg
union = edges_summer_neg | edges_autumn_neg
jaccard_similarity = len(intersection) / len(union)
print(f"Jaccard similarity between Summer and Autumn: {jaccard_similarity:.4f}")

intersection = edges_summer_neg & edges_winter_neg
union = edges_summer_neg | edges_winter_neg
jaccard_similarity = len(intersection) / len(union)
print(f"Jaccard similarity between Summer and Winter: {jaccard_similarity:.4f}")

intersection = edges_autumn_neg & edges_winter_neg
union = edges_autumn_neg | edges_winter_neg
jaccard_similarity = len(intersection) / len(union)
print(f"Jaccard similarity between Autumn and Winter: {jaccard_similarity:.4f}\n")
```

### Graph example visualization

```{code-cell}
fig, axes = plt.subplots(1, 3, figsize=(18, 6))

for ax, season in zip(axes, ['Summer', 'Autumn', 'Winter']):
    G = network_results[season]['graph']
    colors = nx.get_edge_attributes(G, 'color')
    pos = nx.spring_layout(G, k=0.2, iterations=50, seed=42)
    nx.draw_networkx_nodes(G, pos, ax=ax, alpha=0.2, node_color='grey', node_size=15)
    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.2, edge_color=list(colors.values()))
    ax.set_title(season)
    ax.axis('off')

plt.tight_layout()
plt.show()
```

## Summary

None of this answers any of the possible question you want to ask, but with the growing literature applying complex network methods to metagenomics data in combination of unique unified sampling and sequencing approach by EMO-BON, this dataset is ideal to push the boundaries of monitoring using the above described tool set.
