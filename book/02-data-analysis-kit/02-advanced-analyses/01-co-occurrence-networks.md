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
David Palecek, August 8, 2025
:::

Complex networks provide metrics and to identify kestone species bridges/bottleneck species within the community. EMO-BON data are unified in methods but highly diverse in respect to the metadata (sampling station, temperature, season etc.), which can be interpreted as control/treatment. To refrase a particular question in respect to seasonality, is there structural difference in taxonomy buildup in the communities between spring, summer, autumn and winter.

We will construct association matrices between taxa over the series of sampling events grouped by the factor (season in this case).

## Recipe

1. Load data
2. Remove high taxa (non-identified sequences)
3. Pivot table from sampling events to taxa.
4. Remove low abundance taxa
5. Rarefy, or normalize
6. Remove replicates
7. Split to groups on chosen `factor`
8. Calculate associations (Bray-curits dissimilarity, Spearman's correlation, etc.)
9. False discovery rate correction
10. Build and analyse network per group

## 1. Loading data and metadata

```{code-cell}
:class: dropdown
import os
from typing import Dict
import numpy as np
import pandas as pd

from scipy.stats import spearmanr
import networkx as nx

from momics.taxonomy import (
    remove_high_taxa,
    pivot_taxonomic_data,
    normalize_abundance,
    prevalence_cutoff,
    rarefy_table,
    split_metadata,
    split_taxonomic_data,
    split_taxonomic_data_pivoted,
    compute_bray_curtis,
    fill_taxonomy_placeholders,
    fdr_pvals,
)
from momics.utils import load_and_clean
from momics.networks import (
    interaction_to_graph,
    interaction_to_graph_with_pvals,
    pairwise_jaccard_lower_triangle,
    build_interaction_graphs,
)
from momics.plotting import plot_network
from momics.stats import (
    spearman_from_taxonomy,
    plot_fdr,
    plot_association_histogram,
)

from mgo.udal import UDAL

# All low level functions are imported from the momics package
from momics.loader import load_parquets_udal
from momics.metadata import get_metadata_udal, enhance_metadata

import seaborn as sns
import matplotlib.pyplot as plt
```

All the methods which are not part of `marine-omics` package are defined below. Documnetation and context for the `marine-omics (momics)` methods can be found on [readthedocs.io](https://marine-omics-methods.readthedocs.io/en/latest/index.html).

```{code-cell}
:label: methods
:class: dropdown

def get_valid_samples():
    url = "https://raw.githubusercontent.com/emo-bon/momics-demos/refs/heads/main/data/shipment_b1b2_181.csv"
    df_valid = pd.read_csv(
        url, header=0, index_col=0,
    )
    return df_valid
```

```{code-cell}
valid_samples = get_valid_samples()

# load, filter and enhance metadata
full_metadata, mgf_parquet_dfs = load_and_clean(valid_samples=valid_samples)

ssu = mgf_parquet_dfs['ssu']
ssu.head()
```

```{code-cell}
assert full_metadata.shape[0] == 181

print(full_metadata.shape)
```

Batch 1+2 contain 181 valid samples, which should be the first value printed above. Now we convert `object` columns to categorical variables

```{code-cell}
# Identify object columns
categorical_columns = sorted(full_metadata.select_dtypes(include=['object', "boolean"]).columns)

# Convert them all at once to category
full_metadata = full_metadata.astype({col: 'category' for col in categorical_columns})

if not isinstance(full_metadata['season'].dtype, pd.CategoricalDtype):
        raise ValueError(f"Column 'season' is not categorical (object dtype).")
```

## 2. Remove high taxa

First infer missing higher taxa if some lower taxon was identified using `unclassified_` placeholder. Example of filled row of taxonomy can look like this:

kingdom phylum	class	order	family	genus	species
Archaea unclassified_Thaumarchaeota	Thaumarchaeota	unclassified_Cenarchaeales	Cenarchaeales	Cenarchaeaceae	None	None

```{code-cell}
taxonomy_cols = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
ssu_filled = fill_taxonomy_placeholders(ssu, taxonomy_cols)
ssu_filled.head(10)
```

Because many sequences are not well classified, information that the sequence is from superKingdom Bacteria is not helpful and does not add information value to our analysis. Here we remove all taxa which are only identified on `phylum` level and higher.

```{code-cell}
print("Original DataFrame shape:", ssu.shape)
ssu_filled = remove_high_taxa(ssu_filled, taxonomy_ranks=taxonomy_cols, tax_level='phylum')
print("Filtered DataFrame shape:", ssu.shape)
```

## 3. Pivot tables

Pivoting converts taxonomy table which rows contain each taxon identified in every sample to table with rows of unique taxa with columns representing samples and values of abundances per sample.

```{code-cell}
ssu_pivot = pivot_taxonomic_data(ssu_filled)
ssu_pivot.head()
```

:::{note} Tidy data 💡
:class: dropdown
This is the step, when we violate tidy data principles, because the abundance is no longer single column (feature) and it is not even clear what the values in the table are now.
:::

## 4. Remove low prevalence taxa

The cutoff selected here is 10 %, which means that all taxa which do not appear at least in 10 % of all the samples in the taxonomy table will be removed. This doees not refer to low abundance taxa within single samples.

```{code-cell}
ssu_filt = prevalence_cutoff(ssu_pivot, percent=10, skip_columns=0)
ssu_filt.head()
```

## 5. Rarefy or normalize

Many normalization techniques are readily available to mornalize taxonomy data. For the purposes of co-occurrence networks, it is reccomended to rarefy [refs]. Please correct if wrong.

```{code-cell}
ssu_rarefied = rarefy_table(ssu_filt, depth=None, axis=1)
ssu_rarefied.head()
```

## 6. Remove replicates

Anything which has `replicate info` count > 1 is replicate of some sort. Note that this information comes from the enhanced metadata information, because `replicate info` is not in the original metadata. It is a combination of `observatory`, `sample type`, `date`, and `filter` information.

It should be considered to filter replicates before steps 4. and 5. Identify replicates:

```{code-cell}
full_metadata['replicate info'].value_counts()
```

Remove the replicates:

```{code-cell}
filtered_metadata = full_metadata.drop_duplicates(subset='replicate info', keep='first')
filtered_metadata.shape
```

## 7. Split to groups on chosen `factor`

Split metadata to groups based on unique values in your factor.

::::{caution}
Even for the case of selected factor of `season`, which has only 4 valid values, the rule of thumb of having 25 samples per group is broken.

Below we only check and remove groups which have only 2 samples, which is too little to calculate any meaningful statistics.
:class: dropdown
::::

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

Split pivoted taxonomy data according to groups of `ref_code`s defined above:

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

## 8. Calculate associations (Bray-curtis and Spearman's)

Bray-curtis dissimilarity:

```{code-cell}
bray_taxa = {}
for factor, df in split_taxonomy.items():
    bray_taxa[factor] = compute_bray_curtis(df, skip_cols=0, direction='taxa')
```

Spearman correlation:

```{code-cell}
spearman_taxa = spearman_from_taxonomy(split_taxonomy)

spearman_taxa['Summer']['correlation'].head()
```

## 9. False discovery rate (FDR) correction

### Calculate FDR

We have more than 1000 taxa in the tables, which means 1000*1000 / 2 association values. Too many to statistically trust that any low p-value is not just statistical coincidence. Therefore False Discovery Rate correction is absolutely essential here.

```{code-cell}
PVAL = 0.05  # Define a p-value threshold for FDR correction
for factor, d in spearman_taxa.items():
    pvals_fdr = fdr_pvals(d['p_vals'], pval_cutoff=PVAL)
    spearman_taxa[factor]['p_vals_fdr'] = pvals_fdr
```

Let's crosscheck one particular set of DFs for one factor value, which should all have the same shape:

```{code-cell}
factor = 'Summer'
spearman_taxa[factor]['correlation'].shape, spearman_taxa[factor]['p_vals'].shape, spearman_taxa[factor]['p_vals_fdr'].shape
```

### Display histograms FDR corrected p-values

For Bray-curtis dissimilarity we plot only the histogram. Does it make sense to do FDR? How?

```{code-cell}
plot_association_histogram(bray_taxa)
```

We performed FDR on Spearman's correlations, so both correlation values and FDR curves are displayed.

```{code-cell}
plot_association_histogram(spearman_taxa)
plot_fdr(spearman_taxa, PVAL)
```

The FDR curve shows huge inflation of significant associations if taken from the raw data. The above print show the stats on how associations got discarded.

## 10. Build network for each factor value group

### Nodes and edges for positive and negative associations (Example)

Code example for generating a list of nodes and significant positive/negative edges to construct the graph from.

```{code-cell}
for factor, dict_df in spearman_taxa.items():
    print(f"factor: {factor}")
    nodes, edges_pos, edges_neg = interaction_to_graph_with_pvals(
        dict_df['correlation'],
        dict_df['p_vals_fdr'],
        pos_cutoff=0.7,
        neg_cutoff=-0.5,
        p_val_cutoff=0.05
    )
    break
```

### Network generation and metrics calculation

At the same time as we generate the graph, we calculate several descriptive metrics of themetwork (implement saving of the network, if you need the graphs for later). We are going to look at:

- `degree centrality`, which measures the number of direct connections a node (taxon) has. High degree centrality suggests a taxon is highly connected and may play a central role in the community. Therefore represents a so called `keystone species`.
- `betweenness`, which represents how often a node appears on the shortest paths between other nodes. Taxa with high/low betweenness may act as bridges or bottlenecks in the network, respectively.

```{code-cell}
network_results = build_interaction_graphs(spearman_taxa)
```

### Jaccard similarity of edge sets

Jaccard similarity evaluates pairwise similarity of constructed networks.

For positive edges:

```{code-cell}
df_jaccard = pairwise_jaccard_lower_triangle(network_results, edge_type='edges_pos')
df_jaccard
```

For negative edges:

```{code-cell}
df_jaccard = pairwise_jaccard_lower_triangle(network_results, edge_type='edges_neg')
df_jaccard
```

For all edges altogether:

```{code-cell}
df_jaccard = pairwise_jaccard_lower_triangle(network_results, edge_type='all')
df_jaccard
```

Compile dataframes from the dictionary. The dataframe structure is not ideal yet, feel free to open a PR [here](https://github.com/fair-ease/book-marine-omics-observation/pulls).

Here are two functions to generate quite silly dataframes of results:

```{code-cell}
def network_results_df(network_results, factor_name):
    out = pd.DataFrame(columns=[factor_name, 'centrality', 'top_betweenness', 'bottom_betweenness', 'total_nodes', 'total_edges'])
    for factor, dict_results in network_results.items():
        out = pd.concat([out, pd.DataFrame([{
            factor_name: factor,
            'centrality': dict_results['degree_centrality'],
            'top_betweenness': dict_results['top_betweenness'],
            'bottom_betweenness': dict_results['bottom_betweenness'],
            'total_nodes': dict_results['total_nodes'],
            'total_edges': dict_results['total_edges']
        }])], ignore_index=True)

    return out


def create_centrality_dataframe(network_results, factor_name):
    df_centrality = pd.DataFrame(columns=[factor_name, 'node', 'centrality'])

    for factor, dict_results in network_results.items():
        for node, centrality in dict_results['degree_centrality']:
            df_centrality = pd.concat([df_centrality, pd.DataFrame({
            factor_name: [factor],
            'node': [node],
            'centrality': [centrality]
        })], ignore_index=True)

    return df_centrality
```

Generate the dataframes

```{code-cell}
df_results = network_results_df(network_results, FACTOR)
df_centrality = create_centrality_dataframe(network_results, FACTOR)
```

The node identifiers correspond to the NCBI tax_id. Is any node shared in the high centrality table?

```{code-cell}
sum(df_centrality['node'].value_counts() > 1)
```

Whole table looks like this

```{code-cell}
df_results.head()
```

### Graph example visualization

```{code-cell}
plot_network(network_results, spearman_taxa)
```

If you do not do the FDR and remove high taxa, the networks will be hard to distinguis visually, mostly because of all the false positive association hits. Quite surprisingly, out of 1000*1000 tables, we get a clear visual separation of the community interaction patterns.

## Streamline the whole NB

It is now simple to define high-level method which takes dictionary of all the parametres as input and runs all of the above functions.

```{code-cell}
def network_analysis(params):
    # prepare metadata
    if params['drop_metadata_duplicates']:
        full_metadata = params['metadata'].drop_duplicates(subset='replicate info', keep='first')
    else:
        full_metadata = params['metadata']

    ssu_filled = fill_taxonomy_placeholders(params["data"], params["taxonomy_cols"])

    print("Original DataFrame shape:", ssu_filled.shape)
    ssu_filled = remove_high_taxa(ssu_filled,
                                  taxonomy_ranks=params["taxonomy_cols"],
                                  tax_level=params["remove_high_taxa_level"])
    print("Filtered DataFrame shape:", ssu_filled.shape)

    ssu_pivot = pivot_taxonomic_data(ssu_filled)

    del ssu_filled
    ssu_filt = prevalence_cutoff(ssu_pivot, percent=params['prevalence_cutoff'], skip_columns=0)
    del ssu_pivot
    ssu_rarefied = rarefy_table(ssu_filt, depth=params['rarefy_depth'], axis=1)
    del ssu_filt

    groups = split_metadata(full_metadata, params['factor'])

    for groups_key in list(groups.keys()):
        print(f"{groups_key}: {len(groups[groups_key])} samples")
        if len(groups[groups_key]) < 3:
            del groups[groups_key]
            print(f"Warning: {groups_key} has less than 3 samples, therefore removed.")

    split_taxonomy = split_taxonomic_data_pivoted(ssu_rarefied, groups)

    spearman_taxa = spearman_from_taxonomy(split_taxonomy)

    for factor, d in spearman_taxa.items():
        pvals_fdr = fdr_pvals(d['p_vals'], pval_cutoff=params['fdr_cutoff'])
        spearman_taxa[factor]['p_vals_fdr'] = pvals_fdr

    plot_fdr(spearman_taxa, params['fdr_cutoff'])
    plot_association_histogram(spearman_taxa)

    network_results = build_interaction_graphs(
        spearman_taxa,
        pos_cutoff=params["network_pos_cutoff"],
        neg_cutoff=params["network_neg_cutoff"],
        p_val_cutoff=params["fdr_cutoff"],
    )

    print('Jaccard Similarity (Positive Edges):')
    display(pairwise_jaccard_lower_triangle(network_results, edge_type='edges_pos'))
    print('Jaccard Similarity (Negative Edges):')
    display(pairwise_jaccard_lower_triangle(network_results, edge_type='edges_neg'))
    print('Jaccard Similarity (All Edges):')
    display(pairwise_jaccard_lower_triangle(network_results, edge_type='all'))

    res_df = network_results_df(network_results, params['factor'])
    display(res_df)

    df_centrality = create_centrality_dataframe(network_results, params['factor'])
    display(df_centrality)

    plot_network(network_results, spearman_taxa)
```

Define your inputs and run the whole pipeline

```{code-cell}
params = {
    "data": ssu,
    "metadata": full_metadata,
    "taxonomy_cols": ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'],
    "remove_high_taxa_level": 'phylum',
    'prevalence_cutoff': 10,
    'rarefy_depth': None,
    'drop_metadata_duplicates': True,
    'factor': 'environment (material)',
    'fdr_cutoff': 0.05,
    "network_pos_cutoff": 0.5,
    "network_neg_cutoff": -0.3,
}

# run everything
# network_analysis(params)
```

## Summary

None of this answers any of the possible question you want to ask, but with the growing literature applying complex network methods to metagenomics data in combination of unique unified sampling and sequencing approach by EMO-BON, this dataset is ideal to push the boundaries of monitoring using the above described tool set.

Factors can be any categorical variable from the metadata table. However alternative approach to color and evaluate the associations is to color the graph but some knowledge about the taxa, such r- or K- community species classification etc.
