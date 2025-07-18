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

# EMO-BON vs other campaign at MGnify

:::{note} Last update 👈
:class: dropdown
David Palecek, July 17, 2025
:::

Many campaigns have been analysed with MGnify pipeline. EMO-BON data use metaGOflow pipeline, which is a `reads` subworkflow of the MGnify. Therefore the summury tables available at MGnify are directly comparable with the EMO-BON metaGOflow outputs.

Here we (a) select several campaigns related to the European region at MGnify, (b) query taxonomic sumaries for the whole campaign based on the `study ID` and (c) provide an overview visualization.

Overall, this is a starting point for taxonomic universal 

::::{note} Even though the function to query the data via the [API](https://www.ebi.ac.uk/metagenomics/api/v1/) is provided below, we use pre-downloaded files here.
:class: dropdown
::::

## Setup and data loading

```{code-cell}
:class: dropdown
import os

from functools import partial

import pandas as pd
import panel as pn
from dotenv import load_dotenv
load_dotenv()

# All low level functions are imported from the momics package
from momics.loader import load_parquets_udal

from momics.taxonomy import (
    fill_taxonomy_placeholders,
    pivot_taxonomic_data,
    prevalence_cutoff,
)
```

Below is the method used to query and save. Run it in the loop if you have a list of studies to download

```{code-cell}
:label: API call method
:class: dropdown

from jsonapi_client import Session as APISession
from urllib.request import urlretrieve

def retrieve_summary(studyId: str, matching_string: str = 'Taxonomic assignments SSU') -> None:
    """
    Retrieve summary data for a given analysis ID and save it to a file. Matching strings 
    are substrings of for instance:
    - Phylum level taxonomies SSU
    - Taxonomic assignments SSU
    - Phylum level taxonomies LSU
    - Taxonomic assignments LSU
    - Phylum level taxonomies ITSoneDB
    - Taxonomic assignments ITSoneDB
    - Phylum level taxonomies UNITE
    - Taxonomic assignments UNITE

    Example usage:
    retrieve_summary('MGYS00006680', matching_string='Taxonomic assignments SSU')

    Args:
        studyId (str): The ID of the analysis to retrieve. studyId is the MGnify study ID, used
            also to save the output .tsv file.
        matching_string (str): The string to match in the download description label.
    
    Returns:
        None
    """
    with APISession("https://www.ebi.ac.uk/metagenomics/api/v1") as session:
        for download in session.iterate(f"studies/{studyId}/downloads"):
            if download.description.label == matching_string:
                print(f"Downloading {download.alias}...")
                urlretrieve(download.links.self.url, f'{studyId}.tsv')
```

```python
# example of execution
retrieve_summary('MGYS00006680', matching_string='Taxonomic assignments SSU')
```

SSU table is used for this demonstration, but the methods shows are agnostic to taxonomic table origin. Metadata are not treated in this example, because their harmonization requires tedious manual work.

```{code-cell}
# load data 
mgf_parquet_dfs = load_parquets_udal()
ssu = mgf_parquet_dfs['ssu']
ssu.head()
```

Define a taxonomy rank strings as a global variable

```{code-cell}
TAXONOMY_RANKS = ['superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
```

## Pivot EMO-BON data

EMO-BON taxonomy tables

:::{caution} DP, please unify for the future versions.
:class: dropdown
This originates in taxonomy methods in marine-omics-methods package. And can be changed.
:::

```{code-cell}
def clean_tax_row(row):
    """
    Cleans the taxonomic row by removing empty strings and replacing spaces with underscores.
    """
    # replace string with underscores
    row = row.replace('_', '__')
    split_row = row.split(';')
    res = [split_row[1]]

    for tax in split_row[3:]:
        if tax[-1] == '_':
            break
        res.append(tax)
    
    return ';'.join(res)

# fill missing higher taxa (usually happens for tentative taxa)
ssu_filt = fill_taxonomy_placeholders(ssu, TAXONOMY_RANKS)

## pivot the abundance table
ssu_filt = pivot_taxonomic_data(ssu_filt, normalize=None, rarefy_depth=None)

# remove tax id
ssu_filt = ssu_filt.drop(columns=['ncbi_tax_id'])

ssu_filt['taxonomic_concat'] = ssu_filt['taxonomic_concat'].apply(clean_tax_row)
# rename columns, to agree with the MGnify tables
ssu_filt = ssu_filt.rename(columns={
    'taxonomic_concat': '#SampleID',
})
ssu_filt.head()
```

## Load MGnify datasets

The list of selected datasets is far from exhausted but includes:

| Study ID | Description |
| :------- | ---------: |
| MGYS00006608 | 16S rRNA amplicon sequencing from the Ocean Sampling Day (OSD), June 2018 |
| MGYS00006607 | 16S rRNA amplicon sequencing from the Ocean Sampling Day (OSD), June 2019 |
| MGYS00000492 | Amplicon sequencing of Tara Oceans DNA samples corresponding to size fractions for prokaryotes or protist |
| MGYS00006680 | SOLA sampling point Raw sequence reads |
| MGYS00006682 | Vertical stratification of environmental DNA in the open ocean captures ecological patterns and behavior of deep-sea fishes |
| MGYS00006678 | Dataset on spatiotemporal variation of microbial plankton communities in the Baltic Sea |
| MGYS00006675 | 16S rRNA gene amplicon time-series in Blanes Bay Microbial Observatory (BBMO) |
| MGYS00003725 | Arctic microbiome along Svalbard Cross Shelf transects |
| MGYS00006686 | Environmental DNA and zooplankton samples taken at Helgoland Roads in June 2019 |
| MGYS00006714 | Regional and vertical patterns in microbial communities across Fram Strait (2015-2019) |

```{code-cell}
ds_list = {
    'OSD-2018': 'mgnify_data/ERP124424_taxonomy_abundances_SSU_v5.0.tsv',
    'OSD-2019': 'mgnify_data/ERP124431_taxonomy_abundances_SSU_v5.0.tsv',
    'Tara': 'mgnify_data/ERP003634_taxonomy_abundances_SSU_v5.0.tsv',
    'Sola': 'mgnify_data/SRP237882_taxonomy_abundances_SSU_v5.0.tsv',
    'Biscay': 'mgnify_data/SRP334933_taxonomy_abundances_SSU_v5.0.tsv',
    'Baltic': 'mgnify_data/ERP140185_taxonomy_abundances_SSU_v5.0.tsv',
    'BBMO': 'mgnify_data/ERP122219_taxonomy_abundances_SSU_v5.0.tsv',
    'Svalbard': 'mgnify_data/ERP106348_taxonomy_abundances_SSU_v5.0.tsv',
    'Helgoland': 'mgnify_data/ERP144826_taxonomy_abundances_SSU_v5.0.tsv',
    'Fram': 'mgnify_data/ERP151329_taxonomy_abundances_SSU_v5.0.tsv',
}

# loading tsv tables
ds = {}
for key, value in ds_list.items():
    url = f"https://raw.githubusercontent.com/fair-ease/doc-jupyterbook-mgo/refs/heads/main/book/assets/data/{value}"
    df = pd.read_csv(url, sep="\t", header=0, index_col=0)
    # df = pd.read_csv(os.path.join(data_folder, value), sep='\t')
    df['#SampleID'] = df['#SampleID'].apply(clean_tax_row_mgnify)
    print(key, df.shape)
    ds[key] = df
```

:::{note}
:class: dropdown
loading the data from GitHub is a consequence of inability of the book executed on binder to reference "local files" using relative references.
:::


## Future prospects

We are working on merging and harmonizing the metadata tables as well, as they serve as factors for any meaningful statistical analysis.
