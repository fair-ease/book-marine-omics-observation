---
title: Virtual Research Environment
authors:
  - F-E
---

## VRE

:::{warning}
Expected deploymnet at BC 2026 is expected during autumn 2025.
:::

## General Principles

The VRE uses UDAL to search and sub-set the data based on the RDF triple store. The IDDAS/UDAL provides access to data sets, and subsets, or co-locates data from other data spaces into VRE space. The VRE contains a number of Jupyter Notebooks that provide default, continuously integrated, data representations and analyses (e.g. alpha and beta diversity analyses between taxonomic and functional attributes of samples), plus a number of workflows demonstrating potential analyses that can be conducted on the EMO BON data (e.g. biosynthetic gene cluster analyses).

Some analyses in Jupyter notebooks use the Galaxy API to run on the Galaxy back-end that automatically designates a specific computing node within the Pulsar network that is local to the data.
