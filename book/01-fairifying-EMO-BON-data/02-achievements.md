---
title: Deliverables
authors:
  - David Palecek, Cymon Cox
---

## Delivered

- Definition of an ontology and namespace describing the EMO BON data model
- Development of sampling event (meta)data validation (Pydantic library framework) and quality control procedures using explicit rules and GitHub actions
- Procedures and software creation to perform semantic uplift of EMO BON data using the data model to RDF triples
- `RO-crate` specification and generation software for the metaGOflow data products
- Jupyter notebooks of data visualisations of taxonomic and genetic `alpha` and `beta` diversities
- PyPI distributed [marine-omics-methods](https://github.com/emo-bon/marine-omics-methods) package of reusable methods implemented for the Jupyter notebooks
- Jupyter notebook (Python) workflow for [biosynthetic gene cluster](https://github.com/emo-bon/momics-demos/blob/main/README.md#wf3-biosynthetic-gene-clusters-bgcs) identification using software on the Galaxy back-end via their API
Galaxy workflow execution using newly installed Galaxy Pulsar node
- Execute workflow via OGC:WPS (on top of Galaxy specific API)
- Metadata RDF triple generation
- DCAT catalogue templating descriptions of EMO BON data assets, i.e. IDDAS ingestion framework for the metadata triple store
- Storage solution for MGF data products - S3 Object Store in testbed infrastructure
- IDDAS ingestion framework for the metadata triple database and the MGF data product RO-crates
- RO-crate viewer to browse and download specific files from the sample and MGF RO-crates
- Incorporate our RO-crates and triple store into the UDAL code, provide named queries
- Initial VRE construction on the BlueCloud platform

## To be delivered

- SPARQL endpoint construction for the RDF triples of sample metadata and metaGOflow data products
- Documentation and landing pages from the VRE on BlueCloud2026 infrastructure.
- Build FAIR metadata records for all our datasets/triple store and have them assessed through tools like F-UJI
- Ensure provenance is provided out of our Jupyter Notebooks in the VRE
