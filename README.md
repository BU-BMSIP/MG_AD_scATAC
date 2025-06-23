# Microglial State Transition Integration Pipeline

A pipeline for integrating snATAC-seq and snRNA-seq data to uncover regulatory mechanisms underlying microglial state transitions in Alzheimer’s disease.

## Project Overview

Microglia, the brain's resident immune cells, adopt diverse transcriptional states during aging and Alzheimer's disease (AD) progression. This project aims to investigate the epigenomic regulation underlying these transitions by:

- Performing dimensionality reduction and clustering on snATAC-seq data using ArchR
- Comparing marker genes derived from gene scores (ATAC-seq) and RNA expression (RNA-seq)
- Linking transcription factors (TFs), cis-regulatory elements (CREs), and genes via peak-to-gene co-accessibility
- Visualizing UMAP plots, differential peaks, and regulatory heatmaps to interpret regulatory programs

## Tools and Technologies

- **ArchR** for snATAC-seq preprocessing, LSI/Harmony integration, clustering, and gene score analysis
- **Seurat** for snRNA-seq processing
- **MACS2** for peak calling
- **Fisher’s Exact Test** for marker gene overlap heatmaps
- **GO enrichment** for functional annotation of ATAC-seq marker genes

## Related References

- Sun et al., *Cell*, 2023 – snRNA-seq of human microglia in AD  
- Xiong et al., *Cell*, 2023 – Epigenomic analysis and TF–CRE–gene linking  
- ArchR documentation – https://www.archrproject.com/bookdown/