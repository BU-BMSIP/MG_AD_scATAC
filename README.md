# Microglial State Transition Integration Pipeline

Understanding the Epigenomic Landscape Underlying Microglial State Transitions: A pipeline for integrating snATAC-seq and snRNA-seq data to uncover regulatory mechanisms underlying microglial state transitions in Alzheimer’s disease.

## Project Overview

Microglia, the brain's resident immune cells, adopt diverse transcriptional states during aging and Alzheimer's disease (AD) progression. This project aims to investigate the epigenomic regulation underlying these transitions by:

- Performing dimensionality reduction and clustering on snATAC-seq data using ArchR
- Comparing marker genes derived from gene scores (ATAC-seq) and RNA expression (RNA-seq)
- Linking transcription factors (TFs), cis-regulatory elements (CREs), and genes via peak-to-gene co-accessibility
- Visualizing UMAP plots, differential peaks, and regulatory heatmaps to interpret regulatory programs

## Tools and Technologies

| Tool         | Purpose                                                   |
|--------------|-----------------------------------------------------------|
| **ArchR**    | snATAC-seq preprocessing, LSI/Harmony integration, clustering |
| **Seurat**   | snRNA-seq clustering and marker detection                 |
| **clusterProfiler** | GO term enrichment and annotation of marker genes        |
| **ggplot2**  | Custom visualizations (bubble plots, heatmaps, UMAPs)     |
| **Fisher's Exact Test** | Marker gene overlap between ATAC and RNA clusters    |
| **MACS2**    | Peak calling for ATAC data (external preprocessing)       |


## Analysis Workflow

### 1. **Preprocessing & Clustering**
- Load raw ATAC ArchR project
- Filter out suspected doublets (e.g., C1/C8)
- Apply Iterative LSI and Harmony batch correction
- Perform clustering on Harmony-reduced dimensions
- Generate UMAPs for sample, region, cluster visualization

### 2. **Marker Gene Analysis**
- Use `getMarkerFeatures()` on the GeneScoreMatrix
- Apply relaxed thresholds for sparse clusters:
  - Default: `FDR <= 0.01 & Log2FC >= 1.25`
  - Relaxed: `FDR <= 0.1 & Log2FC >= 0.5`
- Save `markerList.rds` and print the number of markers per cluster

### 3. **Gene Enrichment Analysis**
- Combine all markers, rank by `Log2FC`, select top 100 genes
- Convert to Entrez ID using `bitr()`
- Run `enrichGO()` (ontology: BP)
- Visualize top terms using ggplot2 bubble plot:
  - X-axis: GeneRatio
  - Y-axis: GO term
  - Size: Gene count
  - Color: Adjusted p-value (log scale)

### 4. **ATAC-RNA Marker Overlap**
- Use RNA marker list from Seurat clustering
- Perform Fisher’s exact test across ATAC vs RNA clusters
- Generate odds ratio heatmap with significance asterisks
- Standardize background gene set for fair comparison

### 5. **Modality Alignment**
- Calculate proportion of nearest neighbors from RNA to ATAC clusters (C1–C5)
- Evaluate modality concordance after Harmony integration

### 6. **Visualization Outputs**
- Gene score heatmaps (top marker genes, microglia markers)
- UMAP gene score plots (with/without imputation)
- Browser track plots of selected genes
- GO enrichment bubble chart
- Overlap heatmap (ATAC vs RNA)


## Related References

- Sun et al., *Cell*, 2023 – snRNA-seq of human microglia in AD  
- Xiong et al., *Cell*, 2023 – Epigenomic analysis and TF–CRE–gene linking  
- ArchR documentation – https://www.archrproject.com/bookdown/
