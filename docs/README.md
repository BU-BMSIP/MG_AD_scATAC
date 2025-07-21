# Project Documentation

This directory contains the full set of documentation, code, and interpretation notes for the project:

### "Understanding the Epigenomic Landscape Underlying Microglial State Transitions"

Our analysis aims to integrate **single-nucleus ATAC-seq (snATAC-seq)** and **single-nucleus RNA-seq (snRNA-seq)** data to uncover transcriptional and regulatory features associated with microglial state changes in the aging or disease-affected human brain.

---

## Contents Overview

Each file in this directory corresponds to a specific module of the analytical pipeline. These include:

### 1. Preprocessing of snATAC-seq and snRNA-seq Data
- **Scripts**: `01_preprocessing_snATAC.R`, `01_preprocessing_snRNA.R`
- Covers cell filtering, quality control (e.g. TSS enrichment, nucleosome signal), doublet removal, and gene activity matrix generation.
- Applies `ArchR` and `Seurat` pipelines respectively.
- **Outputs**: high-quality filtered cell matrices, metadata, and QC visualizations.

### 2. Dimensionality Reduction and Clustering
- **Scripts**: `02_dimensionality_reduction.R`, `03_integration_harmony.R`
- Performs LSI for snATAC-seq, PCA for snRNA-seq, and Harmony integration across batches or donors.
- Includes UMAP visualization and clustering (e.g., Louvain or Leiden algorithm).
- **Outputs**: cluster annotations, integrated embeddings, and cluster-level QC plots.

### 3. Marker Gene Identification and Cross-Modality Comparison
- **Scripts**: `04_marker_gene_identification.R`, `05_marker_overlap_analysis.R`
- Identifies differentially accessible peaks and expressed genes for each cluster.
- Compares marker genes across ATAC and RNA modalities using Fisher’s exact test.
- **Visualizations**: dot plots and heatmaps for cluster correspondence.

### 4. Peak Calling and TF–CRE–Gene Linking
- **Scripts**: `06_peak_calling_linking.R`
- Performs peak calling using MACS2 within ArchR.
- Links candidate cis-regulatory elements (CREs) to gene promoters using co-accessibility and peak-to-gene scores.
- Identifies potential transcription factor regulators using motif enrichment.

### 5. Gene Ontology (GO) and Pathway Enrichment Analysis
- **Scripts**: `07_go_enrichment_analysis.R`
- Applies `clusterProfiler` to perform GO and KEGG enrichment on gene sets from both RNA and ATAC clusters.
- **Visualizations**: bar plots, dot plots, and bubble charts.
- Helps interpret cluster-specific biological processes (e.g., immune activation, synapse pruning).

### 6. Trajectory Inference of Microglial State Transitions
- **Scripts**: `08_trajectory_analysis.R`
- Constructs pseudotime trajectories using RNA and/or ATAC data to model dynamic microglial transitions.
- Identifies key genes and regulatory features that vary along pseudotime.
- **Visualizations**: trajectory plots, heatmaps of dynamic genes, and branch-specific marker analysis.
- Supports understanding of progressive cell state changes during aging or disease.

---

## Reproducibility

- Each `.R` script is self-contained and can be executed in sequence with the provided input files.
- Key outputs (tables, plots) are saved in the corresponding `results/` or `plots/` subfolders.
- Scripts include embedded comments and section headers to guide interpretation and reuse.

---

## Notes

- Some scripts require large datasets or compute-intensive steps (e.g., MACS2 peak calling, Harmony integration); these should be run on an HPC or local server.
- Custom plotting themes are defined in `utils/plotting_functions.R`.
- For any discrepancies or data format issues, refer to the `README.md` inside each subfolder.