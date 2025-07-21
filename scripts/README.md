# Project Scripts

This directory contains all scripts used in the analysis of microglial state transitions using single-nucleus multiomic data. The scripts are organized to reflect the computational pipeline described in the original publication.

### Original Paper  
**Sun, Victor et al., Cell, 2023**  
[Link to paper](https://www.cell.com/cell/fulltext/S0092-8674(23)00971-6)

### Data Download  
Raw and processed data are available from:  
🔗 [https://personal.broadinstitute.org/cboix/sun_victor_et_al_data/](https://personal.broadinstitute.org/cboix/sun_victor_et_al_data/)

- ATAC fragment files were obtained directly from the Na et al. dataset.  
- Peak-to-gene linking results were extracted from the supplementary materials of the Cell publication:  
  [Supplementary Data](https://www.cell.com/cell/fulltext/S0092-8674(23)00971-6#mmc1)

---

## Methods Overview

### ATAC-seq Data Processing

We followed the computational workflow used in the accompanying manuscript by **Xiong et al.**:

1. **Preprocessing**
   - FASTQ files were generated via demultiplexing using `cellranger-atac` v1.1.0.
   - Reads were aligned to the GRCh38 human genome with `cellranger-atac count` to produce fragment files.

2. **Quality Control & Filtering**
   - Analysis was performed using **ArchR** v1.0.1.
   - Doublets were filtered using `filterDoublets`.
   - Cells were retained if they met both:
     - TSS enrichment score > 6
     - Number of fragments between 1,000 and 100,000

3. **Dimensionality Reduction & Clustering**
   - Iterative LSI was performed on 500 bp tiles with parameters:
     - `iterations = 4`
     - `resolution = 0.2`
     - `varFeatures = 50,000`
   - UMAP was used for visualization.
   - Gene scores were computed using ArchR.
   - Cell types were annotated using well-known brain markers.

4. **Microglia Subsetting**
   - Only clusters annotated as **microglia/immune cells** were retained for downstream analysis.

---

### Full Dataset

1. **Cluster Annotation**
   - Cell clusters were annotated using canonical markers for major brain cell types:
     - Excitatory/Inhibitory neurons
     - Astrocytes
     - Oligodendrocytes
     - OPCs
     - Microglia
     - Vascular cells
   - Annotation was validated using enrichment analysis of a broader set of literature-derived markers.

2. **Microglia/Immune Cell Selection**
   A cell was retained as microglia/immune if **all three** of the following were true:
   - Belonged to a cluster annotated as microglia/immune;
   - Had the highest cell type score for microglia/immune;
   - Microglia/immune score was **≥ 2×** higher than the second highest score.

3. **Downstream Processing**
   - Selected microglia/immune cells were re-processed using the same dimensionality reduction and clustering pipeline.
   - Differentially expressed genes were identified using **Wilcoxon rank-sum test** in Seurat with:
     - `min.pct = 0.25`
     - `logfc.threshold = 0.25`
=