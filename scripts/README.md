# Project scripts

All project scripts are in this directory

orignal paper https://www.cell.com/cell/fulltext/S0092-8674(23)00971-6

download from https://personal.broadinstitute.org/cboix/sun_victor_et_al_data/

Fragments from Na direclty

peak-to-gene Linking from supplementary table, https://www.cell.com/cell/fulltext/S0092-8674(23)00971-6?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS0092867423009716%3Fshowall%3Dtrue#mmc1


## method

We processed snATAC-seq data using the same computational pipeline as Xiong et al.21 accompanied manuscript. Specifically, we firstly generated raw data of FASTQ for each sample by demultiplexing the reads with cellranger-atac software(v1.1.0),70 and then mapped the reads to human reference genome version GRCh38 using “cellranger-atac count” to obtain the fragment file for each sample. We processed the snATAC-seq data using ArchR (v1.0.1).57 We removed the potential doublets using the “filterDoublets” function in ArchR. We kept the cells with TSS enrichment more than 6 and the number of fragments between 1000 and 100,000 for further analysis. We performed Iterative LSI dimension reduction and clustering using the matrix of 500 bp tile-based, with parameters “iterations = 4, resolution = 0.2, varFeat = 50000”. The UMAP was used to visualize cell embedding for all cells. We generated the gene score matrix using ArchR, and annotated the cell type for each cluster based on the gene score of well-known markers in the brain.11 We integrated the clusters that were annotated as microglia/immune cells for further analysis.



For the full datasets with all cell types (2.8 million cells), we first annotated the cell type for each cluster based on three widely-used canonical markers of major cell types in the brain (including excitatory and inhibitory neurons, astrocytes, oligodendrocytes, OPCs, microglia and vascular cells)11 and a list of markers for immune cells.11,71 We also tested the enrichment of a large set of markers72 in highly expressed genes for each cluster to confirm the annotation based on several marker genes. We next calculated the cell type scores (i.e., astrocyte, oligodendrocyte, microglia, etc) for each cell, which were represented by the average expression of a group of markers for each cell type.72 The cells were then selected as microglia/immune cells for further integrative analysis if and only if (1) the clusters that the cells belong to were annotated as microglia/immune cells; and (2) the cells had the highest score for microglia/immune cells, and 3) the score for microglia/immune cells was 2-fold higher than the second highest score. For the selected microglia/immune cells, we followed the same pipeline to perform dimensional reduction and clustering with the same parameters as full datasets. We used the Wilcoxon rank-sum test in Seurat with customized parameters (min.pct = 0.25, logfc.threshold = 0.25) to identify highly expressed genes for each cluster compared to all cells from other clusters.