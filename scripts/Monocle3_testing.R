# The basic pipeline of Monocle3

# Installation
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'ggrastr'))

install.packages("devtools")
remotes::install_github("bnprks/BPCells/r")

devtools::install_github('cole-trapnell-lab/monocle3')

library(monocle3)
library(Seurat)