# Project notebooks

---
title: "Trajectory Analysis"
---

```{r}
# Load libraries
library(Matrix)
library(monocle3)

# Load data
expr_mat <- readRDS("expression_matrix.Rds")
gene_meta <- read.csv("gene_annotation.csv", header = TRUE)
cell_meta <- read.csv("cell_metadata.csv", header = TRUE)

# Construct cell_data_set
cds <- new_cell_data_set(expr_mat,
                         cell_metadata = cell_meta,
                         gene_metadata = gene_meta)

# Preprocess and reduce dimensions
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, reduction_method = "UMAP")
cds <- cluster_cells(cds)

# Learn trajectory graph
cds <- learn_graph(cds)

# Order cells (set root cell type/cluster here)
cds <- order_cells(cds)

# Plot trajectory
plot_cells(cds, color_cells_by = "pseudotime")
```
