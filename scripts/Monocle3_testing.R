# ===============================
# Complete Monocle3 Trajectory Analysis Pipeline
# ===============================

# --- Step 1: Install and load required packages (run once) ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(version = "3.20")

BiocManager::install(c(
  'BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
  'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
  'SummarizedExperiment', 'batchelor', 'HDF5Array',
  'ggrastr'
))

install.packages("devtools")
install.packages("remotes")
remotes::install_github("bnprks/BPCells/r")
devtools::install_github('cole-trapnell-lab/monocle3')

library(monocle3)
library(Seurat)

# --- Step 2: Convert Seurat object to Monocle3 CellDataSet (CDS) ---
# Replace 'seurat_obj' with your actual Seurat object variable name
DefaultAssay(seurat_obj) <- "RNA"  # Or "integrated" if you use integration assay

cds <- as.cell_data_set(seurat_obj)

# Transfer cluster and UMAP embeddings from Seurat to CDS
cds@clusters$UMAP$clusters <- seurat_obj$seurat_clusters
cds@int_colData@listData$reducedDims$UMAP <- Embeddings(seurat_obj, "umap")

# --- Step 3: Learn trajectory graph ---
cds <- learn_graph(cds)

# --- Step 4: Set root cells for pseudotime ordering ---
# Option 1: Interactive selection (uncomment to use)
# plot_cells(cds,
#            label_groups_by_cluster = FALSE,
#            label_leaves = TRUE,
#            label_branch_points = TRUE)
# cds <- order_cells(cds)

# Option 2: Automatic root cell assignment (recommended for automation)
root_cell <- colnames(cds)[1]  # Replace with biologically meaningful cell name if possible
cds <- order_cells(cds, root_cells = root_cell)

# --- Step 5: Plot cells colored by pseudotime ---
plot_cells(
  cds,
  color_cells_by = "pseudotime",
  label_groups_by_cluster = FALSE,
  label_leaves = TRUE,
  label_branch_points = TRUE
)

# --- Step 6: Identify genes that change over pseudotime ---
gene_fits <- fit_models(cds, model_formula_str = "~pseudotime")
fit_coefs <- coefficient_table(gene_fits)
sig_genes <- subset(fit_coefs, term == "pseudotime" & q_value < 0.05)$gene_short_name

# --- Step 7: Plot expression of top pseudotime-dependent genes ---
plot_genes_in_pseudotime(cds[sig_genes[1:6], ])