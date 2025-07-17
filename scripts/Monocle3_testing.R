# Step 1: Install and load packages (run once)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.20")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'ggrastr'))

install.packages("devtools")
install.packages("remotes")
remotes::install_github("bnprks/BPCells/r")
devtools::install_github('cole-trapnell-lab/monocle3')

library(monocle3)
library(Seurat)

# Step 2: Convert Seurat to Monocle3 CDS
DefaultAssay(seurat_obj) <- "RNA"  # Or "integrated", if you're using integration
cds <- as.cell_data_set(seurat_obj)
cds@clusters$UMAP$clusters <- seurat_obj$seurat_clusters
cds@int_colData@listData$reducedDims$UMAP <- Embeddings(seurat_obj, "umap")

# Step 3: Learn trajectory graph
cds <- learn_graph(cds)

# Step 4: Define pseudotime root (choose one option)

# Option A: Manually select root cell
# plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = TRUE, label_branch_points = TRUE)
# cds <- order_cells(cds)

# Option B: Set root by cell name (recommended for automation)
root_cell <- colnames(cds)[1]  # You can replace with a biologically meaningful cell ID
cds <- order_cells(cds, root_cells = root_cell)

# Step 5: Visualize pseudotime
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_leaves = TRUE,
           label_branch_points = TRUE)

# Step 6: Find pseudotime-dependent genes
gene_fits <- fit_models(cds, model_formula_str = "~pseudotime")
fit_coefs <- coefficient_table(gene_fits)
sig_genes <- subset(fit_coefs, term == "pseudotime" & q_value < 0.05)$gene_short_name

# Step 7: Plot gene expression across pseudotime
plot_genes_in_pseudotime(cds[sig_genes[1:6], ])