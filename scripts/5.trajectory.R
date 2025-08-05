# ===== 0. Set up library paths and install dependencies =====
.libPaths("/projectnb/cepinet/users/vhe/R_4.4.0_libs_monocle3")

# Install package managers
install.packages("BiocManager", repos = "https://cloud.r-project.org")
install.packages(c("devtools", "remotes"))

# Install required Bioconductor packages
BiocManager::install(c("SingleCellExperiment", "BiocGenerics", "DelayedArray", "matrixStats"))

# Install general dependencies
install.packages(c("Matrix", "ggplot2", "RcppEigen", "igraph", "viridis", "rcpp"))

# Install Monocle3 from GitHub
remotes::install_github("cole-trapnell-lab/monocle3")

# ===== 1. Load libraries =====
library(Seurat)
library(monocle3)
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)  # For combining plots

# ===== 2. Load Seurat object =====
seurat_obj <- readRDS("/projectnb/cepinet/data/scRNA/cell-2023-Sun/ROSMAP.Microglia.6regions.seurat.harmony.selected.deidentified.rds")

# Visualize Seurat clustering
p1 <- DimPlot(seurat_obj, group.by = "seurat_clusters", label = TRUE)
ggsave("1.UMAP_seurat_clusters.pdf", plot = p1, width = 6.5, height = 5)

# ===== 3. Construct Monocle3 CellDataSet =====
counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
cell_metadata <- seurat_obj@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))

cds <- new_cell_data_set(
  expression_data = counts,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

# ===== 4. Preprocessing and dimensionality reduction =====
cds <- preprocess_cds(cds, method = 'PCA', num_dim = 50)
cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = 'PCA')

# UMAP from Monocle3
p1 <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = "seurat_clusters", show_trajectory_graph = FALSE) + ggtitle("Monocle3 UMAP")

# ===== 5. Import Seurat UMAP into Monocle3 =====
seurat_umap <- Embeddings(seurat_obj, reduction = "umap")
seurat_umap <- seurat_umap[rownames(cds@int_colData$reducedDims$UMAP), ]
cds@int_colData$reducedDims$UMAP <- seurat_umap

# UMAP using imported Seurat embedding
p2 <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = "seurat_clusters", show_trajectory_graph = FALSE) + ggtitle("Seurat UMAP")

# Compare the two embeddings
p <- p1 | p2
ggsave("2.Reduction_Compare.pdf", plot = p, width = 10, height = 5)

# ===== 6. Clustering and partitioning =====
cds <- cluster_cells(cds, reduction_method = "UMAP")
p1 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + ggtitle("Partitioning")
ggsave("3.cluster_Partition.pdf", plot = p1, width = 6, height = 5)

# ===== 7. Learn trajectory graph =====
cds <- learn_graph(cds, learn_graph_control = list(euclidean_distance_ratio = 0.8))
p <- plot_cells(cds, color_cells_by = "partition", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
ggsave("4.Trajectory.pdf", plot = p, width = 6, height = 5)

# ===== 8. Define root and order cells by pseudotime =====
cds <- order_cells(cds)
p <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
ggsave("5.Trajectory_Pseudotime.pdf", plot = p, width = 8, height = 6)