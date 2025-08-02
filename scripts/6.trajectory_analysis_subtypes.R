# Set your personal library path
.libPaths("/projectnb/cepinet/users/vhe/R_4.4.0_libs_monocle3")

# Install package managers
install.packages("BiocManager", repos = "https://cloud.r-project.org")
install.packages("devtools")
install.packages("remotes")

# Install required Bioconductor packages
BiocManager::install(c("SingleCellExperiment", "BiocGenerics", "DelayedArray", "matrixStats"))

# Install general dependencies
install.packages(c("Matrix", "ggplot2", "RcppEigen", "igraph", "viridis", "rcpp"))

# Install monocle3 from GitHub
remotes::install_github("cole-trapnell-lab/monocle3")

# ===== 0. Load libraries =====
library(Seurat)
library(monocle3)
library(SingleCellExperiment)
library(ggplot2)

# ===== 1. Load Seurat object =====
seurat_obj <- readRDS("/projectnb/cepinet/data/scRNA/cell-2023-Sun/ROSMAP.Microglia.6regions.seurat.harmony.selected.deidentified.rds")

colnames(seurat_obj@meta.data)

p1 = DimPlot ( seurat_obj, group.by = "seurat_clusters", label = T )
ggsave ( "1.UMAP_seurat_clusters.pdf", plot = p1, width = 6.5, height =5 )

#-------------------------------------------------------------------------------------------------------
seurat_subset_non <- readRDS("/projectnb/cepinet/users/vhe/Na_Cell_2023_MG/seurat_subset_nonAD.rds")
seurat_subset_ear <- readRDS("/projectnb/cepinet/users/vhe/Na_Cell_2023_MG/seurat_subset_earlyAD.rds")
seurat_subset_lat <- readRDS("/projectnb/cepinet/users/vhe/Na_Cell_2023_MG/seurat_subset_lateAD.rds")

colnames(seurat_subset_non@meta.data)

p1 = DimPlot ( seurat_subset_non, group.by = "seurat_clusters", label = T )
ggsave ( "1.UMAP_seurat_clusters.pdf", plot = p1, width = 6.5, height =5 )

# ===== 2. Construct Monocle3 CellDataSet object =====
counts <- GetAssayData(seurat_subset_non, assay = "RNA", slot = "counts")
cell_metadata <- seurat_subset_non@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))

cds <- new_cell_data_set(
  expression_data = counts,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

# ===== 3. Preprocess and use Seurat's UMAP embeddings =====
cds <- preprocess_cds(cds, method = 'PCA', num_dim = 50)

cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = 'PCA')

p1 <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = "seurat_clusters", show_trajectory_graph = FALSE) + ggtitle('cds.umap')

# Use Seurat UMAP embeddings
#umap_embeddings <- seurat_obj@reductions$umap@cell.embeddings
#umap_embeddings <- umap_embeddings[colnames(cds), ]
#reducedDims(cds)$UMAP <- umap_embeddings

#将seurat对象的UMAP导入
int.embed <- Embeddings(seurat_subset_non, reduction = "umap")
#排序
int.embed <- int.embed [rownames(cds@int_colData$reducedDims$UMAP),]
#导入
cds@int_colData$reducedDims$UMAP <- int.embed
#画图
p2 <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = "seurat_clusters", show_trajectory_graph = FALSE) + ggtitle('seurat.umap')
p =p1|p2
ggsave("2.Reduction_Compare.pdf",plot = p, width = 10, height = 5)


# ===== 4. Use Seurat clustering (optional) =====
#cds@clusters$UMAP$clusters <- as.character(seurat_obj@meta.data$seurat_clusters)

# ===== 5. Build cell clusters graph (required step) =====
cds <- cluster_cells(cds, reduction_method = "UMAP")

p1 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + ggtitle("partition")
ggsave("3.cluster_Partition.pdf",plot = p1, width = 6, height = 5)

# ===== 6. Learn trajectory graph =====
#cds <- learn_graph(cds)
cds <- learn_graph(cds, learn_graph_control = list(euclidean_distance_ratio =0.8))

p = plot_cells(cds, color_cells_by = "partition", label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
ggsave("4.Trajectory.pdf", plot = p, width = 6, height = 5)

# ===== 7. Define root cluster and order cells by pseudotime =====
cds <- order_cells(cds)

p = plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE)
ggsave("5.Trajectory_Pseudotime.pdf", plot = p, width = 8, height = 6)