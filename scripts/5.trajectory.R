# Set your personal library path
.libPaths("/projectnb/cepinet/users/vhe/R_4.4.0_libs_monocle3")

# Install managers
install.packages("BiocManager", repos = "https://cloud.r-project.org")
install.packages("devtools")
install.packages("remotes")

# Required Bioconductor packages
BiocManager::install(c("SingleCellExperiment", "BiocGenerics", "DelayedArray", "matrixStats"))

# General dependencies
install.packages(c("Matrix", "ggplot2", "RcppEigen", "igraph", "viridis", "rcpp"))

# Now install monocle3
remotes::install_github("cole-trapnell-lab/monocle3")

# ===== 0. 加载包 =====
library(Seurat)
library(monocle3)
library(SingleCellExperiment)
library(ggplot2)

# ===== 1. 读取 Seurat 对象 =====
seurat_obj <- readRDS("ROSMAP.Microglia.6regions.seurat.harmony.selected.deidentified.rds")

# ===== 2. 构建 Monocle3 CellDataSet 对象 =====
counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
cell_metadata <- seurat_obj@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))

cds <- new_cell_data_set(
  expression_data = counts,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

# ===== 3. 预处理 + 用 Seurat 的 UMAP 结果 =====
cds <- preprocess_cds(cds, num_dim = 50)

# 3. 用 Seurat 的 UMAP 嵌入
umap_embeddings <- seurat_obj@reductions$umap@cell.embeddings
umap_embeddings <- umap_embeddings[colnames(cds), ]
reducedDims(cds)$UMAP <- umap_embeddings

# 4. 用 Seurat 聚类（可选）
cds@clusters$UMAP$clusters <- as.character(seurat_obj@meta.data$seurat_clusters)

# 🔧 插入这一步！！构建邻接图（必须的）
cds <- cluster_cells(cds, reduction_method = "UMAP")

# 5. 构建 trajectory graph
cds <- learn_graph(cds)

# ===== 6. 指定起始 cluster（比如 cluster 0）并排序 pseudotime =====
get_earliest_principal_node <- function(cds, cluster){
  cell_ids <- colnames(cds)[cds@clusters$UMAP$clusters == cluster]
  closest_vertex <- cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  tab <- table(closest_vertex[cell_ids, ])
  principal_node <- names(tab)[which.max(tab)]
  return(principal_node)
}

cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds, cluster = "0"))

# ===== 7. 可视化 pseudotime =====
p <- plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = TRUE,
           label_leaves = TRUE,
           label_branch_points = TRUE)

ggsave("pseudotime_w_cluster.png", plot = p, width = 10, height = 8, dpi = 300)