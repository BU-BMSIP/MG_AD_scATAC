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
library(dplyr)
library(pheatmap)

# ===== 1. Load Seurat object =====
seurat_obj <- readRDS("/projectnb/cepinet/data/scRNA/cell-2023-Sun/ROSMAP.Microglia.6regions.seurat.harmony.selected.deidentified.rds")

colnames(seurat_obj@meta.data)

p1 = DimPlot ( seurat_obj, group.by = "seurat_clusters", label = T )
ggsave ( "1.UMAP_seurat_clusters.pdf", plot = p1, width = 6.5, height =5 )

#-------------------------------------------------------------------------------------------------------
seurat_subset_non <- readRDS("/projectnb/cepinet/users/vhe/Na_Cell_2023_MG/seurat_subset_nonAD.rds")
seurat_subset_ear <- readRDS("/projectnb/cepinet/users/vhe/Na_Cell_2023_MG/seurat_subset_earlyAD.rds")
seurat_subset_lat <- readRDS("/projectnb/cepinet/users/vhe/Na_Cell_2023_MG/seurat_subset_lateAD.rds")

colnames(seurat_subset_lat@meta.data)

p1 = DimPlot ( seurat_subset_ear, group.by = "seurat_clusters", label = T )
ggsave ( "1.UMAP_seurat_clusters.pdf", plot = p1, width = 6.5, height =5 )

# ===== 2. Construct Monocle3 CellDataSet object =====
counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
cell_metadata <- seurat_obj@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))

#counts <- readMM("/projectnb/cepinet/users/Xiaojiang/monocle_TF_CRE/expression_matrix.mtx")
#gene_metadata <- read.csv("/projectnb/cepinet/users/Xiaojiang/monocle_TF_CRE/gene_annotation.csv", header = TRUE, stringsAsFactors = FALSE)
#cell_metadata <- read.csv("/projectnb/cepinet/users/Xiaojiang/monocle_TF_CRE/cell_metadata.csv", header = TRUE, stringsAsFactors = FALSE)

cds <- new_cell_data_set(
  expression_data = counts,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

#table(colData(cds)$seurat_clusters)
#cds <- cds[, colData(cds)$seurat_clusters != "12"]
#table(colData(cds)$seurat_clusters)

# ===== 3. Preprocess and use Seurat's UMAP embeddings =====
cds <- preprocess_cds(cds, method = 'PCA', num_dim = 50)

cds <- reduce_dimension(cds, reduction_method = "UMAP", preprocess_method = 'PCA')

p1 <- plot_cells(cds, reduction_method = "UMAP", color_cells_by = "seurat_clusters", show_trajectory_graph = FALSE) + ggtitle('cds.umap')

# Use Seurat UMAP embeddings
#umap_embeddings <- seurat_obj@reductions$umap@cell.embeddings
#umap_embeddings <- umap_embeddings[colnames(cds), ]
#reducedDims(cds)$UMAP <- umap_embeddings

#将seurat对象的UMAP导入
int.embed <- Embeddings(seurat_obj, reduction = "umap")
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

# ===== 1. 确保 pseudotime 已经计算 =====
colData(cds)$pseudotime <- pseudotime(cds)





library(dplyr)
library(ggplot2)

# ===== 1. 提取 pseudotime 并合并 meta =====
df <- as.data.frame(colData(cds))
df$seurat_clusters <- as.character(df$seurat_clusters)

# ===== 2. 去掉 cluster12 =====
df <- df %>% filter(seurat_clusters != "12")

# ===== 3. 去掉无效 pseudotime =====
df <- df %>% filter(is.finite(pseudotime))

# ===== 4. 计算每个 cluster 的平均 pseudotime（不分组） =====
cluster_order <- df %>%
  group_by(seurat_clusters) %>%
  summarise(mean_pseudotime = mean(pseudotime, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_pseudotime) %>%
  pull(seurat_clusters)

# ===== 5. 设置 cluster 顺序（factor levels） =====
df$seurat_clusters <- factor(df$seurat_clusters, levels = cluster_order)

# ===== 6. 绘制 boxplot（按照 pseudotime 排序后的 cluster） =====
p <- ggplot(df, aes(x = seurat_clusters, y = pseudotime, fill = ADdiag3types)) +
  geom_boxplot(outlier.size = 0.5, position = position_dodge(width = 0.8)) +
  theme_bw() +
  labs(x = "Cluster (ordered by avg pseudotime)", 
       y = "Pseudotime", 
       title = "Pseudotime by Cluster and ADdiag3types") +
  scale_fill_manual(values = c("blue", "pink", "red")) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Cluster_pseudotime_boxplot_sorted.png", plot = p, 
       width = 10, height = 6, dpi = 300)










# ===== 2. 转成 dataframe 并确保 cluster 是字符型 =====
df <- as.data.frame(colData(cds))
df$seurat_clusters <- as.character(df$seurat_clusters)

# ===== 3. 排除 cluster12 =====
df <- df[df$seurat_clusters != "12", ]

# ===== 4. 计算每个 cluster 的平均 pseudotime（去掉 Inf / NaN） =====
cluster_means <- tapply(df$pseudotime, df$seurat_clusters, function(x) {
  x <- x[is.finite(x)]        # 去掉 Inf / NaN
  mean(x, na.rm = TRUE)
})

# ===== 5. 构建 heatmap 矩阵（行列都是 cluster，值为平均 pseudotime） =====
cluster_ids <- sort(unique(df$seurat_clusters))
heat_mat <- outer(cluster_ids, cluster_ids, function(x, y) cluster_means[x])
rownames(heat_mat) <- cluster_ids
colnames(heat_mat) <- cluster_ids

# ===== 6. 绘制 heatmap 并保存 =====
p <- pheatmap(heat_mat,
              cluster_rows = FALSE,
              cluster_cols = FALSE,
              color = colorRampPalette(c("blue", "white", "red"))(100),
              main = "Cluster average pseudotime - latAD",
              silent = TRUE)

ggsave(filename = "Cluster_average_pseudotime_latAD_noCluster12.png",
       plot = grid::grid.grabExpr(grid::grid.draw(p$gtable)),
       width = 8, height = 6, dpi = 300)


library(ggplot2)

cluster_means <- tapply(df$pseudotime, df$seurat_clusters, function(x) mean(x[is.finite(x)], na.rm = TRUE))
cluster_means_df <- data.frame(cluster = names(cluster_means), avg_pseudotime = cluster_means)

# 按 pseudotime 排序（从小到大）
cluster_means_df$cluster <- factor(cluster_means_df$cluster, 
                                   levels = cluster_means_df$cluster[order(cluster_means_df$avg_pseudotime)])

ggplot(cluster_means_df, aes(x = cluster, y = avg_pseudotime, fill = avg_pseudotime)) +
  geom_bar(stat = 'identity') +
  scale_fill_gradient(low = 'blue', high = 'red') +
  theme_minimal() +
  ggtitle('Average pseudotime per cluster - latAD') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 旋转X轴标签方便看

ggsave("Cluster_avg_pseudotime_bar_earlyAD_non12_sorted.png", width = 8, height = 6, dpi = 300)