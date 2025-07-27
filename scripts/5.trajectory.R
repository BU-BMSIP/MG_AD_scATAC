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
seurat_obj <- readRDS("ROSMAP.Microglia.6regions.seurat.harmony.selected.deidentified.rds")

# ===== 2. Construct Monocle3 CellDataSet object =====
counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")
cell_metadata <- seurat_obj@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))

cds <- new_cell_data_set(
  expression_data = counts,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

# ===== 3. Preprocess and use Seurat's UMAP embeddings =====
cds <- preprocess_cds(cds, num_dim = 50)

# Use Seurat UMAP embeddings
umap_embeddings <- seurat_obj@reductions$umap@cell.embeddings
umap_embeddings <- umap_embeddings[colnames(cds), ]
reducedDims(cds)$UMAP <- umap_embeddings

# ===== 4. Use Seurat clustering (optional) =====
cds@clusters$UMAP$clusters <- as.character(seurat_obj@meta.data$seurat_clusters)

# ===== 5. Build cell clusters graph (required step) =====
cds <- cluster_cells(cds, reduction_method = "UMAP")

# ===== 6. Learn trajectory graph =====
cds <- learn_graph(cds)

# ===== 7. Define root cluster and order cells by pseudotime =====
get_earliest_principal_node <- function(cds, cluster){
  cell_ids <- colnames(cds)[cds@clusters$UMAP$clusters == cluster]
  closest_vertex <- cds@principal_graph_aux$UMAP$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  tab <- table(closest_vertex[cell_ids, ])
  principal_node <- names(tab)[which.max(tab)]
  return(principal_node)
}

cds <- order_cells(cds, root_pr_nodes = get_earliest_principal_node(cds, cluster = "0"))

# ===== 8. Visualize pseudotime =====
p <- plot_cells(cds,
                color_cells_by = "pseudotime",
                label_groups_by_cluster = TRUE,
                label_leaves = TRUE,
                label_branch_points = TRUE)

# Plot expression of selected genes without trajectory graph
p1 <- plot_cells(
  cds,
  genes = c("CX3CR1", "PPARG", "ITGA4"),
  label_cell_groups = FALSE,
  show_trajectory_graph = FALSE,
  min_expr = 1
)

# Save the gene expression plot
ggsave("pseudotime_three_genes.png", plot = p1, width = 10, height = 8, dpi = 300)