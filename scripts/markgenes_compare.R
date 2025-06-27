# ===========================================
# Marker Gene Overlap: ATAC vs RNA
# ===========================================

# Prerequisites
library(pheatmap)
library(dplyr)
library(RColorBrewer)

# Save marker gene lists and combined data frames
saveRDS(rna_marker_list, file = "rna_markers_list.rds")
write.table(rna_markers_df, file = "rna_markers_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)

saveRDS(atac_marker_list, file = "atac_markers_list.rds")
write.table(atac_markers_df, file = "atac_markers_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Load precomputed marker lists
atac_markers <- readRDS("atac_markers_list.rds")  # List: cluster → gene vector
rna_markers  <- readRDS("rna_markers_list.rds")   # List: cluster → gene vector

# Filter out empty ATAC clusters
atac_markers_filt <- atac_markers[sapply(atac_markers, length) > 0]
atac_df <- do.call(rbind, lapply(names(atac_markers_filt), function(clu) {
  data.frame(cluster = clu, gene = atac_markers_filt[[clu]])
}))
write.table(atac_df, file = "atac_markers.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Filter out empty RNA clusters
rna_markers_filt <- rna_markers[sapply(rna_markers, length) > 0]
rna_df <- do.call(rbind, lapply(names(rna_markers_filt), function(clu) {
  data.frame(cluster = clu, gene = rna_markers_filt[[clu]])
}))
write.table(rna_df, file = "rna_markers.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Ensure uppercase gene symbols for consistency
atac_markers <- lapply(atac_markers, toupper)
rna_markers  <- lapply(rna_markers, toupper)

# Create universe of genes (background set for Fisher test)
background_genes <- union(unlist(atac_markers), unlist(rna_markers))

# Prepare output matrices
atac_clusters <- names(atac_markers)
rna_clusters  <- names(rna_markers)

odds_ratio_mat <- matrix(0, nrow = length(atac_clusters), ncol = length(rna_clusters),
                         dimnames = list(atac_clusters, rna_clusters))
pval_mat <- matrix(1, nrow = length(atac_clusters), ncol = length(rna_clusters),
                   dimnames = list(atac_clusters, rna_clusters))

# Loop over all ATAC × RNA cluster pairs
for (i in atac_clusters) {
  for (j in rna_clusters) {
    genes_i <- atac_markers[[i]]
    genes_j <- rna_markers[[j]]

    # Construct contingency table for Fisher's Exact Test
    a <- length(intersect(genes_i, genes_j))                         # Overlap count
    b <- length(setdiff(genes_i, genes_j))                           # ATAC-only genes
    c <- length(setdiff(genes_j, genes_i))                           # RNA-only genes
    d <- length(setdiff(background_genes, union(genes_i, genes_j)))  # Genes in neither

    mat <- matrix(c(a, b, c, d), nrow = 2)

    # Perform Fisher's Exact Test
    test <- fisher.test(mat)
    odds_ratio_mat[i, j] <- test$estimate
    pval_mat[i, j]       <- test$p.value
  }
}

# Log-transform odds ratio matrix for visualization
log_or_mat <- log2(odds_ratio_mat)

# Cap extreme values to improve color scaling
log_or_mat[log_or_mat > 4]  <- 4
log_or_mat[log_or_mat < -2] <- -2

# Custom color scale: Blue → White → Red
col_fun <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))

# Generate significance star matrix (* / ** / ***)
sig_mat <- matrix("", nrow = nrow(pval_mat), ncol = ncol(pval_mat))
sig_mat[pval_mat < 0.05]  <- "*"
sig_mat[pval_mat < 0.01]  <- "**"
sig_mat[pval_mat < 0.001] <- "***"
rownames(sig_mat) <- rownames(log_or_mat)
colnames(sig_mat) <- colnames(log_or_mat)

# Plot heatmap
pheatmap(
  mat             = log_or_mat,
  color           = col_fun(100),
  display_numbers = sig_mat,          # Show significance stars
  number_color    = "black",
  fontsize_number = 12,
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  main            = "Marker Gene Overlap: scATAC vs scRNA (log2 Odds Ratio)",
  na_col          = "grey90",
  border_color    = "grey60",
  filename        = "marker_overlap_heatmap.png"
)