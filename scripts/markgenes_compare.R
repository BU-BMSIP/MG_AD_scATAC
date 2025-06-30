# ===========================================
# Marker Gene Overlap: scATAC vs scRNA
# ===========================================

# Load libraries
library(pheatmap)
library(dplyr)
library(RColorBrewer)

# -------------------------------
# Step 1: Load & Save Input Data
# -------------------------------

# Save full marker gene lists
saveRDS(rna_marker_list,  "rna_markers_list.rds")
saveRDS(atac_marker_list, "atac_markers_list.rds")
write.table(rna_markers_df,  "rna_markers_all.txt",  sep = "\t", quote = FALSE, row.names = FALSE)
write.table(atac_markers_df, "atac_markers_all.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------------
# Step 2: Clean Marker Lists
# -------------------------------

# Load from file
atac_markers <- readRDS("atac_markers_list.rds")
rna_markers  <- readRDS("rna_markers_list.rds")

# Ensure uppercase gene names
atac_markers <- lapply(atac_markers, toupper)
rna_markers  <- lapply(rna_markers, toupper)

# Filter out empty clusters
atac_markers <- atac_markers[sapply(atac_markers, length) > 0]
rna_markers  <- rna_markers[sapply(rna_markers, length) > 0]

# Save flattened cluster × gene tables
atac_df <- bind_rows(lapply(names(atac_markers), function(clu) data.frame(cluster = clu, gene = atac_markers[[clu]])))
rna_df  <- bind_rows(lapply(names(rna_markers),  function(clu) data.frame(cluster = clu, gene = rna_markers[[clu]])))

write.table(atac_df, "atac_markers.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rna_df,  "rna_markers.txt",  sep = "\t", quote = FALSE, row.names = FALSE)

# -------------------------------
# Step 3: Fisher's Exact Test
# -------------------------------

# Background: all unique genes from both modalities
background_genes <- union(unlist(atac_markers), unlist(rna_markers))

# Cluster names
atac_clusters <- names(atac_markers)
rna_clusters  <- names(rna_markers)

# Initialize result matrices
odds_ratio_mat <- matrix(0, nrow = length(atac_clusters), ncol = length(rna_clusters),
                         dimnames = list(atac_clusters, rna_clusters))
pval_mat <- matrix(1, nrow = length(atac_clusters), ncol = length(rna_clusters),
                   dimnames = list(atac_clusters, rna_clusters))

# Fisher's test helper
run_fisher_test <- function(g1, g2, background) {
  a <- length(intersect(g1, g2))
  b <- length(setdiff(g1, g2))
  c <- length(setdiff(g2, g1))
  d <- length(setdiff(background, union(g1, g2)))
  mat <- matrix(c(a, b, c, d), nrow = 2)
  tryCatch(fisher.test(mat), error = function(e) list(estimate = NA, p.value = NA))
}

# Compute overlap statistics
for (i in atac_clusters) {
  for (j in rna_clusters) {
    result <- run_fisher_test(atac_markers[[i]], rna_markers[[j]], background_genes)
    odds_ratio_mat[i, j] <- result$estimate
    pval_mat[i, j]       <- result$p.value
  }
}

# -------------------------------
# Step 4: Visualization Prep
# -------------------------------

# Log2(odds ratio), with caps
log_or_mat <- log2(odds_ratio_mat)
log_or_mat[log_or_mat > 4]  <- 4
log_or_mat[log_or_mat < -2] <- -2

# Significance star matrix
sig_mat <- matrix("", nrow = nrow(pval_mat), ncol = ncol(pval_mat))
sig_mat[pval_mat < 0.05]  <- "*"
sig_mat[pval_mat < 0.01]  <- "**"
sig_mat[pval_mat < 0.001] <- "***"
rownames(sig_mat) <- rownames(log_or_mat)
colnames(sig_mat) <- colnames(log_or_mat)

# Color palette: blue → white → red
col_fun <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))

# -------------------------------
# Step 5: Plot Heatmap
# -------------------------------

pheatmap(
  mat             = log_or_mat,
  color           = col_fun(100),
  display_numbers = sig_mat,
  number_color    = "black",
  fontsize_number = 12,
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  main            = "Marker Gene Overlap: scATAC vs scRNA (log2 Odds Ratio)",
  na_col          = "grey90",
  border_color    = "grey60",
  filename        = "marker_overlap_heatmap.png",
  width           = 8,
  height          = 6
)