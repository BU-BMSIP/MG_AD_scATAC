##-------------try to print the top 10 genes for each clusters that don't have much marker genes------------

markersGS <- readRDS("brain.microglia.filter/markersGS_7.2.rds")
markerList <- getMarkers(markersGS)

target_clusters <- c("C1", "C6")
top_n <- 100
top_markers_list <- list()

for (cluster in target_clusters) {
  df <- markerList[[cluster]]
  if (is.null(df)) {
    cat("Cluster", cluster, "is NULL. Skipped.\n")
    next
  }
  
  # 自动选择排序列
  if ("FDR" %in% colnames(df)) {
    sort_col <- "FDR"
  } else if ("Pval" %in% colnames(df)) {
    sort_col <- "Pval"
  } else {
    stop("No FDR or Pval found in", cluster)
  }
  
  df_sorted <- df[order(df[[sort_col]]), ]
  top_genes <- head(df_sorted, min(nrow(df_sorted), top_n))
  top_markers_list[[cluster]] <- top_genes
  
  cat("Top genes for", cluster, ":\n")
  print(head(top_genes[, c("name", sort_col)], 10))
}

##-------------how to see all the genes-------------

for (cluster in names(markerList)) {
  df <- markerList[[cluster]]
  if (is.null(df)) {
    cat("Cluster", cluster, "is NULL. Skipped.\n")
    next
  }
  
  # Get top 10 genes (no sorting)
  top_genes <- head(df$name, 10)
  
  cat("\nTop 10 genes for", cluster, ":\n")
  print(top_genes)
}