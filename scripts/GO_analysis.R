# ===============================================
# GO Biological Process Annotation Per Cluster
# ===============================================

# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# Load marker gene list
marker_list <- readRDS("atac_markers_list.rds")

# Ensure gene symbols are uppercase
marker_list <- lapply(marker_list, toupper)

# Run GO enrichment per cluster
go_results <- list()

for (clu in names(marker_list)) {
  genes <- marker_list[[clu]]

  # SYMBOL → ENTREZ
  gene.df <- bitr(
    genes,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )

  if (!is.null(gene.df) && nrow(gene.df) > 0) {
    ego <- enrichGO(
      gene          = gene.df$ENTREZID,
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",           # Biological Process
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    )

    go_results[[clu]] <- ego
  } else {
    go_results[[clu]] <- NULL
  }
}

# Create output folder
dir.create("GO_plots", showWarnings = FALSE)

# Plot and save GO barplots (Top 10 terms per cluster)
for (clu in names(go_results)) {
  ego <- go_results[[clu]]
  
  if (!is.null(ego) && nrow(ego) > 0) {
    p <- barplot(ego, showCategory = 10, title = paste("Cluster", clu, "GO:BP"))

    ggsave(
      filename = paste0("GO_plots/Cluster_", clu, "_GO_BP.pdf"),
      plot     = p,
      width    = 8,
      height   = 6
    )
  }
}

# Save full GO tables per cluster
for (clu in names(go_results)) {
  ego <- go_results[[clu]]

  if (!is.null(ego) && nrow(ego) > 0) {
    write.table(
      as.data.frame(ego),
      file      = paste0("GO_plots/Cluster_", clu, "_GO_BP.txt"),
      sep       = "\t",
      quote     = FALSE,
      row.names = FALSE
    )
  }
}

# Save entire result list for future reuse
saveRDS(go_results, file = "GO_results_rna.rds")