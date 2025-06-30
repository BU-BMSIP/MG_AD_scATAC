# ===============================================
# GO Biological Process Annotation Per Cluster
# ===============================================

# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# Load marker gene list (each element is a character vector of gene symbols)
marker_list <- readRDS("atac_markers_list.rds")

# Convert all gene symbols to uppercase
marker_list <- lapply(marker_list, toupper)

# Initialize list to store enrichment results
go_results <- list()

# Run GO enrichment for each cluster
for (clu in names(marker_list)) {
  message("Processing cluster: ", clu)

  genes <- marker_list[[clu]]
  
  # Convert gene symbols to Entrez IDs
  gene_df <- tryCatch({
    bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    warning("Conversion failed for cluster ", clu, ": ", e$message)
    return(NULL)
  })

  # If valid conversion results, perform GO enrichment
  if (!is.null(gene_df) && nrow(gene_df) > 0) {
    ego <- enrichGO(
      gene          = gene_df$ENTREZID,
      OrgDb         = org.Hs.eg.db,
      ont           = "BP",           # Biological Process
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2,
      readable      = TRUE
    )

    go_results[[clu]] <- ego
  } else {
    message("No valid genes for cluster: ", clu)
    go_results[[clu]] <- NULL
  }
}

# Create output folder if not exists
if (!dir.exists("GO_plots")) dir.create("GO_plots")

# Plot top 10 enriched GO:BP terms per cluster
for (clu in names(go_results)) {
  ego <- go_results[[clu]]

  if (!is.null(ego) && nrow(ego) > 0) {
    p <- barplot(ego, showCategory = 10, title = paste("Cluster", clu, "GO:BP"))

    ggsave(
      filename = file.path("GO_plots", paste0("Cluster_", clu, "_GO_BP.pdf")),
      plot     = p,
      width    = 8,
      height   = 6
    )
  }
}

# Save GO term tables per cluster
for (clu in names(go_results)) {
  ego <- go_results[[clu]]

  if (!is.null(ego) && nrow(ego) > 0) {
    write.table(
      as.data.frame(ego),
      file      = file.path("GO_plots", paste0("Cluster_", clu, "_GO_BP.tsv")),
      sep       = "\t",
      quote     = FALSE,
      row.names = FALSE
    )
  }
}

# Save all GO enrichment results
saveRDS(go_results, file = "GO_results_rna.rds")