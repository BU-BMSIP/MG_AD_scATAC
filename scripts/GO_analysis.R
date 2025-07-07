# Setting up environment ===================================================

# Clean environment
rm(list = ls(all.names = TRUE)) # Clear all objects, including hidden ones
gc() # Free up memory and report usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # Avoid truncated output and scientific notation

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# Load markerList
markerList <- readRDS("brain.microglia.filter/markerList_7.2.rds")

# Convert all gene names to uppercase
markerList <- lapply(markerList, function(df) {
  df$name <- toupper(df$name)
  df
})

# Create directory to save results
dir.create("PEA/Results", recursive = TRUE, showWarnings = FALSE)

# Define background gene set and convert to ENTREZID
background_genes <- unique(unlist(lapply(markerList, function(df) df$name)))
background_entrez <- bitr(background_genes,
                          fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)

mapped_genes <- background_entrez$SYMBOL
unmapped_genes <- setdiff(background_genes, mapped_genes)
cat("Unmapped background genes:", length(unmapped_genes), "\n")

# Loop through each cluster
for (cluster in names(markerList)) {
  
  cat("\nProcessing cluster:", cluster, "\n")
  
  # Top 100 marker genes
  top_genes <- head(markerList[[cluster]]$name, 100)
  
  if (length(top_genes) == 0 || all(is.na(top_genes))) {
    cat("No genes available for cluster", cluster, "- skipping.\n")
    next
  }
  
  # Convert SYMBOL to ENTREZID
  gene_entrez <- bitr(top_genes,
                      fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)
  
  if (nrow(gene_entrez) == 0) {
    cat("No valid Entrez IDs for cluster", cluster, "- skipping.\n")
    next
  }
  
  # GO BP enrichment analysis
  ego <- enrichGO(
    gene = gene_entrez$ENTREZID,
    universe = background_entrez$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    keyType = "ENTREZID",
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    readable = TRUE
  )
  
  # Simplify redundant GO terms
  ego_simplified <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
  
  # Plotting logic: use simplified result if available, otherwise use original
  if (!is.null(ego_simplified) && nrow(ego_simplified) > 0) {
    plot_obj <- dotplot(ego_simplified, showCategory = 20) + ggtitle(paste("Cluster", cluster, "- GO BP"))
  } else if (!is.null(ego) && nrow(ego) > 0) {
    plot_obj <- dotplot(ego, showCategory = 20) + ggtitle(paste("Cluster", cluster, "- GO BP (unsimplified)"))
  } else {
    cat("No enriched GO BP terms for cluster", cluster, "- skipping plot.\n")
    next
  }
  
  # Save plot
  ggsave(paste0("PEA/Results/", cluster, "_GO_BP_dotplot.png"), plot = plot_obj, width = 10, height = 8)
  
  # Clean up variables to prevent contamination in the next loop
  rm(ego, ego_simplified, gene_entrez, plot_obj)
}