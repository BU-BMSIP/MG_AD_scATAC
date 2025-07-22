# Load required libraries
library(ArchR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)

# Define a function to detach all non-base R packages
detach_all_libs <- function() {
  loaded_pkgs <- setdiff(
    grep("^package:", search(), value = TRUE),
    paste0("package:", c("base", "stats", "graphics", "grDevices", "utils", "datasets", "methods", "tools"))
  )
  for (pkg in loaded_pkgs) {
    try(detach(pkg, character.only = TRUE, unload = TRUE), silent = TRUE)
  }
  message("All non-base packages detached.")
}

# Load marker gene score and marker list objects
markersGS <- readRDS("brain.microglia.filter/markersGS_7.2.rds")
markerList <- readRDS("brain.microglia.filter/markerList_7.10.rds")

# Convert all gene names to uppercase
markerList <- lapply(markerList, function(df) {
  df$name <- toupper(df$name)
  df
})

# Create a directory to save results
dir.create("PEA_7.11/Results", recursive = TRUE, showWarnings = FALSE)

# Use all genes from the markerGS object as the background gene set
background_genes <- rowData(markersGS)$name
background_genes <- toupper(background_genes)

# Convert SYMBOL to ENTREZID
background_entrez <- bitr(background_genes,
                          fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)

# Optional: report the number of unmapped genes
mapped_genes <- background_entrez$SYMBOL
unmapped_genes <- setdiff(background_genes, mapped_genes)
cat("Unmapped background genes:", length(unmapped_genes), "\n")

# Perform GO enrichment analysis for each cluster
for (cluster in names(markerList)) {
  
  cat("\nProcessing cluster:", cluster, "\n")
  
  # Get top 100 marker gene names
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
  
  # Perform GO Biological Process enrichment analysis
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
  
  # Plot results; if simplification yields nothing, use original
  if (!is.null(ego_simplified) && nrow(ego_simplified) > 0) {
    plot_obj <- dotplot(ego_simplified, showCategory = 20) + 
      ggtitle(paste("Cluster", cluster, "- GO BP"))
  } else if (!is.null(ego) && nrow(ego) > 0) {
    plot_obj <- dotplot(ego, showCategory = 20) + 
      ggtitle(paste("Cluster", cluster, "- GO BP (unsimplified)"))
  } else {
    cat("No enriched GO BP terms for cluster", cluster, "- skipping plot.\n")
    next
  }
  
  # Save the plot to file
  ggsave(paste0("PEA_7.11/Results/", cluster, "_GO_BP_dotplot.png"), 
         plot = plot_obj, width = 10, height = 8)
  
  # Clean up variables to avoid contamination
  rm(ego, ego_simplified, gene_entrez, plot_obj)
}