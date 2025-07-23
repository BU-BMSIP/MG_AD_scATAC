library(ArchR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)

# Function to detach all non-base R packages (optional cleanup)
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

# Load precomputed marker results
markersGS <- readRDS("brain.microglia.filter/markersGS_7.2.rds")
markerList <- readRDS("brain.microglia.filter/markerList_7.10.rds")

# Convert all gene names to uppercase to ensure consistent ID mapping
markerList <- lapply(markerList, function(df) {
  df$name <- toupper(df$name)
  df
})

# Create output directory for enrichment results
dir.create("PEA_7.11/Results", recursive = TRUE, showWarnings = FALSE)

# Use all marker genes from markersGS as background universe
background_genes <- rowData(markersGS)$name
background_genes <- toupper(background_genes)

# Map background gene symbols to Entrez IDs
background_entrez <- bitr(background_genes,
                          fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)

# Optional: Print number of unmapped background genes
mapped_genes <- background_entrez$SYMBOL
unmapped_genes <- setdiff(background_genes, mapped_genes)
cat("Unmapped background genes:", length(unmapped_genes), "\n")

# Iterate over all clusters and perform GO enrichment analysis
for (cluster in names(markerList)) {
  
  cat("\nProcessing cluster:", cluster, "\n")
  
  # Get top 100 marker gene symbols for the current cluster
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
  
  # Perform GO Biological Process enrichment
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
  
  # Plot the result using simplified or original enrichment
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
  
  # Save the dotplot to file
  ggsave(paste0("PEA_7.11/Results/", cluster, "_GO_BP_dotplot.png"), 
         plot = plot_obj, width = 10, height = 8)
  
  # Clean up temporary variables
  rm(ego, ego_simplified, gene_entrez, plot_obj)
}