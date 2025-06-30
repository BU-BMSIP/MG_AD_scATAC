# ===============================================
# GO Enrichment Analysis for Top 300 Marker Genes
# ===============================================

# Load required packages
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

# 1. Load and filter marker gene list
markerList <- readRDS("brain.microglia.filter/markerList.rds")
filtered_markerList <- Filter(function(x) !is.null(x) && nrow(x) > 0, markerList)

# 2. Combine all clusters into one data.frame
combined_df <- bind_rows(lapply(names(filtered_markerList), function(clu) {
  df <- filtered_markerList[[clu]]
  df$cluster <- clu
  df
})) %>% as.data.frame()

# 3. Select top 300 marker genes by Log2FC
top300 <- combined_df %>%
  arrange(desc(Log2FC)) %>%
  distinct(name, .keep_all = TRUE) %>%
  slice_head(n = 300)

# 4. Convert gene symbols to Entrez IDs
gene.df <- tryCatch({
  bitr(
    top300$name,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
}, error = function(e) {
  stop("Gene ID conversion failed: ", e$message)
})

# Filter to keep only successfully converted genes
valid_genes <- top300 %>% filter(name %in% gene.df$SYMBOL)

# 5. GO enrichment (Biological Process)
ego <- enrichGO(
  gene          = gene.df$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pvalueCutoff  = 0.1,
  qvalueCutoff  = 0.2,
  readable      = TRUE
)

# 6. Process and filter GO result
ego_df <- as.data.frame(ego)

top_terms <- ego_df %>%
  arrange(p.adjust) %>%
  slice_head(n = 20)

# 7. Bubble chart plot
p <- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count, color = p.adjust), alpha = 0.8) +
  scale_color_gradient(low = "red", high = "blue", trans = "log10") +
  labs(
    title  = "GO Enrichment (Top 300 Marker Genes)",
    x      = "Gene Ratio",
    y      = "GO Term",
    color  = "Adjusted p-value",
    size   = "Gene Count"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.y = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

# 8. Save plot
ggsave("top300_marker_gene_GO_bubble.pdf", plot = p, width = 8, height = 6)