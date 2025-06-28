library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(dplyr)

# Load your marker gene list
markerList <- readRDS("brain.microglia.filter/markerList.rds")
filtered_markerList <- Filter(function(x) !is.null(x) && nrow(x) > 0, markerList)

# Combine all cluster markers into one data.frame
combined_df <- do.call(rbind, lapply(names(filtered_markerList), function(clu) {
  df <- filtered_markerList[[clu]]
  df$cluster <- clu
  return(df)
})) %>% as.data.frame()

# Select top 300 genes ranked by Log2FC
top300 <- combined_df %>%
  arrange(desc(Log2FC)) %>%
  distinct(name, .keep_all = TRUE) %>%
  head(300)

# Convert gene symbols to Entrez IDs
gene.df <- bitr(top300$name, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

# Run GO enrichment on Biological Process (BP)
ego <- enrichGO(gene = gene.df$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pvalueCutoff = 0.1,
                qvalueCutoff = 0.2,
                readable = TRUE)

# Convert to data.frame for plotting
ego_df <- as.data.frame(ego)

# Select top 20 terms based on adjusted p-value
top_terms <- ego_df %>%
  arrange(p.adjust) %>%
  head(20)

# Plot bubble chart
p <- ggplot(top_terms, aes(x = GeneRatio, y = reorder(Description, GeneRatio))) +
  geom_point(aes(size = Count, color = p.adjust)) +
  scale_color_gradient(low = "red", high = "blue", trans = "log10") +
  labs(title = "GO Enrichment (Top 300 Marker Genes)",
       x = "Gene Ratio",
       y = "GO Term",
       color = "Adjusted p-value",
       size = "Gene Count") +
  theme_minimal(base_size = 12)

# Save the plot
ggsave("top300_marker_gene_GO_bubble.pdf", plot = p, width = 8, height = 6)