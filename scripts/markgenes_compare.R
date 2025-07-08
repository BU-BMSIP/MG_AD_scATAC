# -------------------------------
# Step 5: Dot Plot of Cluster Overlap (scATAC vs scRNA)
# -------------------------------

library(ggplot2)
library(ggrepel)

# Prepare long-format data for ggplot
plot_df <- expand.grid(
  atac_cluster = rownames(log_or_mat),
  rna_cluster = colnames(log_or_mat),
  stringsAsFactors = FALSE
)

plot_df$logOR   <- as.vector(log_or_mat)
plot_df$pvalue  <- as.vector(pval_mat)

# Compute -log10(p) and handle infinite/NA
plot_df$negLogP <- -log10(plot_df$pvalue)
plot_df$logOR[!is.finite(plot_df$logOR)] <- NA
plot_df$negLogP[!is.finite(plot_df$negLogP)] <- NA

# Cap values to improve visualization
logOR_cap <- c(-2, 4)
negLogP_cap <- 10

plot_df$logOR   <- pmax(pmin(plot_df$logOR, logOR_cap[2]), logOR_cap[1])
plot_df$negLogP <- pmin(plot_df$negLogP, negLogP_cap)

# Plot
dotplot <- ggplot(plot_df, aes(x = rna_cluster, y = atac_cluster)) +
  geom_point(aes(size = negLogP, color = logOR), alpha = 0.8) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = logOR_cap, oob = scales::squish, name = expression(log[2]*"(OR)")
  ) +
  scale_size_continuous(
    name = expression(-log[10]*"(p-value)"),
    range = c(1, 8),
    limits = c(0, negLogP_cap)
  ) +
  coord_fixed(ratio = 1) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    legend.key.height = unit(0.8, "cm")
  ) +
  labs(
    title = "Overlap of Marker Genes Between scATAC and scRNA Clusters",
    x = "scRNA Clusters",
    y = "scATAC Clusters"
  )

# Save plot
ggsave("marker_overlap_dotplot.png", dotplot, width = 10, height = 8, dpi = 300)