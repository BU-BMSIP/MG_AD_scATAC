# -------------------------------
# Step 5: Plot Dot Plot (DotMap)
# -------------------------------

library(ggplot2)
library(ggrepel)

# Prepare data frame for ggplot
plot_df <- expand.grid(atac = rownames(log_or_mat), rna = colnames(log_or_mat))
plot_df$logOR   <- as.vector(log_or_mat)
plot_df$pvalue  <- as.vector(pval_mat)
plot_df$negLogP <- -log10(plot_df$pvalue)

# Replace infinite or NA values
plot_df$logOR[is.infinite(plot_df$logOR)] <- NA
plot_df$negLogP[is.infinite(plot_df$negLogP)] <- NA

# Cap values for display
plot_df$logOR <- pmin(pmax(plot_df$logOR, -2), 4)
plot_df$negLogP[plot_df$negLogP > 10] <- 10  # cap -log10(p) at 10

# Plot dot plot
p <- ggplot(plot_df, aes(x = rna, y = atac)) +
  geom_point(aes(size = negLogP, color = logOR), alpha = 0.8) +
  scale_color_gradient2(
    low = "blue", mid = "white", high = "red", midpoint = 0,
    limits = c(-2, 4), name = "log2(OR)"
  ) +
  scale_size_continuous(name = "-log10(p)", range = c(1, 8)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey90")
  ) +
  labs(
    title = "Marker Gene Overlap: scATAC vs scRNA",
    x = "scRNA Clusters",
    y = "scATAC Clusters"
  )

# Save plot
ggsave("marker_overlap_dotplot.png", p, width = 10, height = 8)