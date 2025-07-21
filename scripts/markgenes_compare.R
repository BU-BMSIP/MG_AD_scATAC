# =================== Step 0: Setup Environment ===================
rm(list = ls(all.names = TRUE))  # Clear all variables
gc()  # Free up memory
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE)

library(dplyr)
library(pheatmap)
library(RColorBrewer)

# =================== Step 1: Load Marker Gene Lists ===================
atac_markers_raw <- readRDS("/projectnb/cepinet/users/vhe/Na_Cell_2023_MG/brain.microglia.filter/markerList_7.10.rds")
rna_markers      <- readRDS("/projectnb/cepinet/users/vhe/Na_Cell_2023_MG/6.20_meeting_materials/rna_markers_list.rds")

# Extract gene names from ATAC cluster markers (column `name`)
atac_markers <- lapply(atac_markers_raw, function(df) {
  if ("name" %in% colnames(df)) {
    as.character(df$name)
  } else {
    character(0)
  }
})

# =================== Step 2: Build Background Gene Set ===================
# Background is the union of all marker genes across ATAC and RNA clusters
all_genes <- unique(c(unlist(atac_markers), unlist(rna_markers)))

clusters_atac <- names(atac_markers)
clusters_rna  <- names(rna_markers)

# Initialize matrices to store Fisher's test results and overlap counts
odds_matrix <- matrix(NA, nrow = length(clusters_atac), ncol = length(clusters_rna),
                      dimnames = list(clusters_atac, clusters_rna))
pval_matrix <- matrix(NA, nrow = length(clusters_atac), ncol = length(clusters_rna),
                      dimnames = list(clusters_atac, clusters_rna))
overlap_count_matrix <- matrix(0, nrow = length(clusters_atac), ncol = length(clusters_rna),
                               dimnames = list(clusters_atac, clusters_rna))

# =================== Step 3: Fisher's Exact Test ===================
for (i in clusters_atac) {
  genes_i <- atac_markers[[i]]
  for (j in clusters_rna) {
    genes_j <- rna_markers[[j]]
    
    a <- length(intersect(genes_i, genes_j))                      # shared genes
    b <- length(setdiff(genes_i, genes_j))                        # ATAC-only
    c <- length(setdiff(genes_j, genes_i))                        # RNA-only
    d <- length(setdiff(all_genes, union(genes_i, genes_j)))      # neither

    # One-sided Fisher’s Exact Test (alternative = "greater")
    fisher_res <- fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
    
    odds_matrix[i, j] <- fisher_res$estimate
    pval_matrix[i, j] <- fisher_res$p.value
    overlap_count_matrix[i, j] <- a
  }
}

# =================== Step 4: Add Significance Stars Based on Raw P-values ===================
get_significance_stars <- function(pval) {
  if (is.na(pval)) return("")
  else if (pval < 0.001) return("***")
  else if (pval < 0.01) return("**")
  else if (pval < 0.05) return("*")
  else return("")
}

# Apply star annotation to each p-value
stars_matrix <- matrix("", nrow = nrow(pval_matrix), ncol = ncol(pval_matrix),
                       dimnames = dimnames(pval_matrix))

for (i in rownames(pval_matrix)) {
  for (j in colnames(pval_matrix)) {
    stars_matrix[i, j] <- get_significance_stars(pval_matrix[i, j])
  }
}

# =================== Step 5: Heatmap Plot ===================
# Log-transform odds ratio and handle infinite/undefined values
log_odds <- log2(odds_matrix)
log_odds[!is.finite(log_odds)] <- NA

# Define heatmap color scale
breaks <- seq(floor(min(log_odds, na.rm = TRUE)), ceiling(max(log_odds, na.rm = TRUE)), length.out = 100)
colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdYlBu")))(length(breaks) - 1)

# Combine overlap counts and significance stars for display
display_labels <- matrix("",
                         nrow = nrow(log_odds),
                         ncol = ncol(log_odds),
                         dimnames = dimnames(log_odds))

for (i in rownames(display_labels)) {
  for (j in colnames(display_labels)) {
    count <- overlap_count_matrix[i, j]
    star  <- stars_matrix[i, j]
    if (count == 0 && star == "") {
      display_labels[i, j] <- ""
    } else if (count == 0 && star != "") {
      display_labels[i, j] <- star
    } else {
      display_labels[i, j] <- paste0(count, star)
    }
  }
}

# =================== Step 6: Save Heatmap ===================
png("Overlap_Heatmap_logOdds_withStars.png", width = 2000, height = 1600, res = 300)
pheatmap(log_odds,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         color = colors,
         breaks = breaks,
         display_numbers = display_labels,
         number_color = "black",
         fontsize_number = 12,
         main = "Heatmap of marker gene overlap (log2 Odds Ratio)",
         fontsize_row = 10,
         fontsize_col = 10,
         legend = TRUE,
         na_col = "lightgray")
dev.off()