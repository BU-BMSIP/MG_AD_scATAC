library(ArchR)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)

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


markersGS <- readRDS("brain.microglia.filter/markersGS_7.2.rds")
markerList <- readRDS("brain.microglia.filter/markerList_7.10.rds")

# 所有基因名转为大写
markerList <- lapply(markerList, function(df) {
  df$name <- toupper(df$name)
  df
})

# 创建结果保存目录
dir.create("PEA_7.11/Results", recursive = TRUE, showWarnings = FALSE)

# 使用 MarkerGS 结果中的所有基因作为背景集
background_genes <- rowData(markersGS)$name
background_genes <- toupper(background_genes)

# 转换为 ENTREZID
background_entrez <- bitr(background_genes,
                          fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)

# 可选：输出无法映射的基因数
mapped_genes <- background_entrez$SYMBOL
unmapped_genes <- setdiff(background_genes, mapped_genes)
cat("Unmapped background genes:", length(unmapped_genes), "\n")

# 2️⃣ 遍历每个 cluster 做 GO 富集分析
for (cluster in names(markerList)) {
  
  cat("\nProcessing cluster:", cluster, "\n")
  
  # 获取 top 100 marker gene 名称
  top_genes <- head(markerList[[cluster]]$name, 100)
  
  if (length(top_genes) == 0 || all(is.na(top_genes))) {
    cat("No genes available for cluster", cluster, "- skipping.\n")
    next
  }
  
  # SYMBOL to ENTREZID
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
  
  # 简化冗余 GO term
  ego_simplified <- simplify(ego, cutoff = 0.7, by = "p.adjust", select_fun = min)
  
  # 画图逻辑，简化后没结果则用原始结果
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
  
  # 保存图片
  ggsave(paste0("PEA_7.11/Results/", cluster, "_GO_BP_dotplot.png"), 
         plot = plot_obj, width = 10, height = 8)
  
  # 清理变量，防止污染
  rm(ego, ego_simplified, gene_entrez, plot_obj)
}