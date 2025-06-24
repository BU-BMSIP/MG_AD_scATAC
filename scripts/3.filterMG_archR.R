.libPaths("/projectnb/cepinet/libs/R_4.4.0_libs")
library(ArchR)
library(pheatmap)

addArchRThreads(threads = 16)
addArchRGenome("hg38")

# ==================== Marker File Paths ====================
files.in.marker <- c(
  MG.marker       = "/projectnb/cepinet/users/vhe/Na_Cell_2023_MG/suppTable/microglia.markers.human.txt",
  MG.state.marker = "/projectnb/cepinet/users/vhe/Na_Cell_2023_MG/suppTable/ROSMAP.Microglia.6regions.seurat.harmony.selected.clusterDEGs.txt"
)

# ==================== Sample Info ====================
df <- read.table(
  '/projectnb/cepinet/data/scATAC/Na_Cell_2023_MG/All.ATAC.samp.info.txt',
  header = TRUE,
  sep = '\t'
)
table(df$region)
df <- df[df$region != 'MB', ]

# ==================== Output Paths ====================
dir.ou.brain <- "brain.microglia.filter/"
dir.create(dir.ou.brain, showWarnings = FALSE, recursive = TRUE)
file.ou.brain.RDS  <- paste0(dir.ou.brain, "brain.rds")
file.ou.brain2.RDS <- paste0(dir.ou.brain, "brain2.rds")

# ==================== Load & Filter ArchRProject ====================
brain <- readRDS("brain.microglia/brain2.rds")
meta <- brain@cellColData
table(meta$Clusters)

remove <- c("C1", "C8")
keep <- !(brain$Clusters %in% remove)
brain <- brain[keep, ]

brain@projectMetadata$outputDirectory <- dir.ou.brain

# ==================== Dimensionality Reduction ====================
if (!file.exists(file.ou.brain.RDS)) {
  brain <- addIterativeLSI(
    ArchRProj     = brain,
    useMatrix     = "TileMatrix",
    name          = "IterativeLSI",
    iterations    = 10,
    clusterParams = list(resolution = 1, sampleCells = 10000, n.start = 10),
    varFeatures   = 25000,
    dimsToUse     = 1:30,
    force         = TRUE
  )

  brain <- addHarmony(
    ArchRProj    = brain,
    reducedDims  = "IterativeLSI",
    name         = "Harmony",
    groupBy      = "Sample",
    force        = TRUE
  )

  brain <- addClusters(
    input        = brain,
    reducedDims  = "IterativeLSI",
    method       = "Seurat",
    name         = "Clusters",
    resolution   = 1,
    force        = TRUE
  )

  brain <- addClusters(
    input        = brain,
    reducedDims  = "Harmony",
    method       = "Seurat",
    name         = "Harmony.Clusters",
    resolution   = 0.5
  )

  cM <- confusionMatrix(paste0(brain$Clusters), paste0(brain$Sample))
  cM <- cM / Matrix::rowSums(cM)

  cM.harmony <- confusionMatrix(paste0(brain$Harmony.Clusters), paste0(brain$Sample))
  cM.harmony <- cM.harmony / Matrix::rowSums(cM.harmony)

  pdf(paste0(dir.ou.brain, "Plots/samples_cluster.pheatmap.pdf"), width = 10, height = 8)
  pheatmap::pheatmap(as.matrix(cM), color = paletteContinuous("whiteBlue"), border_color = "black", main = "Clusters")
  pheatmap::pheatmap(as.matrix(cM.harmony), color = paletteContinuous("whiteBlue"), border_color = "black", main = "harmony.Clusters")
  dev.off()

  saveRDS(brain, file = file.ou.brain.RDS)

  # ========== UMAP ==========
  brain <- addUMAP(brain, reducedDims = "IterativeLSI", name = "UMAP", nNeighbors = 30, minDist = 0.5, metric = "cosine", force = TRUE)
  plotPDF(
    plotEmbedding(brain, "cellColData", "Sample", "UMAP"),
    plotEmbedding(brain, "cellColData", "Clusters", "UMAP"),
    plotEmbedding(brain, "cellColData", "region", "UMAP"),
    name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = brain, addDOC = FALSE, width = 6, height = 6
  )

  brain <- addUMAP(brain, reducedDims = "Harmony", name = "UMAPHarmony", nNeighbors = 30, minDist = 0.5, metric = "cosine", force = TRUE)
  plotPDF(
    plotEmbedding(brain, "cellColData", "Sample", "UMAPHarmony"),
    plotEmbedding(brain, "cellColData", "Harmony.Clusters", "UMAPHarmony"),
    plotEmbedding(brain, "cellColData", "region", "UMAPHarmony"),
    name = "Plot-UMAP2Harmony-Sample-Clusters.pdf", ArchRProj = brain, addDOC = FALSE, width = 6, height = 6
  )
  saveRDS(brain, file = file.ou.brain.RDS)
} else {
  brain2 <- readRDS(file.ou.brain.RDS)
}

# ==================== Marker Gene Analysis ====================
nm2embAndClu <- list(
  LSI     = c(umap = "UMAP", cluster = "Clusters"),
  Harmony = c(umap = "UMAPHarmony", cluster = "Harmony.Clusters")
)
nm.slct <- "Harmony"

markersGS <- getMarkerFeatures(
  ArchRProj  = brain2,
  useMatrix  = "GeneScoreMatrix",
  groupBy    = nm2embAndClu[[nm.slct]]["cluster"],
  bias       = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
saveRDS(markerList, file = "brain.microglia.filter/markerList.rds")

markerGenes.mic <- read.table(files.in.marker["MG.marker"], sep = '\t')$V1
rnadeg <- read.table(files.in.marker["MG.state.marker"], header = TRUE, sep = '\t')
markerGenes <- unique(as.character(rnadeg$gene))

heatmapGS <- markerHeatmap(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25", labelMarkers = markerGenes, transpose = TRUE)
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 12, height = 6, ArchRProj = brain2, addDOC = FALSE)

heatmapGS.mic <- markerHeatmap(markersGS, cutOff = "FDR <= 0.1 & Log2FC >= 0.5", labelMarkers = markerGenes.mic, transpose = TRUE)
plotPDF(heatmapGS.mic, name = "GeneScores-micMarker-Heatmap", width = 12, height = 6, ArchRProj = brain2, addDOC = FALSE)

heatmapGS.both <- markerHeatmap(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1", labelMarkers = c(markerGenes, markerGenes.mic), transpose = TRUE)
plotPDF(heatmapGS.both, name = "GeneScores-bothMarker-Heatmap", width = 20, height = 6, ArchRProj = brain2, addDOC = FALSE)

atacgenes <- colnames(heatmapGS@matrix)
ovgenes <- intersect(atacgenes, markerGenes)

p <- plotEmbedding(brain2, colorBy = "GeneScoreMatrix", name = ovgenes, embedding = nm2embAndClu[[nm.slct]]["umap"], quantCut = c(0.01, 0.95), imputeWeights = NULL)
plotPDF(p, name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf", ArchRProj = brain2, addDOC = FALSE, width = 5, height = 5)

brain2 <- addImputeWeights(brain2)
saveRDS(brain2, file = file.ou.brain2.RDS)

p <- plotEmbedding(brain2, colorBy = "GeneScoreMatrix", name = ovgenes, embedding = nm2embAndClu[[nm.slct]]["umap"], imputeWeights = getImputeWeights(brain2))
plotPDF(p, name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", ArchRProj = brain2, addDOC = FALSE, width = 5, height = 5)

p <- plotBrowserTrack(brain2, groupBy = nm2embAndClu[[nm.slct]]["cluster"], geneSymbol = ovgenes, upstream = 50000, downstream = 50000)
plotPDF(p, name = "Plot-Tracks-Marker-Genes.pdf", ArchRProj = brain2, addDOC = FALSE, width = 5, height = 5)

markerGenes.mic <- markerGenes.mic[-2]
p <- plotEmbedding(brain2, colorBy = "GeneScoreMatrix", name = markerGenes.mic, embedding = nm2embAndClu[[nm.slct]]["umap"], quantCut = c(0.01, 0.95), imputeWeights = NULL)
plotPDF(p, name = "Plot-UMAP-Mic_Marker-Genes-WO-Imputation.pdf", ArchRProj = brain2, addDOC = FALSE, width = 5, height = 5)

p <- plotEmbedding(brain2, colorBy = "GeneScoreMatrix", name = markerGenes.mic, embedding = nm2embAndClu[[nm.slct]]["umap"], imputeWeights = getImputeWeights(brain2))
plotPDF(p, name = "Plot-UMAP-Mic_Marker-Genes-W-Imputation.pdf", ArchRProj = brain2, addDOC = FALSE, width = 5, height = 5)

p <- plotBrowserTrack(brain2, groupBy = nm2embAndClu[[nm.slct]]["cluster"], geneSymbol = markerGenes.mic, upstream = 50000, downstream = 50000)
plotPDF(p, name = "Plot-Tracks-Mic_Marker-Genes.pdf", ArchRProj = brain2, addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = brain, outputDirectory = "Save-brain", load = FALSE)
saveRDS(brain2, file = file.ou.brain2.RDS)