# date: "2025-06-20"

.libPaths("/projectnb/cepinet/libs/R_4.4.0_libs")
library(ArchR)

# ==================== Sample Info ====================
df <- read.table(
  '/projectnb/cepinet/data/scATAC/Na_Cell_2023_MG/All.ATAC.samp.info.txt',
  header = TRUE,
  sep = '\t'
)
table(df$region)
df <- df[df$region != 'MB', ]

# ==================== Output Directory ====================
dir.ou.brain1 <- "brain.microglia/"
dir.create(dir.ou.brain1, showWarnings = FALSE, recursive = TRUE)
file.ou.brain1.pdf <- paste0(dir.ou.brain1, "/brain1.pdf")

# ==================== Arrow Files ====================
ArrowFiles <- list.files(
  path = "/projectnb/cepinet/data/scATAC/Na_Cell_2023_MG/arrowFiles/",
  pattern = 'arrow',
  full.names = TRUE
)

addArchRThreads(threads = 16)
addArchRGenome("hg38")

brain1 <- ArchRProject(
  ArrowFiles       = ArrowFiles, 
  outputDirectory  = dir.ou.brain1,
  copyArrows       = FALSE
)
brain1

paste0("Memory Size = ", round(object.size(brain1) / 10^6, 3), " MB")
getAvailableMatrices(brain1)

# ==================== Metadata Matching ====================
head(brain1$cellNames)
head(brain1$Sample)

rownames(df) <- df$SampID
libs <- as.matrix(unname(as.data.frame(strsplit(brain1$cellNames, '#'))[1, ]))[1, ]
mymeta <- df[libs, ]
rownames(mymeta) <- brain1$cellNames

for (i in 2:ncol(mymeta)) {
  brain1 <- addCellColData(
    ArchRProj = brain1,
    data      = mymeta[, i],
    name      = colnames(mymeta)[i],
    cells     = brain1$cellNames,
    force     = TRUE
  )
}

str(brain1@cellColData)

# ==================== Sample QC Plots ====================
p1 <- plotGroups(
  ArchRProj  = brain1, 
  groupBy    = "Sample", 
  colorBy    = "cellColData", 
  name       = "TSSEnrichment",
  plotAs     = "ridges"
)

p2 <- plotGroups(
  ArchRProj   = brain1, 
  groupBy     = "Sample", 
  colorBy     = "cellColData", 
  name        = "TSSEnrichment",
  plotAs      = "violin",
  alpha       = 0.4,
  addBoxPlot  = TRUE
)

p3 <- plotGroups(
  ArchRProj  = brain1, 
  groupBy    = "Sample", 
  colorBy    = "cellColData", 
  name       = "log10(nFrags)",
  plotAs     = "ridges"
)

p4 <- plotGroups(
  ArchRProj   = brain1, 
  groupBy     = "Sample", 
  colorBy     = "cellColData", 
  name        = "log10(nFrags)",
  plotAs      = "violin",
  alpha       = 0.4,
  addBoxPlot  = TRUE
)

plotPDF(
  p1, p2, p3, p4,
  name       = "QC-Sample-Statistics.pdf",
  ArchRProj  = brain1,
  addDOC     = FALSE,
  width      = 30,
  height     = 30
)

# ==================== Dimensionality Reduction ====================
brain2 <- brain1

brain2 <- addIterativeLSI(
  ArchRProj     = brain2,
  useMatrix     = "TileMatrix", 
  name          = "IterativeLSI", 
  iterations    = 10, 
  clusterParams = list(
    resolution   = 0.5, 
    sampleCells  = 10000, 
    n.start      = 10
  ), 
  varFeatures   = 25000, 
  dimsToUse     = 1:30
)

# Batch correction with Harmony
brain2 <- addHarmony(
  ArchRProj    = brain2,
  reducedDims  = "IterativeLSI",
  name         = "Harmony",
  groupBy      = "projid",
  force        = TRUE
)

# ==================== Clustering ====================
brain2 <- addClusters(
  input        = brain2,
  reducedDims  = "IterativeLSI",
  method       = "Seurat",
  name         = "Clusters",
  resolution   = 0.5
)

table(brain2$Clusters)

cM <- confusionMatrix(paste0(brain2$Clusters), paste0(brain2$Sample))

library(pheatmap)
cM <- cM / Matrix::rowSums(cM)

p <- pheatmap::pheatmap(
  mat          = as.matrix(cM), 
  color        = paletteContinuous("whiteBlue"), 
  border_color = "black"
)

plotPDF(
  p,
  name       = "samples_cluster.pheatmap.pdf",
  ArchRProj  = brain2,
  addDOC     = FALSE,
  width      = 10,
  height     = 8
)

# ==================== UMAP ====================
brain2 <- addUMAP(
  ArchRProj    = brain2, 
  reducedDims  = "IterativeLSI", 
  name         = "UMAP", 
  nNeighbors   = 30, 
  minDist      = 0.5, 
  metric       = "cosine"
)

p1 <- plotEmbedding(brain2, colorBy = "cellColData", name = "Sample",  embedding = "UMAP")
p2 <- plotEmbedding(brain2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(brain2, colorBy = "cellColData", name = "region",  embedding = "UMAP")
p4 <- plotEmbedding(brain2, colorBy = "cellColData", name = "projid",  embedding = "UMAP")

plotPDF(
  p1, p2, p3, p4,
  name       = "Plot-UMAP-Sample-Clusters.pdf",
  ArchRProj  = brain2,
  addDOC     = FALSE,
  width      = 6,
  height     = 6
)

# ========== UMAP After Harmony ==========
brain2 <- addUMAP(
  ArchRProj    = brain2, 
  reducedDims  = "Harmony", 
  name         = "UMAPHarmony", 
  nNeighbors   = 30, 
  minDist      = 0.5, 
  metric       = "cosine",
  force        = TRUE
)

p1 <- plotEmbedding(brain2, colorBy = "cellColData", name = "Sample",  embedding = "UMAPHarmony")
p2 <- plotEmbedding(brain2, colorBy = "cellColData", name = "Clusters", embedding = "UMAPHarmony")
p3 <- plotEmbedding(brain2, colorBy = "cellColData", name = "region",  embedding = "UMAPHarmony")

plotPDF(
  p1, p2, p3,
  name       = "Plot-UMAP2Harmony-projid-Clusters.pdf",
  ArchRProj  = brain2,
  addDOC     = FALSE,
  width      = 6,
  height     = 6
)

# Save final ArchR project object
saveRDS(brain2, file = "brain.microglia/brain2.rds")