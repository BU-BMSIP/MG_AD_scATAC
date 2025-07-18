#============================#
#         Set Up            #
#============================#
.libPaths("/projectnb/cepinet/libs/R_4.4.0_libs")
library(ArchR)
library(pheatmap)

addArchRThreads(threads = 16)
addArchRGenome("hg38")

#----------------------------#
#       File Paths          #
#----------------------------#
file.in.RNA.rds      <- "/projectnb/cepinet/data/scRNA/cell-2023-Sun/ROSMAP.Microglia.6regions.seurat.harmony.selected.deidentified.rds"
file.in.brain2.rds   <- "brain.microglia.filter/brain2.rds"
file.ou.brain3.rds   <- "brain.microglia.filter/brain3.Integ.RNA.rds"
file.ou.brain4.rds   <- "brain.microglia.filter/brain4.integRNA.filt.rds"

files.in.marker <- c(
  MG.marker        = "/projectnb/cepinet/users/vhe/Na_Cell_2023_MG/suppTable/microglia.markers.human.txt",
  MG.state.marker  = "/projectnb/cepinet/users/vhe/Na_Cell_2023_MG/suppTable/ROSMAP.Microglia.6regions.seurat.harmony.selected.clusterDEGs.txt"
)

#============================#
#      Load Data            #
#============================#
brain.rna <- readRDS(file.in.RNA.rds)
brain.rna$groupRNA <- paste("MG", brain.rna$seurat_clusters, sep = "")

brain2 <- readRDS(file.in.brain2.rds)

#============================#
#     Integration Step      #
#============================#
if (!file.exists(file.ou.brain3.rds)) {
  brain3 <- addGeneIntegrationMatrix(
    ArchRProj     = brain2,
    useMatrix     = "GeneScoreMatrix",
    matrixName    = "GeneIntegrationMatrix",
    reducedDims   = "IterativeLSI",
    seRNA         = brain.rna,
    addToArrow    = FALSE,
    groupRNA      = "groupRNA",
    nameCell      = "predictedCell_Un",
    nameGroup     = "predictedGroup_Un",
    nameScore     = "predictedScore_Un"
  )

  # Save intermediate results
  saveRDS(brain3, file.ou.brain3.rds)

  # Plot UMAPs
  p1 <- plotEmbedding(brain3, colorBy = "cellColData", name = "predictedGroup_Un")
  p2 <- plotEmbedding(brain3, colorBy = "cellColData", name = "Clusters")
  plotPDF(p1, p2, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = brain3, addDOC = FALSE)

  # Add imputation weights
  brain3 <- addImputeWeights(brain3)

  # Save metadata and project
  write.table(brain3@cellColData, file = "brain.microglia.filter/brain3.meta.txt", sep = '\t', quote = FALSE)
  saveRDS(brain3, file = file.ou.brain3.rds)
} else {
  brain3 <- readRDS(file.ou.brain3.rds)
}

#============================#
#     Filter & Relabel      #
#============================#
remove_clusters <- c("C2", "C3", "C14")
keep_cells <- !(brain3$Clusters %in% remove_clusters)
brain4 <- brain3[keep_cells,]

cluster_map <- c(
  "C1" = "C1",  "C4" = "C2",  "C5" = "C3",  "C6" = "C4",
  "C7" = "C5",  "C8" = "C6",  "C9" = "C7",  "C10" = "C8",
  "C11" = "C9", "C12" = "C10","C13" = "C11"
)
brain4$Clusters <- mapLabels(brain4$Clusters, newLabels = cluster_map, oldLabels = names(cluster_map))

#============================#
#      Visualization        #
#============================#
p1 <- plotEmbedding(brain4, colorBy = "cellColData", name = "Sample")
p2 <- plotEmbedding(brain4, colorBy = "cellColData", name = "Clusters")
p3 <- plotEmbedding(brain4, colorBy = "cellColData", name = "region")
p4 <- plotEmbedding(brain4, colorBy = "cellColData", name = "projid")
plotPDF(p1, p2, p3, p4, name = "Plot-UMAP-Sample-Clusters.filtered.pdf", ArchRProj = brain4)

#============================#
#   Marker Gene Detection   #
#============================#
markersGS <- getMarkerFeatures(
  ArchRProj  = brain4,
  useMatrix  = "GeneScoreMatrix",
  groupBy    = "Clusters",
  bias       = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
saveRDS(markerList, "brain.microglia.filter/markerList.rds")
saveRDS(brain4, file.ou.brain4.rds)

#============================#
#     Marker Heatmaps       #
#============================#
markerGenes.mic <- as.character(read.table(files.in.marker["MG.marker"], sep = '\t')$V1)
markerGenes      <- unique(as.character(read.table(files.in.marker["MG.state.marker"], header = TRUE, sep = '\t')$gene))

# Heatmap 1: AD state marker genes
heatmapGS <- markerHeatmap(
  markersGS, "FDR <= 0.01 & Log2FC >= 1.25",
  labelMarkers = markerGenes,
  transpose = TRUE
)
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 12, height = 6, ArchRProj = brain4)

# Heatmap 2: General microglia marker genes
heatmapGS.mic <- markerHeatmap(
  markersGS, "FDR <= 0.1 & Log2FC >= 0.5",
  labelMarkers = markerGenes.mic,
  transpose = TRUE
)
plotPDF(heatmapGS.mic, name = "GeneScores-micMarker-Heatmap", width = 12, height = 6, ArchRProj = brain4)

# Heatmap 3: Combined
heatmapGS.both <- markerHeatmap(
  markersGS, "FDR <= 0.05 & Log2FC >= 1",
  labelMarkers = c(markerGenes, markerGenes.mic),
  transpose = TRUE
)
plotPDF(heatmapGS.both, name = "GeneScores-bothMarker-Heatmap", width = 20, height = 6, ArchRProj = brain4)

#============================#
#     Imputed Gene Scores   #
#============================#
markerGenes.mic <- markerGenes.mic[-2]  # remove one gene if needed

mymat <- getMatrixFromProject(brain4, "GeneScoreMatrix")
brain4 <- addImputeWeights(brain4)
matGS <- imputeMatrix(assay(mymat), getImputeWeights(brain4))

dgc_imput_mat <- as(matGS, "dgCMatrix")
rownames(dgc_imput_mat) <- mymat@elementMetadata$name
saveRDS(dgc_imput_mat, "brain.microglia.filter/brain4.imput_genescoremat.rds")

#============================#
#        UMAPs & Tracks      #
#============================#
umaps <- getEmbedding(brain4, "UMAP", returnDF = TRUE)
saveRDS(umaps, "brain.microglia.filter/brain4.integRNA.filt.umap.rds")

# Plot marker gene UMAPs (raw)
p <- plotEmbedding(brain4, colorBy = "GeneScoreMatrix", name = markerGenes.mic, embedding = "UMAP", quantCut = c(0.01, 0.95))
plotPDF(p, name = "Plot-UMAP-Mic_Marker-Genes-WO-Imputation.pdf", ArchRProj = brain4)

# Plot marker gene UMAPs (imputed)
p <- plotEmbedding(brain4, colorBy = "GeneScoreMatrix", name = markerGenes.mic, embedding = "UMAP", imputeWeights = getImputeWeights(brain4))
plotPDF(p, name = "Plot-UMAP-Mic_Marker-Genes-W-Imputation.pdf", ArchRProj = brain4)

# Genome tracks
p <- plotBrowserTrack(brain4, groupBy = "Clusters", geneSymbol = markerGenes.mic, upstream = 50000, downstream = 50000)
plotPDF(p, name = "Plot-Tracks-Mic_Marker-Genes.pdf", ArchRProj = brain4)