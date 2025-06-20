.libPaths("~/cEpiNet/libs/R_4.4.0_libs/")
library(ArchR)
library(pheatmap)

addArchRThreads(threads = 16)
addArchRGenome("hg38") 



########marker files
files.in.marker=c(MG.marker="/projectnb/cepinet/users/vhe/Na_Cell_2023_MG/suppTable/microglia.markers.human.txt",
                  MG.state.marker="/projectnb/cepinet/users/vhe/Na_Cell_2023_MG/suppTable/ROSMAP.Microglia.6regions.seurat.harmony.selected.clusterDEGs.txt")


########
df=read.table('/projectnb/cepinet/data/scATAC/Na_Cell_2023_MG/All.ATAC.samp.info.txt',header = T,sep = '\t')
table(df$region)
df=df[df$region!='MB',]

dir.ou.brain="brain.microglia.filter/"
dir.create(dir.ou.brain, showWarnings=F, recursive=T)
#file.ou.brain1.pdf=paste0(dir.ou.brain1, "/brain1.pdf")
file.ou.brain.RDS=paste0(dir.ou.brain, "brain.rds")
file.ou.brain2.RDS=paste0(dir.ou.brain, "brain2.rds")
#########
brain=readRDS('brain.microglia/brain2.rds')
meta=brain@cellColData
table(meta$Clusters)

remove<-c("C1","C8") ## remove potential doublet clusters based on whether it is far away from major chunk
keep<-!(brain$Clusters %in% remove)
sum(keep=="TRUE") ## the number of cell to be kept
brain=brain[keep,]

brain@projectMetadata$outputDirectory="/projectnb/cepinet/users/vhe/Na_Cell_2023_MG/brain.microglia.filter/"
#dir.create("brain.microglia.filter/Embeddings", recursive = TRUE, showWarnings = FALSE)
#################################################################
##### Dimensionality Reduction with ArchR########################
#################################################################


if(!file.exists(file.ou.brain.RDS))
{
  brain <- addIterativeLSI(
    ArchRProj = brain,
    useMatrix = "TileMatrix", 
    name = "IterativeLSI", 
    iterations = 10, 
    clusterParams = list( #See Seurat::FindClusters
      resolution = c(1), 
      sampleCells = 10000, 
      n.start = 10
    ), 
    varFeatures = 25000, 
    dimsToUse = 1:30,
    force=TRUE
  )

  ## batch correction with Harmony
  brain <- addHarmony(
    ArchRProj = brain,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample",
    force = TRUE
  )

  #################################################################
  #####  ################################################
  #################################################################
  ## 
  print("Clustering")
  brain <- addClusters(
    input = brain,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 1,
    force = TRUE
  )

  #4.22.2025,  # Use batch-corrected space
  brain <- addClusters(
    input = brain,
    reducedDims = "Harmony", 
    method = "Seurat",
    name = "Harmony.Clusters",
    resolution = 0.5
  )

  table(brain$Clusters)
  cM <- confusionMatrix(paste0(brain$Clusters), paste0(brain$Sample))
  cM <- cM / Matrix::rowSums(cM)

  cM.harmony <-  confusionMatrix(paste0(brain$Harmony.Clusters), paste0(brain$Sample))
  cM.harmony <- cM.harmony / Matrix::rowSums(cM.harmony)

  pdf(paste0(dir.ou.brain, "Plots/samples_cluster.pheatmap.pdf"), width = 10, height = 8)
  pheatmap::pheatmap(
    mat = as.matrix(cM), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black",
    main="Clusters"
  )
  pheatmap::pheatmap(
    mat = as.matrix(cM.harmony), 
    color = paletteContinuous("whiteBlue"), 
    border_color = "black",
    main="harmony.Clusters"
  )
  dev.off()

  saveRDS(brain,file='brain.microglia.filter/brain.rds')
  #################################################################
  ##### UMAP ######################################################
  #################################################################
  ### UMAP
  print("UMAP")
  brain <- addUMAP(
    ArchRProj = brain, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE
  )

  p1 <- plotEmbedding(ArchRProj = brain, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
  p2 <- plotEmbedding(ArchRProj = brain, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
  p3 <- plotEmbedding(ArchRProj = brain, colorBy = "cellColData", name = "region", embedding = "UMAP")
  p4 <- plotEmbedding(ArchRProj = brain, colorBy = "cellColData", name = "projid", embedding = "UMAP")


  plotPDF(p1,p2,p3, name = "Plot-UMAP-Sample-Clusters.pdf", ArchRProj = brain, addDOC = FALSE, width = 6, height = 6)
  saveRDS(brain,file='brain.microglia.filter/brain.rds')

  ## Dimensionality Reduction After Harmony

  brain <- addUMAP(
    ArchRProj = brain, 
    reducedDims = "Harmony", 
    name = "UMAPHarmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine",
    force = TRUE
  )

  p1 <- plotEmbedding(ArchRProj = brain, colorBy = "cellColData", name = "Sample", embedding = "UMAPHarmony")
  p2 <- plotEmbedding(ArchRProj = brain, colorBy = "cellColData", name = "Harmony.Clusters", embedding = "UMAPHarmony")
  p3 <- plotEmbedding(ArchRProj = brain, colorBy = "cellColData", name = "region", embedding = "UMAPHarmony")

  plotPDF(p1,p2,p3, name = "Plot-UMAP2Harmony-Sample-Clusters.pdf", ArchRProj = brain, addDOC = FALSE, width = 6, height = 6)
  
  saveRDS(brain, file=file.ou.brain.RDS)
}else brain2=readRDS(file.ou.brain.RDS)

#################################################################
##### identify marker genes ############################################
##############################################################
### 

nm2embAndClu=list()
nm2embAndClu[["LSI"]]=c(umap="UMAP",cluster="Clusters")
nm2embAndClu[["Harmony"]]=c(umap="UMAPHarmony",cluster="Harmony.Clusters")

nm.slct="Harmony"

print("identify marker genes")
markersGS <- getMarkerFeatures(
  ArchRProj = brain2, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = nm2embAndClu[[nm.slct]]["cluster"], #Cluster
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
saveRDS(markerList,file='brain.microglia.filter/markerList.rds')

df=read.table(files.in.marker["MG.marker"],sep = '\t')
markerGenes.mic = as.character(df$V1)
rnadeg=read.table(files.in.marker["MG.state.marker"],header = T,sep = '\t')
markerGenes=unique(as.character(rnadeg$gene))

######
heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 12, height = 6, ArchRProj = brain2, addDOC = FALSE)

##
heatmapGS.mic <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5", 
  labelMarkers = markerGenes.mic,
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS.mic, name = "GeneScores-micMarker-Heatmap", width = 12, height = 6, ArchRProj = brain2, addDOC = FALSE)
############
heatmapGS.both <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1", 
  labelMarkers = c(markerGenes,markerGenes.mic),
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS.both, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapGS.both, name = "GeneScores-bothMarker-Heatmap", width = 20, height = 6, ArchRProj = brain2, addDOC = FALSE)

## overlap between snRNA subtype markers and snATAC subtype markers
atacgenes=colnames(heatmapGS@matrix)
ovgenes=intersect(atacgenes,markerGenes)
rnadeg.sel=rnadeg[rnadeg$gene %in% ovgenes,]

######### ovgenes UMAP ##
print("ovgenes UMAP")
p <- plotEmbedding(
  ArchRProj = brain2, 
  colorBy = "GeneScoreMatrix", 
  name = ovgenes, 
  embedding = nm2embAndClu[[nm.slct]]["umap"], #"UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)

plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes-WO-Imputation.pdf", 
        ArchRProj = brain2, 
        addDOC = FALSE, width = 5, height = 5)

######### Marker Genes Imputation with MAGIC
brain2 <- addImputeWeights(brain2)
saveRDS(brain2,file=file.ou.brain2.RDS)

p <- plotEmbedding(
  ArchRProj = brain2, 
  colorBy = "GeneScoreMatrix", 
  name = ovgenes, 
  embedding = nm2embAndClu[[nm.slct]]["umap"], #"UMAP",
  imputeWeights = getImputeWeights(brain2)
)

plotPDF(plotList = p, 
        name = "Plot-UMAP-Marker-Genes-W-Imputation.pdf", 
        ArchRProj = brain2, 
        addDOC = FALSE, width = 5, height = 5)

## Track Plotting with ArchRBrowser
p <- plotBrowserTrack(
  ArchRProj = brain2, 
  groupBy = nm2embAndClu[[nm.slct]]["cluster"], #"Clusters", 
  geneSymbol = ovgenes, 
  upstream = 50000,
  downstream = 50000
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes.pdf", 
        ArchRProj = brain2, 
        addDOC = FALSE, width = 5, height = 5)

##################################################################
######### markerGenes.mic UMAP ##
##################################################################
print("markerGenes.mic UMAP")

markerGenes.mic
markerGenes.mic=markerGenes.mic[-2]
p <- plotEmbedding(
  ArchRProj = brain2, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes.mic, 
  embedding = nm2embAndClu[[nm.slct]]["umap"], #"UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)

plotPDF(plotList = p, 
        name = "Plot-UMAP-Mic_Marker-Genes-WO-Imputation.pdf", 
        ArchRProj = brain2, 
        addDOC = FALSE, width = 5, height = 5)


######### Marker Genes Imputation with MAGIC
print("Marker Genes Imputation with MAGIC")

p <- plotEmbedding(
  ArchRProj = brain2, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes.mic, 
  embedding = nm2embAndClu[[nm.slct]]["umap"], #"UMAP",
  imputeWeights = getImputeWeights(brain2)
)

plotPDF(plotList = p, 
        name = "Plot-UMAP-Mic_Marker-Genes-W-Imputation.pdf", 
        ArchRProj = brain2, 
        addDOC = FALSE, width = 5, height = 5)

## Track Plotting with ArchRBrowser
p <- plotBrowserTrack(
  ArchRProj = brain2, 
  groupBy = nm2embAndClu[[nm.slct]]["cluster"], #"Clusters", 
  geneSymbol = markerGenes.mic, 
  upstream = 50000,
  downstream = 50000
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Mic_Marker-Genes.pdf", 
        ArchRProj = brain2, 
        addDOC = FALSE, width = 5, height = 5)

saveArchRProject(ArchRProj = brain, outputDirectory = "Save-brain", load = FALSE)
saveRDS(brain2,file=file.ou.brain2.RDS)