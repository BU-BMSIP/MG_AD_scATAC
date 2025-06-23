
.libPaths("/projectnb/cepinet/libs/R_4.4.0_libs")
library(ArchR)
library(pheatmap)

addArchRThreads(threads = 16)
addArchRGenome("hg38") 


file.in.RNA.rds="/projectnb/cepinet/data/scRNA/cell-2023-Sun/ROSMAP.Microglia.6regions.seurat.harmony.selected.deidentified.rds"
file.in.brain2.rds="brain.microglia.filter/brain2.rds"
files.in.marker=c(MG.marker="/projectnb/cepinet/users/vhe/Na_Cell_2023_MG/suppTable/microglia.markers.human.txt",
                  MG.state.marker="/projectnb/cepinet/users/vhe/Na_Cell_2023_MG/suppTable/ROSMAP.Microglia.6regions.seurat.harmony.selected.clusterDEGs.txt")

file.ou.brain3.rds="brain.microglia.filter/brain3.Integ.RNA.rds"
file.ou.brain4.rds="brain.microglia.filter/brain4.integRNA.filt.rds"



##################################################################
######### cross-platform linkage of scATAC-seq cells with scRNA-seq cells
################################################################## 


brain.rna=readRDS(file.in.RNA.rds)
table(brain.rna@meta.data$seurat_clusters)
brain.rna$groupRNA=paste('MG',brain.rna@meta.data$seurat_clusters,sep='')


brain2=readRDS(file.in.brain2.rds)


if(!file.exists(file.ou.brain3.rds))
{

    brain3 <- addGeneIntegrationMatrix(
    ArchRProj = brain2, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = brain.rna,
    addToArrow = F,
    groupRNA = "groupRNA",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
  )
  saveRDS(brain3, file.ou.brain3.rds)


  cM <- as.matrix(confusionMatrix(brain3$Clusters, brain3$predictedGroup_Un))
  preClust <- colnames(cM)[apply(cM, 1 , which.max)]
  cbind(preClust, rownames(cM))
  table(brain3$Clusters)

  unique(brain3$predictedGroup_Un)

  p1 <- plotEmbedding(
    brain3, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un", 
  )
  p2 <- plotEmbedding(
    brain3, 
    colorBy = "cellColData", 
    name = "Clusters", 
  )


  plotPDF(p1, p2, name = "Plot-UMAP-RNA-Integration.pdf", ArchRProj = brain3, addDOC = FALSE, width = 5, height = 5)

  getAvailableMatrices(brain3)
  brain3 <- addImputeWeights(brain3)

  ##### 

  meta=brain3@cellColData
  write.table(meta,file='brain.microglia.filter/brain3.meta.txt',sep='\t',quote=F)
  saveRDS(brain3,file=file.ou.brain3.rds)


}else 
{
  brain3=readRDS(file.ou.brain3.rds)
  meta=brain3@cellColData
}

table(meta$Clusters)
remove<-c("C2","C3", "C14") ## remove clusters with very few cells
keep<-!(brain3$Clusters %in% remove)
sum(keep=="TRUE") ## the number of cell to be kept
brain4=brain3[keep,] ## remove small clusters

remapClust <- c("C1" = "C1","C4" = "C2","C5" = "C3","C6" = "C4","C7" = "C5","C8" = "C6","C9" = "C7","C10" = "C8","C11" = "C9","C12" = "C10","C13" = "C11")
labelNew <- mapLabels(names(remapClust), oldLabels = names(remapClust), newLabels = remapClust)
labelNew

brain4$Clusters <- mapLabels(brain4$Clusters, newLabels = labelNew, oldLabels = names(remapClust))


p1 <- plotEmbedding(ArchRProj = brain4, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = brain4, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = brain4, colorBy = "cellColData", name = "region", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = brain4, colorBy = "cellColData", name = "projid", embedding = "UMAP")

plotPDF(p1,p2,p3,p4, name = "Plot-UMAP-Sample-Clusters.filtered.pdf", ArchRProj = brain4, addDOC = FALSE, width = 6, height = 6)


#### plot marker genes 

##### identify marker genes ############################################
##############################################################
### 
markersGS <- getMarkerFeatures(
  ArchRProj = brain4, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)
saveRDS(brain4,file=file.ou.brain4.rds)
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

##################################################################
######### markerGenes.mic UMAP ##
##################################################################
markerGenes.mic
markerGenes.mic=markerGenes.mic[-2]
mymat=getMatrixFromProject(
  ArchRProj = brain4,
  useMatrix = "GeneScoreMatrix")
brain4 <- addImputeWeights(brain4)
matGS <- imputeMatrix(assay(mymat), getImputeWeights(brain4))
dgc_imput_mat <- as(matGS, "dgCMatrix")
str(dgc_imput_mat)
rownames(dgc_imput_mat)=mymat@elementMetadata$name
saveRDS(dgc_imput_mat,file='brain.microglia.filter/brain4.imput_genescoremat.rds')

umaps=getEmbedding(ArchRProj = brain4, embedding = "UMAP", returnDF = TRUE)

saveRDS(brain4,file=file.ou.brain4.rds)
saveRDS(umaps,file="brain.microglia.filter/brain4.integRNA.filt.umap.rds")


###################

p <- plotEmbedding(
  ArchRProj = brain4, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes.mic, 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)

plotPDF(plotList = p, 
        name = "Plot-UMAP-Mic_Marker-Genes-WO-Imputation.pdf", 
        ArchRProj = brain4, 
        addDOC = FALSE, width = 5, height = 5)
######### Marker Genes Imputation with MAGIC
p <- plotEmbedding(
  ArchRProj = brain4, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes.mic, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(brain4)
)

plotPDF(plotList = p, 
        name = "Plot-UMAP-Mic_Marker-Genes-W-Imputation.pdf", 
        ArchRProj = brain4, 
        addDOC = FALSE, width = 5, height = 5)

## Track Plotting with ArchRBrowser
p <- plotBrowserTrack(
  ArchRProj = brain4, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes.mic, 
  upstream = 50000,
  downstream = 50000
)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Mic_Marker-Genes.pdf", 
        ArchRProj = brain4, 
        addDOC = FALSE, width = 5, height = 5)