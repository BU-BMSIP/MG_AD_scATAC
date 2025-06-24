
#integration and peak calling
#and generate pseudobulk signal for differential peak calling for sample groups

.libPaths("~/cEpiNet/libs/R_4.4.0_libs/")
library(ArchR)
library(pheatmap)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(DESeq2)

addArchRThreads(threads = 16)
addArchRGenome("hg38") 


#file.in.RNA.rds="../../scRNA/cell-2023-Sun/ROSMAP.Microglia.6regions.seurat.harmony.selected.deidentified.rds"
file.in.brain2.rds="brain.microglia.filter/brain2.rds"
files.in.marker=c(MG.marker="./suppTable/microglia.markers.human.txt",
                  MG.state.marker="./suppTable/ROSMAP.Microglia.6regions.seurat.harmony.selected.clusterDEGs.txt")
file.in.430.meta="/projectnb/cepinet/data/sc_eQTL/AD430_PFC.eQTLs.Oct2023/1.annotation/patients.metadata.withADdiag.txt"



file.ou.brain3.rds="brain.microglia.filter/brain3.Integ.RNA.rds"
file.ou.brain3.addCoverage.rds="brain.microglia.filter/brain3.integRNA.addCoverageBySamp.rds"



if(!file.exists(file.ou.brain3.addCoverage.rds))
{
    brain3=readRDS(file.ou.brain3.rds)
  
    ## Making Pseudo-bulk Replicates
    brain3 <- addGroupCoverages(ArchRProj = brain3, 
                                groupBy = "Sample",
                                useLabels = FALSE,         # DO NOT create further replicates by grouping by sample
                                minReplicates = 2,
                                maxReplicates = 2,         # Only one replicate per sample
                                minCells = 40,             # Keep this to prevent tiny samples
                                maxCells = 1e6,            # Set very high so no cells are excluded
                                maxFragments = 1e9,        # Also very high to avoid fragment downsampling
                                threads=16, 
                                force=F)

    saveRDS(brain3,file.ou.brain3.addCoverage.rds)
}else
{
    brain3=readRDS(file.ou.brain3.addCoverage.rds)
}
#

#########################################################
####### Calling Peaks w/ Macs2
##########################################################
## 
pathToMacs2 <- findMacs2()

brain3 <- addReproduciblePeakSet(ArchRProj = brain3, 
                                groupBy = "Sample", 
                                pathToMacs2 = pathToMacs2)
saveRDS(brain3, file.ou.brain3.addCoverage.rds)
getPeakSet(brain3)

saveRDS(getPeakSet(brain3),file='brain.microglia.filter/brain3.integRNA.addCoverageBySamp.getPeakSet.rds')

## Add Peak Matrix ??????
brain3 <- addPeakMatrix(brain3)
#getAvailableMatrices(brain3)

cell.meta <- getCellColData(brain3)
# Step 2: Modify projid to character
projid.new <- as.character(cell.meta$projid)
names(projid.new) <- rownames(cell.meta)  # must be named by cell IDs
brain3 <- addCellColData(
  ArchRProj = brain3,
  data = projid.new,
  cells = names(projid.new),
  name = "projid",
  force = TRUE
)



brain3.pseudobulk.mat <- getGroupSE(
  ArchRProj = brain3,
  divideN=F,
  useMatrix = "PeakMatrix",
  groupBy = "projid"
)
saveRDS(brain3.pseudobulk.mat,file='brain.microglia.filter/brain3.integRNA.addCoverageBySamp.pseudobulk.peakMat.rds')

#read meta data
meta.df=read.table(file.in.430.meta, sep="\t", header=T, row.names=1)



# dds <- DESeqDataSet(se.pb, design = ~ 1)  # use real design formula here
# dds <- estimateSizeFactors(dds)
# dds <- estimateDispersions(dds)
# dds <- DESeq(dds)

# res <- results(dds)