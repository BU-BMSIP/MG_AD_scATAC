# peak calling based on MG cluster
# and generate pseudobulk signal for differential peak calling for sample groups

.libPaths("/projectnb/cepinet/libs/R_4.4.0_libs")
library(ArchR)
library(pheatmap)
library(GenomicRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(DESeq2)
library(GenomicRanges)

addArchRThreads(threads = 16)
addArchRGenome("hg38")

# file.in.RNA.rds = "../../scRNA/cell-2023-Sun/ROSMAP.Microglia.6regions.seurat.harmony.selected.deidentified.rds"
file.ou.brain4.rds  <- "brain.microglia.filter/brain4.integRNA.filt.rds"  # pseudo bulk profile based on clusters
file.in.430.meta     <- "/projectnb/cepinet/data/sc_eQTL/AD430_PFC.eQTLs.Oct2023/1.annotation/patients.metadata.withADdiag.txt"
file.in.pkSet.rds    <- "./brain.microglia.filter/brain4.integRNA.filt.PeakSet.byClusters.rds"

# Load filtered ArchR project
brain4 <- readRDS(file.ou.brain4.rds)

# Update metadata: projid should be character
cell.meta     <- getCellColData(brain4)
projid.new    <- as.character(cell.meta$projid)
names(projid.new) <- rownames(cell.meta)

brain4 <- addCellColData(
  ArchRProj = brain4,
  data      = projid.new,
  cells     = names(projid.new),
  name      = "projid",
  force     = TRUE
)

# Load peak set and add peak matrix
pkSet  <- readRDS(file.in.pkSet.rds)
brain4 <- addPeakSet(brain4, peakSet = pkSet, force = TRUE)
brain4 <- addPeakMatrix(brain4)

# Optionally save updated project
# brain4 <- saveArchRProject(
#   ArchRProj        = brain4,
#   outputDirectory  = "/projectnb/cepinet/data/scATAC/Na_Cell_2023_MG/brain.microglia.filter/",
#   load             = TRUE
# )
# pkSet.test <- getPeakSet(brain4)

saveRDS(brain4, file.ou.brain4.rds)

# Generate pseudobulk matrix by sample (projid)
brain4.pseudobulk.mat <- getGroupSE(
  ArchRProj = brain4,
  divideN   = FALSE,
  useMatrix = "PeakMatrix",
  groupBy   = "projid"  # or your grouping
)

saveRDS(
  brain4.pseudobulk.mat,
  file = 'brain.microglia.filter/brain4.integRNA.filt.pseudoBulkBySamp.peakMat.rds'
)