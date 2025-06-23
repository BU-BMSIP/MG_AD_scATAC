# date: 2025-06-17
.libPaths("/projectnb/cepinet/libs/R_4.4.0_libs")
library(ArchR)
########
df=read.table('/projectnb/cepinet/data/scATAC/Na_Cell_2023_MG/All.ATAC.samp.info.txt',header = T,sep = '\t')
table(df$region)
df=df[df$region!='MB',]


inputFiles=list.files(path = '/projectnb/cepinet/data/scATAC/Na_Cell_2023_MG/fragments',pattern = 'tsv.gz',full.names = T)
sNames=c()
for (f in inputFiles){nm=strsplit(rev(strsplit(f,'/')[[1]])[1],'[.]')[[1]][1];sNames=c(sNames,nm)}
names(inputFiles)=sNames

ov=intersect(sNames,df$SampID)
df=df[df$SampID %in% ov,]
sNames=ov
inputFiles=inputFiles[ov]

dir.ou.arrow="./arrowFiles/"
dir.create(dir.ou.arrow, showWarnings=F, recursive=T)

#####################################################################
#building arrow files
# #
# TileMatrix: binned accessibility (e.g., 500bp bins)
# GeneScoreMatrix: estimated gene activity scores
# Cell-level QC metrics: fragment counts, TSS enrichment, etc.
#####################################################################
## setup 
addArchRThreads(threads = 16)
addArchRGenome("hg38") 
####
ArrowFiles <- createArrowFiles(inputFiles = inputFiles,
                               sampleNames = sNames,
                               outputNames= paste0(dir.ou.arrow, sNames),
                               filterTSS = 6, 
                               filterFrags = 1000, 
                               addTileMat = TRUE,
                               addGeneScoreMat = TRUE)