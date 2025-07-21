library(ArchR)

# 1. Build trajectory using RNA-predicted cell types as ordered states
brain3 <- addTrajectory(
  ArchRProj = brain3,
  name = "MicrogliaTrajectory",
  groupBy = "predictedGroup_Un",                # RNA-based cell type labels
  trajectory = c("Progenitor", "Homeostatic", "DAM"),  # Define cell state order
  embedding = "UMAP"                            # Use UMAP embedding for trajectory
)

# 2. Add motif annotations for TF binding sites (if not done yet)
brain3 <- addMotifAnnotations(
  ArchRProj = brain3,
  motifSet = "cisbp",   # Motif database to use, alternatives: "encode"
  name = "Motif"
)

# 3. Compute chromVAR deviation scores to quantify TF motif activity per cell
brain3 <- addDeviationsMatrix(
  ArchRProj = brain3,
  peakAnnotation = "Motif",    # Use motif annotations added above
  matrixName = "MotifMatrix"   # Name for output deviation matrix
)

# 4. Visualize CTCF motif activity dynamics along the trajectory
plotTrajectory(
  ArchRProj = brain3,
  trajectory = "MicrogliaTrajectory",
  name = "CTCF",               # TF motif to visualize
  useMatrix = "MotifMatrix"    # Use chromVAR deviation scores
)

# 5. Extract and plot CTCF gene expression along the trajectory from integrated RNA
ctcfExpr <- getTrajectory(
  ArchRProj = brain3,
  name = "MicrogliaTrajectory",
  useMatrix = "GeneIntegrationMatrix",  # RNA expression matrix integrated into ArchR
  log2Norm = TRUE                       # Log2 normalization of expression
)

plotTrajectoryHeatmap(ctcfExpr, labelMarkers = "CTCF")

# 6. Add Peak-to-Gene links to connect CREs with target genes (if not done yet)
brain3 <- addPeak2GeneLinks(
  ArchRProj = brain3,
  reducedDims = "IterativeLSI"          # Use LSI-reduced dimensions for linking
)

# 7. Retrieve peak-to-gene links for downstream analysis (e.g., focusing on CTCF motifs)
peak2gene <- getPeak2GeneLinks(brain3)