# TF Dynamics via Trajectory Analysis – Tool Comparison

## Goal
Detect changes in TF activity (e.g., CTCF loss) along microglial state transitions using snATAC + snRNA integration.

---

## 1. Build Cell State Trajectory

| Feature                           | ArchR                             | Monocle3                          | Cicero               |
|----------------------------------|-----------------------------------|-----------------------------------|----------------------|
| Supports trajectory inference    | ✅ `addTrajectory()`              | ✅ `learn_graph()` + `order_cells()` | ❌ Not supported     |
| Pseudotime based on ATAC + RNA   | ✅ RNA integration (Seurat-based) | ✅ RNA-based                      | ❌ No pseudotime     |
| Manual ordering or root cell     | ✅ Required                       | ✅ Optional or interactive        | ❌ Not applicable     |
| Output: ordered cells in time    | ✅ Yes                             | ✅ Yes                             | ❌ No                |

---

## 2. Analyze TF Motif Activity (chromVAR)

| Feature                                   | ArchR                             | Monocle3 | Cicero |
|------------------------------------------|-----------------------------------|----------|--------|
| Motif annotation                         | ✅ `addMotifAnnotations()`        | ❌       | ❌     |
| TF motif activity scoring (chromVAR)     | ✅ `addDeviationsMatrix()`        | ❌       | ❌     |
| Plot motif accessibility over trajectory | ✅ `plotTrajectory()`             | ❌       | ❌     |
| Track CTCF motif signal                  | ✅ Yes                             | ❌       | ❌     |

---

## 3. Monitor TF Gene Expression (e.g., CTCF)

| Feature                            | ArchR                                 | Monocle3                              | Cicero |
|-----------------------------------|---------------------------------------|----------------------------------------|--------|
| RNA integration                    | ✅ `addGeneIntegrationMatrix()`       | ✅ RNA-native or Seurat converted       | ❌     |
| Plot TF expression along pseudotime| ✅ `plotTrajectoryHeatmap()`, `plotTrajectory()` | ✅ `plot_genes_in_pseudotime()` | ❌     |

---

## 4. Link CREs to Genes (e.g., CTCF-bound regions)

| Feature                                  | ArchR                                 | Monocle3 | Cicero                       |
|-----------------------------------------|---------------------------------------|----------|------------------------------|
| CRE–gene linkage                        | ✅ `addPeak2GeneLinks()`              | ❌       | ✅ `run_cicero()`            |
| Supports pseudotime-aware linkage       | ⚠️ Partial (compare CREs by stage)    | ❌       | ❌ (co-accessibility only)   |
| Motif-informed CRE filtering (e.g., CTCF)| ✅ Overlap with CTCF motif annotations | ❌       | ⚠️ Requires manual handling  |

---

## ✅ Recommended Usage Summary

| Task                          | Best Tool     | Notes                                                       |
|-------------------------------|---------------|-------------------------------------------------------------|
| ATAC + RNA trajectory         | **ArchR**     | Uses RNA-predicted groups to define pseudotime              |
| TF motif activity over time   | **ArchR**     | chromVAR integrated                                         |
| RNA expression over time      | ArchR / Monocle3 | Both are useful for checking TF dynamics (e.g., CTCF)    |
| CRE-gene interaction modeling | ArchR / Cicero | ArchR = peak2gene links; Cicero = co-accessibility         |

---

## 🔧 Notes
- **ArchR** is the only tool here that supports **multi-omic integration, TF motif activity scoring, and trajectory**.
- **Monocle3** is excellent for **RNA-only pseudotime and gene expression dynamics**.
- **Cicero** is best suited for **CRE connectivity**, not trajectory.