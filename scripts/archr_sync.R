setwd("~/Documents/PAGAnalgesicATAC")
library(ArchR)
library(parallel)
addArchRThreads(threads = 18) 
set.seed(1) 

projPAG <- loadArchRProject(path = "./PAG_ATAC_Directory4", force = FALSE, showLogo = FALSE)

markersGS <- getMarkerFeatures(
  ArchRProj = projPAG, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 1000
)

# Get markers for each of the clusters
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 2")