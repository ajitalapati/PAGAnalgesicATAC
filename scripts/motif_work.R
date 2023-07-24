setwd("/work/aalapa/pag")
library(ArchR)
library(parallel)
addArchRThreads(threads = 18) 

library(BSgenome.Rnorvegicus.UCSC.rn7)
library(TxDb.Rnorvegicus.UCSC.rn7.refGene)
library(org.Rn.eg.db)

library(Cairo)

projPAG4 <- loadArchRProject(path = "./PAG_ATAC_Directory3", force = FALSE, showLogo = TRUE)

markersPeaks <- getMarkerFeatures(
  ArchRProj = projPAG4, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

heatmapPeaks

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")
plotPDF(heatmapPeaks, name = "heatmapPeaks.pdf", ArchRProj = projPAG4, addDOC = FALSE, width = 5, height = 5)

#section 11.1 in manual
markerTest <- getMarkerFeatures(
  ArchRProj = projPAG4, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList

markerTest <- getMarkerFeatures(
  ArchRProj = projPAG4, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C1",
  bgdGroups = "C2"
)

pma <- plotMarkers(seMarker = markerTest, name = "C1", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
pma

# STILL NEEDS COMPLETION
projPAG4 <- addMotifAnnotations(ArchRProj = projPAG4, motifSet = "cisbp", name = "Motif")

projPAG4 <- addCoAccessibility(
  ArchRProj = projPAG4,
  reducedDims = "IterativeLSI"
)

cA <- getCoAccessibility(
  ArchRProj = projPAG4,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = FALSE
)


