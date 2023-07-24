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

# CHAPTER 12-14
projPAG4 <- addMotifAnnotations(ArchRProj = projPAG4, motifSet = "JASPAR2020", name = "Motif",
       species = getGenome(ArchRProj = projPAG4))

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
#==================================================


# run after integration with RNA Data
projPAG4 <- addPeak2GeneLinks(
  ArchRProj = projPAG4,
  reducedDims = "IterativeLSI"
)

p2g <- getPeak2GeneLinks(
  ArchRProj = projPAG4,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

p1 <- plotEmbedding(ArchRProj = projPAG4, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p1

head(projPAG4$cellNames)

trajectory <- c("Progenitor", "GMP", "Mono")

projPAG4 <- addTrajectory(
  ArchRProj = projPAG4,
  name = "MyeloidU",
  groupBy = "Clusters",
  trajectory = trajectory,
  embedding = "UMAP",
  force = TRUE
)
