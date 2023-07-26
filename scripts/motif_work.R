setwd("~/Documents/PAGAnalgesicATAC")
library(ArchR)
library(parallel)
addArchRThreads(threads = 22) 

library(BSgenome.Rnorvegicus.UCSC.rn7)
library(TxDb.Rnorvegicus.UCSC.rn7.refGene)
library(org.Rn.eg.db)

library(Cairo)

projPAG4 <- loadArchRProject(path = "./PAG_ATAC_Directory3", force = FALSE, showLogo = FALSE)

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
  groupBy = "Sample",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CFA_3DDA",
  bgdGroups = "CFA_ApAP"
)

pma <- plotMarkers(seMarker = markerTest, name = "CFA_3DDA", plotAs = "Volcano")
pma

getCellColData(projPAG4)

# CHAPTER 12-14
library(JASPAR2020)

projPAG4 <- addMotifAnnotations(ArchRProj = projPAG4, motifSet = "encode", name = "Motif")

motifsUp <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = projPAG4,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggUp

motifsDo <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = projPAG4,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)

df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 1.5,
    nudge_x = 2,
    color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggDo

#12.2
enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = projPAG4,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)

ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

#13.1

projPAG4 <- addBgdPeaks(projPAG4)
projPAG4 <- addDeviationsMatrix(
  ArchRProj = projPAG4, 
  peakAnnotation = "Motif",
)

plotVarDev <- getVarDeviations(projPAG4, name = "MotifMatrix", plot = TRUE)

#14.1
motifPositions <- getPositions(projPAG4)
motifs <- c("GATA1", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs

projPAG4 <- addGroupCoverages(ArchRProj = projPAG4, groupBy = "Clusters", force=TRUE)

seFoot <- getFootprints(
  ArchRProj = projPAG4, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Clusters"
)

#14.2
plotFootprints(
  seFoot = seFoot,
  ArchRProj = projPAG4, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias-By-Cluster",
  addDOC = FALSE,
  smoothWindow = 5
)

