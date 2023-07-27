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
  ArchRProj = projPAG5, 
  useMatrix = "PeakMatrix",
  groupBy = "Sample",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "Saline_Veh",
  bgdGroups = "CFA_Veh"
)

pma <- plotMarkers(seMarker = markerTest, name = "Saline_Veh", plotAs = "Volcano")
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

temp <- c("Egr1","Egr2","Elk4","Ezh2","Hdac10","Irf3","Irf6","Mzf1","Nfatc1","Nfatc4","Nfkb2","Notch1","Notch4","Pou6f1","Rfx2","Rreb1","Tcf4","Zfp36l2","Creb3","Grip1")
tfs <- lapply(temp, toupper)

motifs <- c("TCF4", "GRIP1", "NOTCH1", "RREB1", "CAMTA1", "CREBL2", "CREM", "KCNIP3")


markerMotifs <- unlist(lapply(tfs, function(x) grep(x, names(motifPositions), value = TRUE)))

markerMotifs

projPAG4 <- addGroupCoverages(ArchRProj = projPAG4, groupBy = "Sample", force=TRUE)

seFoot <- getFootprints(
  ArchRProj = projPAG4, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Sample"
)

#14.2
plotFootprints(
  seFoot = seFoot,
  ArchRProj = projPAG4, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias-By-Sample-3",
  addDOC = FALSE,
  smoothWindow = 5
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = projPAG4, 
  normMethod = "None",
  plotName = "Footprints-By-Sample-3-No-Norm",
  addDOC = FALSE,
  smoothWindow = 5
)

#14.3
projPAG4 <- addGroupCoverages(ArchRProj = projPAG4, groupBy = "Sample")

seTSS <- getFootprints(
  ArchRProj = projPAG4, 
  positions = GRangesList(TSS = getTSS(projPAG4)), 
  groupBy = "Sample",
  flank = 2000
)

plotFootprints(
  seFoot = seTSS,
  ArchRProj = projPAG4, 
  normMethod = "None",
  plotName = "TSS-No-Normalization",
  addDOC = FALSE,
  flank = 2000,
  flankNorm = 100
)
saveArchRProject(ArchRProj = projPAG4, outputDirectory = "PAG_ATAC_Directory4", load = FALSE)

