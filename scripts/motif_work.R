setwd("~/Documents/PAGAnalgesicATAC")
library(ArchR)
library(parallel)
addArchRThreads(threads = 22) 

library(BSgenome.Rnorvegicus.UCSC.rn7)
library(TxDb.Rnorvegicus.UCSC.rn7.refGene)
library(org.Rn.eg.db)

library(Cairo)

projPAG <- loadArchRProject(path = "./PAG_ATAC_Directory4", force = FALSE, showLogo = FALSE)

markersPeaks <- getMarkerFeatures(
  ArchRProj = projPAG5, 
  useMatrix = "PeakMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon", 
  bgdGroups = c("CFA_Veh"),
  useGroups = c("Saline_Veh", "CFA_Veh", "CFA_ApAP", "CFA_3DDA"))


heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE,
  returnMatrix = TRUE,
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
  ArchRProj = projPAG, 
  useMatrix = "PeakMatrix",
  groupBy = "Sample",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  bgdGroups = "Saline_Veh",
  useGroups = "CFA_Veh"
)

pma <- plotMarkers(seMarker = markerTest, name = "CFA_Veh", plotAs = "Volcano", cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
pma

getCellColData(projPAG4)

# CHAPTER 12-14
library(JASPAR2020)
library(TFBSTools)

PWM <- getMatrixSet(
  x = JASPAR2020,
  opt = list(
    collection = "CORE",
    tax_group = "vertebrates",
    species = "9606",
    all_versions = FALSE,
    matrixtype = "PWM"
  )
)

projPAG <- addMotifAnnotations(ArchRProj = projPAG, motifSet = "homer", name = "Motif", force=TRUE)

motifsUp <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = projPAG,
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
  ArchRProj = projPAG,
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

markersPeaks <- getMarkerFeatures(
  ArchRProj = projPAG5, 
  useMatrix = "PeakMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon", 
  bgdGroups = c("CFA_Veh"),
  useGroups = c("CFA_3DDA", "CFA_ApAP"))

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = projPAG5,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

data.frame(assays(enrichMotifs)$Enrichment)

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)

ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

#12.3


#13.1

projPAG <- addBgdPeaks(projPAG)
projPAG <- addDeviationsMatrix(
  ArchRProj = projPAG, 
  peakAnnotation = "Motif",
)

plotVarDev <- getVarDeviations(projPAG, name = "MotifMatrix", plot = TRUE)
plotVarDev

motifs <- c("ZEB1", "SOX9", "SPI1", "SPIB", "RFX2")

ZEB1 <- c("ZEB1")


markerMotifs <- getFeatures(projPAG5, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs

markerMotifs <- c("z:ZEB1_897", "z:SOX9_1382", "z:SPI1_1022", "z:SPIB_1370", "z:RFX2_1480", "z:RFX2_1488")

markerMotifs <- c("z:EGR1_752", "z:RREB1_877", "z:CREB3_1611", "z:SOX11_1025", "z:ATF3_311")

p <- plotGroups(ArchRProj = projPAG5, 
                groupBy = "Sample", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(projPAG5)
)
p

p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }else{
    p[[x]] + guides(color = FALSE, fill = FALSE) + 
      theme_ArchR(baseSize = 6) +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
      theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
      ) + ylab("")
  }
})
do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))




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
  ArchRProj = projPAG5, 
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

