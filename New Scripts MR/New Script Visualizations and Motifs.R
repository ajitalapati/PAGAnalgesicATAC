library(ArchR)
library(pheatmap)
library(SEtools)
library(presto)
set.seed(1)


projPAG2 <- loadArchRProject(path = "/Volumes/LaCie/PAG_scATAC_Analysis/ArchR_Outs/Save-ProjPAG2", force = FALSE, showLogo = TRUE)

# Arch 11.1/11.3 identifying marker peaks (CFA_Veh vs CFA_ApAP analysis)
# motifs enriched in marker peaks so peaks unique to individual groups
markersPeaks <- getMarkerFeatures(
  ArchRProj = projPAG2, 
  useMatrix = "PeakMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markersPeaks

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "Pval <= 0.001 & Log2FC >= 3",
  transpose = TRUE, pal = paletteContinuous(set = "comet"), nLabel = 5
)

PeakmarkerList <- getMarkers(markersPeaks, cutOff = "Pval <= 0.001 & Log2FC >= 1")
PeakmarkerList$CFA_3DDA
PeakmarkerList$CFA_ApAP
PeakmarkerList$CFA_Veh
PeakmarkerList$Saline_Veh

Peaklist_heatmap <- plotMarkerHeatmap(
  seMarker = PeakmarkerList, 
  cutOff = "Pval <= 0.001 & Log2FC >= 2",
  transpose = TRUE, pal = paletteContinuous(set = "comet"), nLabel = 5
)


heatmapPeaks

markerTest1 <- getMarkerFeatures(
  ArchRProj = projPAG2, 
  useMatrix = "PeakMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  bgdGroups = "Saline_Veh",
  useGroups = "CFA_Veh"
)

pv1 <- plotMarkers(seMarker = markerTest1, name = "CFA_ApAP", plotAs = "Volcano")
pv1 <- pv1 + scale_color_manual(values=c("blue", "gray", "red")) + 
  theme(legend.position = "bottom", legend.direction = "horizontal", 
        legend.text = element_text(size = 8, face = 'bold')) +
  guides(color = guide_legend(override.aes = list(size = 2)))
pv1

markerTest2 <- getMarkerFeatures(
  ArchRProj = projPAG2, 
  useMatrix = "PeakMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  bgdGroups = "CFA_Veh",
  useGroups = "CFA_ApAP"
)


# nothing shoed up as differntial maybe use brians volcano plots for this only

markerTest3 <- getMarkerFeatures(
  ArchRProj = projPAG2, 
  useMatrix = "PeakMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  bgdGroups = "CFA_Veh",
  useGroups = "CFA_3DDA"
)

markerTest4 <- getMarkerFeatures(
  ArchRProj = projPAG2, 
  useMatrix = "PeakMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  bgdGroups = "CFA_ApAP",
  useGroups = "CFA_3DDA"
)


# Motif Work for CFA_Veh vs CFA_ApAP will be doen for each comparison separately
# ArchR 12.1 using homer since jaspar only maps 10 TFs which was used in another paper with rat atac data
devtools::install_github("GreenleafLab/chromVARmotifs")

library(BSgenome.Rnorvegicus.UCSC.rn7)
library(TFBSTools)

getGenome(ArchRProj = projPAG2)

projPAG2 <- addMotifAnnotations(ArchRProj = projPAG2, name = "Motif", motifSet = "homer",
                               species = "BSgenome.Rnorvegicus.UCSC.rn7", force =TRUE)

# Motifs in peaks that are more accessible in CFA ApAP cells vs CFA-Veh
motifsUp1 <- peakAnnoEnrichment(
  seMarker = markerTest1,
  ArchRProj = projPAG2,
  peakAnnotation = "Motif",
  cutOff = "Pval <= 0.05 & Log2FC >= 0.5"
)

motifsUp1

df <- data.frame(TF = rownames(motifsUp1), mlog10Padj = assay(motifsUp1)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))


ggUp1 <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 2.5) +
  ggrepel::geom_label_repel(
    data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 3,
    nudge_x = 2,
    color = "black",
    fontface = "bold",
    max.overlaps = 20
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "purpleOrange"))

ggUp1

# Motifs in peaks that are more accessible in CFA_Veh than in CFA_ApAP
motifsDo1 <- peakAnnoEnrichment(
  seMarker = markerTest1,
  ArchRProj = projPAG2,
  peakAnnotation = "Motif",
  cutOff = "Pval <= 0.05 & Log2FC <= -0.5"
)

motifsDo1

df2 <- data.frame(TF = rownames(motifsDo1), mlog10Padj = assay(motifsDo1)[,1])
df2 <- df2[order(df2$mlog10Padj, decreasing = TRUE),]
df2$rank <- seq_len(nrow(df2))

head(df2)

ggDo1 <- ggplot(df2, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 2.5) +
  ggrepel::geom_label_repel(
    data = df2[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 3,
    nudge_x = 2,
    color = "black",
    fontface = "bold",
    max.overlaps = 20
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "purpleOrange"))

ggDo1

motifsUp2 <- peakAnnoEnrichment(
  seMarker = markerTest2,
  ArchRProj = projPAG2,
  peakAnnotation = "Motif",
  cutOff = "Pval <= 0.05 & Log2FC >= 0.5"
)

motifsUp2

df3 <- data.frame(TF = rownames(motifsUp2), mlog10Padj = assay(motifsUp2)[,1])
df3 <- df3[order(df3$mlog10Padj, decreasing = TRUE),]
df3$rank <- seq_len(nrow(df3))

head(df3)


ggUp2 <- ggplot(df3, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 2.5) +
  ggrepel::geom_label_repel(
    data = df3[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 3,
    nudge_x = 2,
    color = "black",
    fontface = "bold",
    max.overlaps = 20
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "purpleOrange"))

ggUp2

# Motifs in peaks that are more accessible in CFA_Veh than in CFA_ApAP
motifsDo2 <- peakAnnoEnrichment(
  seMarker = markerTest2,
  ArchRProj = projPAG2,
  peakAnnotation = "Motif",
  cutOff = "Pval <= 0.05 & Log2FC <= -0.5"
)

motifsDo2

df4 <- data.frame(TF = rownames(motifsDo2), mlog10Padj = assay(motifsDo2)[,1])
df4 <- df4[order(df4$mlog10Padj, decreasing = TRUE),]
df4$rank <- seq_len(nrow(df4))

head(df4)

ggDo2 <- ggplot(df4, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 2.5) +
  ggrepel::geom_label_repel(
    data = df4[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 3,
    nudge_x = 2,
    color = "black",
    fontface = "bold",
    max.overlaps = 20
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "purpleOrange"))

ggDo2


motifsUp3 <- peakAnnoEnrichment(
  seMarker = markerTest3,
  ArchRProj = projPAG2,
  peakAnnotation = "Motif",
  cutOff = "Pval <= 0.05 & Log2FC >= 0.5"
)

motifsUp3

df5 <- data.frame(TF = rownames(motifsUp3), mlog10Padj = assay(motifsUp3)[,1])
df5 <- df5[order(df5$mlog10Padj, decreasing = TRUE),]
df5$rank <- seq_len(nrow(df5))

head(df5)

ggUp3 <- ggplot(df5, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 2.5) +
  ggrepel::geom_label_repel(
    data = df5[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 3,
    nudge_x = 2,
    color = "black",
    fontface = "bold",
    max.overlaps = 20
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "purpleOrange"))

ggUp3

# Motifs in peaks that are more accessible in CFA_Veh than in CFA_ApAP
motifsDo3 <- peakAnnoEnrichment(
  seMarker = markerTest3,
  ArchRProj = projPAG2,
  peakAnnotation = "Motif",
  cutOff = "Pval <= 0.05 & Log2FC <= -0.5"
)

motifsDo3

df6 <- data.frame(TF = rownames(motifsDo3), mlog10Padj = assay(motifsDo3)[,1])
df6 <- df6[order(df6$mlog10Padj, decreasing = TRUE),]
df6$rank <- seq_len(nrow(df6))

head(df6)

ggDo3 <- ggplot(df6, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 2.5) +
  ggrepel::geom_label_repel(
    data = df6[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
    size = 3,
    nudge_x = 2,
    color = "black",
    fontface = "bold",
    max.overlaps = 20
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "purpleOrange"))

ggDo3




# identifying motifs using getmarker features to plot to see if DAR are correlated with motif enrichment
# we actually ownt use these because the motif matrix is the varDev matrix 
getMatrixFromProject(projPAG2, useMatrix = "MotifMatrix")

markersMotifs_all <- getMarkerFeatures(
  ArchRProj = projPAG2, 
  useMatrix = "MotifMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useSeqnames="z"
)

heatmapEM <- plotMarkerHeatmap(seMarker = markersMotifs_all, cutOff = "abs(MeanDiff) >= 0.1 & Pval <= 0.05",
                               pal = paletteContinuous(set = "comet"))

ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")

# subsetting motif deviation matrix

plotVarDev <- getVarDeviations(projPAG2, name = "MotifMatrix", plot = TRUE)

plotVarDev + geom_point(size = 2.5) + theme_ArchR() + 
  scale_color_gradientn(colors = paletteContinuous(set = "purpleOrange"))

motifs <- c("Sox2", "Sox15", "Sox10", "Sox6", "Sox9", "Sox4", "ELF5","ETS1","RORgt")
markerMotifs <- getFeatures(projPAG2, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- markerMotifs[markerMotifs %ni% "z:Sox9.HMG_265"]
markerMotifs


p <- plotGroups(ArchRProj = projPAG2, 
                groupBy = "Sample", 
                colorBy = "MotifMatrix", 
                name = markerMotifs,
                imputeWeights = getImputeWeights(projPAG2),
                pal = magma_pal
)

p$`z:Sox6.HMG_264`
p$`z:Sox4.HMG_263`
p$`z:Sox2.HMG_261`
p$`z:Sox15.HMG_260`
p$`z:Sox10.HMG_259`
p$`z:RORgt.NR_246`
p$`z:ETS1.ETS_76`


# footprinting is next
motifPositions <- getPositions(projPAG2)
motifPositions
motifPositions@partitioning@NAMES


motifs <- c("Sox2", "Sox15", "Sox10", "Sox6", "Sox9", "Sox4", "ELF5","ETS1","RORgt")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs

library(BSgenome.Rnorvegicus.UCSC.rn7)
projPAG2 <- addGroupCoverages(ArchRProj = projPAG2, groupBy = "Sample")

library(BSgenome.Rnorvegicus.UCSC.rn7)
seFoot <- getFootprints(
  ArchRProj = projPAG2, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Sample"
)

# Finally positive TF regulator identification

seGroupMotif <- getGroupSE(ArchRProj = projPAG2, useMatrix = "MotifMatrix", groupBy = "Sample")

seGroupMotif

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

getFeatures(projPAG2, 'MotifMatrix')

getFeatures(projPAG2, 'GeneScoreMatrix')


getAvailableMatrices(projPAG2)
corGSM_MM <- correlateMatrices(
  ArchRProj = projPAG2,
  useMatrix1 = "GeneScoreMatrix",
  useMatrix2 = "MotifMatrix",
  removeFromName2 = "dot",
  reducedDims = "IterativeLSI"
)

corGSM_MM

corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM$MotifMatrix_name, rowData(seZ)$name), "maxDelta"]

corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,"MotifMatrix_name"]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])

label_data <- data.frame(corGSM_MM)[data.frame(corGSM_MM)$TFRegulator == "YES",]

library(ggrepel)
p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point(size=2) + 
  theme_ArchR() + theme(axis.text = element_text(face="bold")) +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="#D24B27")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05) 
  ) +  geom_label_repel(data = label_data, aes(label = MotifMatrix_matchName), box.padding = 0.2, 
                       point.padding = 0.2, nudge_y = 0.1, label.size = .1, fontface="bold") # Adjust some_column to your desired label
p


TF_labs <- corGSM_MM$MotifMatrix_matchName[which(corGSM_MM$TFRegulator == "YES")]
TF_labs
str(TF_labs)


BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
corGSM_MM@listData[["maxDelta"]]

ArchRPalettes[[set ="purpleOrange"]]

volano_plot <- EnhancedVolcano(corGSM_MM,
                               lab = corGSM_MM$GeneScoreMatrix_matchName,
                               selectLab = "SOX10",
                               pCutoff = 0,
                               FCcutoff = 0,
                               x = "cor",
                               y = "maxDelta",
                               pointSize = 3,
                               legendLabSize = 10,
                               xlim = c(-1, 1 ),
                               ylim = c(0, 2.5),
                               labSize = 2.0,
                               xlab = bquote("Correlation To Gene Score"),
                               ylab = bquote("Max TF Motif Delta"),
                               drawConnectors = TRUE,
                               boxedLabels = TRUE,
                               lengthConnectors = 0.5,
                               arrowheads = FALSE,
                               maxoverlapsConnectors = Inf
                               )
volano_plot



# custom motif seq logos from homer PWM

getPeakAnnotation(projPAG2, "Motif")$motifSummary

library(chromVARmotifs)
library(ggseqlogo)
library(seqLogo)
data("homer_pwms") #motifs from HOMER

print(names(homer_pwms))

PWMatrixToProbMatrix <- function(x){
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
  m <- t(t(m)/colSums(m))
  m
}

ProbMatrices <- lapply(homer_pwms, PWMatrixToProbMatrix)
lapply(ProbMatrices, colSums) %>% range

seqLogo(ProbMatrices$`Sox10(HMG)/SciaticNerve-Sox3-ChIP-Seq(GSE35132)/Homer`)

seqLogo(ProbMatrices$`Sox6(HMG)/Myotubes-Sox6-ChIP-Seq(GSE32627)/Homer`)

seqLogo(ProbMatrices$`ETS1(ETS)/Jurkat-ETS1-ChIP-Seq(GSE17954)/Homer`)

seqLogo(ProbMatrices$`RORgt(NR)/EL4-RORgt.Flag-ChIP-Seq(GSE56019)/Homer`)

plotTSSEnrichment(
  ArchRProj = projPAG2,
  groupBy = "Sample",
  chromSizes = getChromSizes(projPAG2),
  TSS = getTSS(projPAG2),
  flank = 2000,
  norm = 100,
  smooth = 11,
  pal = NULL,
  returnDF = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("plotTSSEnrichment")
)

rm(df_TSS)

# new umap with cell counts on side
p3 <- plotEmbedding(ArchRProj = projPAG2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", pal=umap_pal)  + 
  theme(plot.margin = margin(0, 0, 5, 5)) + theme(legend.position = "bottom", legend.direction = "horizontal", 
                                                  legend.text = element_text(size = 8, face = 'bold')) +
  guides(color = guide_legend(override.aes = list(size = 3)))
p3

counter_hashmap <- table(projPAG2$Clusters)

# Convert the counter hashmap to a data frame
data_df <- data.frame(values = names(counter_hashmap), counts = as.vector(counter_hashmap))

# Sort the data frame by 'values' in alphanumerical order
data_df$values <- factor(data_df$values, levels = sort(unique(data_df$values)))

# Plot the histogram using ggplot2
library(ggpubr)
count_hist <- ggplot(data_df, aes(x = counts, y = values, fill = values)) +
  geom_bar(stat = "identity", fill = umap_pal, color = "black", position = "dodge") + theme_pubr() +
  labs(title = "Number of Cells",
       x = "Frequency", y = "")

count_hist

widths <- c(3, 1.5)

# Combine the plots side by side using cowplot
final_plot <- plot_grid(p3, count_hist + geom_text(aes(label = "")), ncol = 2, rel_widths = widths)

final_plot


