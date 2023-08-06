library(ArchR)
library(pheatmap)
BiocManager::install("SEtools")
library(SEtools)
set.seed(1)


projPAG2 <- readRDS("Save-ArchR-Project.rds")

# ArchR 10.4 ading a peak matrix after calling peaks
projPAG2 <- addPeakMatrix(projPAG2)

getAvailableMatrices(projPAG2)

table(projPAG2$Sample)

getCellColData(projPAG2)

devtools::install_github("immunogenomics/presto")
library(presto)

# clustering and parameters
projPAG2 <- addIterativeLSI(
  ArchRProj = projPAG2,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30
)

# Now we can begin clustering which uses the Seurat FindCLusters function
projPAG2 <- addClusters(
  input = projPAG2,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8
)

# marker genes for cell types
markersGS <- getMarkerFeatures(
  ArchRProj = projPAG2, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Get markers for each of the clusters
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1")

marker_dfC2 <- as.data.frame(markerList$C2)

ranked_dfC1 <- as.data.frame(markerList$C1[order(-markerList$C1$Log2FC), ])

marker_dfC3 <- as.data.frame(markerList$C3)

marker_dfC4 <- as.data.frame(markerList$C4)

marker_dfC5 <- as.data.frame(markerList$C5)

marker_dfC6 <- as.data.frame(markerList$C6)

marker_dfC7 <- as.data.frame(markerList$C7)

marker_dfC8 <- as.data.frame(markerList$C8)

marker_dfC9 <- as.data.frame(markerList$C9)

marker_dfC10 <- as.data.frame(markerList$C10)

marker_dfC11 <- as.data.frame(markerList$C11)

marker_dfC12 <- as.data.frame(markerList$C12)

marker_dfC13 <- as.data.frame(markerList$C13)

marker_dfC14 <- as.data.frame(markerList$C14)

marker_dfC15 <- as.data.frame(markerList$C15)

marker_dfC16 <- as.data.frame(markerList$C16)

# ArchR8.3

projPAG2$Clusters <- mapLabels(projPAG2$Clusters,
                               newLabels = c("Microglia", "Oligodendrocytes","Oligodendrocytes",
                              "Oligodendrocytes","Oligodendrocytes", "Astrocytes", 
                              "OPCs","Ependymal Cells","Glut Neurons","GABAergic Neurons",
                              "GABAergic Neurons", "Glut Neurons","Glut Neurons",
                              "Glut Neurons","Glut Neurons","Glut Neurons"), 
                               oldLabels = c("C1","C2","C3","C4", "C5", "C6", "C7", "C8","C9",
                              "C10","C11","C12", "C13","C14","C15","C16"))

projPAG2$Clusters
pal <- c("GABAergic Neurons" = "#C06CAB","Glut Neurons" = "#E6C2DC","Oligodendrocytes"="#90D5E4", "OPCs"="#F37B7D",
        "Microglia"="#D24B27", "Astrocytes"="#0C727C","Ependymal Cells"="#A6CDE2")

ArchRPalettes$paired

p1 <- plotEmbedding(projPAG2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP", baseSize = 12, discreteSet = "ironMan")
p1

plotPDF(p1, name = "Plot-UMAP-Remap-Clusters.pdf", ArchRProj = projPAG2, addDOC = FALSE, width = 6, height = 6)

# plotting customarkerFeatures using genes of interest

markersGS <- getMarkerFeatures(
  ArchRProj = projPAG2, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  maxCells = 1000
)

# making a data frame from a single cell experiment object
assay_data <- assay(markersGS)
assay_data_df <- as.data.frame(assay_data)
# Extract the row and column annotations
row_data <- as.data.frame(rowData(markersGS))

col_data <- as.data.frame(colData(markersGS))

# Combine with the assay data
final_data <- cbind(row_data, assay_data_df)

final_data <- final_data[,-(1:4)]

final_data <- final_data[,-2]


rownames(final_data) <- final_data$name

subset_df <- final_data[final_data$name %in% 
                          c("Sox11", "Atf3", "Gallr", "Gabbr2", "Grm7", "Gpr158", "Scn9a", 
                            "Scn10", "Kcnip3","Chrm2","Chrm3","Cnr1","Jun","Fos","Il6","Il1a","Sprr1a","Adcyap1","Pdyn",
                            "Adra2a","Kcnj3","Drd2","Scnn11a","Scn8a","Gria2","Ephb1","Vegfa","Aqp1","Tac1"), ]

subset_df <- subset_df[,-1]

my_palette <- ArchRPalettes[["whitePurple"]]

heatmapcustom <- pheatmap::pheatmap(subset_df, color = my_palette, cluster_cols = FALSE)

plotPDF(heatmapcustom, name = "custom_gene_heatmap.pdf", ArchRProj = projPAG2, addDOC = FALSE, width = 6, height = 6)


# Heatmap generic
heatmap_GS1 <- plotMarkerHeatmap(
  seMarker = markersGS,
  log2Norm = TRUE,
  scaleTo = 10^4,
  scaleRows = TRUE,
  plotLog2FC = FALSE,
  limits = c(-2, 2),
  grepExclude = NULL,
  pal = paletteContinuous(set="whitePurple"),
  binaryClusterRows = TRUE,
  clusterCols = TRUE,
  labelMarkers = specificGenes,
  nLabel = 15,
  nPrint = 15,
  labelRows = FALSE,
  returnMatrix = FALSE,
  transpose = FALSE,
  invert = FALSE,
  logFile = createLogFile("plotMarkerHeatmap")
)

plotPDF(heatmap_GS1, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = projPAG2, addDOC = FALSE)



# Arch 11.1/11.3 identifying marker peaks (CFA_Veh vs CFA_ApAP analysis)
markersPeaks_Veh_vs_CFA <- getMarkerFeatures(
  ArchRProj = projPAG2, 
  useMatrix = "PeakMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  bgdGroups = "Saline_Veh",
  useGroups = "CFA_Veh"
)

pma <- plotMarkers(seMarker = markersPeaks_Veh_vs_CFA, name = "CFA_ApAP", plotAs = "MA")
pma

pv <- plotMarkers(seMarker = markersPeaks_Veh_vs_CFA, name = "CFA_ApAP", plotAs = "Volcano")
pv

# ArchR 11.2.3
p <- plotBrowserTrack(
  ArchRProj = projPAG2, 
  groupBy = "Sample", 
  geneSymbol = c("FOS"),
  features =  getMarkers(markersPeaks_Veh_vs_CFA, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["CFA_Veh"],
  upstream = 50000,
  downstream = 50000,
)

grid::grid.newpage()
grid::grid.draw(p$Fos)

# Motif Work for CFA_Veh vs CFA_ApAP will be doen for each comparison seperately
# ArchR 12.1 using homer which was used in another paper with rat atac data
devtools::install_github("GreenleafLab/chromVARmotifs")

projPAG2 <- addMotifAnnotations(ArchRProj = projPAG2, motifSet = "homer", name = "Motif")

motifsUp <- peakAnnoEnrichment(
  seMarker = markersPeaks2v4,
  ArchRProj = projPAG2,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

motifsUp

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

head(df)

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
  scale_color_gradientn(colors = paletteContinuous(set = "beach"))

ggUp

motifsDo <- peakAnnoEnrichment(
  seMarker = markersPeaks2v4,
  ArchRProj = projPAG2,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)

motifsDo

df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

head(df)

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
  scale_color_gradientn(colors = paletteContinuous(set = "beach"))

ggDo

saveArchRProject(ArchRProj = projPAG2, outputDirectory = "/Volumes/LaCie/PAG_scATAC_Analysis/ArchR_Outs/Save-ProjPAG2", load = FALSE)

# adding motifc deviations using ChromVAR

if("Motif" %ni% names(projPAG2@peakAnnotation)){
  projPAG2 <- addMotifAnnotations(ArchRProj = projPAG2, motifSet = "homer", name = "Motif")
}

projPAG2 <- addBgdPeaks(projPAG2)


projPAG2 <- addDeviationsMatrix(
  ArchRProj = projPAG2, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(projPAG2, name = "MotifMatrix", plot = TRUE)
plotVarDev

# clean up

rm(marker_df1,marker_dfC2,marker_dfC3, marker_dfC10, marker_dfC11, marker_df, marker_dfC12, marker_dfC13, 
   marker_dfC14, marker_dfC15,marker_dfC16, marker_dfC4, marker_dfC5, marker_dfC6,marker_dfC7, marker_dfC8, marker_dfC9)
rm(pag.data, matrix_genes, geneScoreMatrix, t, subsetSE, ranked_dfC1, ranked_dfC3, df)


# ArchR 14.1
# Motif Fotprinting

motifPositions <- getPositions(projPAG2)
motifPositions
motifPositions@partitioning@NAMES


motifs <- c("Fos","Jun", "Egr1","Mecp3", "Sox11", "Atf3")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs

library(BSgenome.Rnorvegicus.UCSC.rn7)
seFoot <- getFootprints(
  ArchRProj = projPAG2, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Clusters"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = projPAG2, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = projPAG2, 
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = projPAG2, 
  normMethod = "None",
  plotName = "Footprints-No-Normalization",
  addDOC = FALSE,
  smoothWindow = 5
)


# skiping 14.3 for now until we decide what we want
# moving to co-accessibility in chapter 15
projPAG2 <- addCoAccessibility(
  ArchRProj = projPAG2,
  reducedDims = "IterativeLSI"
)

cA <- getCoAccessibility(
  ArchRProj = projPAG2,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = FALSE
)
cA

cA <- getCoAccessibility(
  ArchRProj = projPAG2,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = TRUE
)
cA

browsertrackgenes <- c("Sox11", 
                       "Atf3",
                       "Gallr", 
                       "Gabbr2", 
                       "Grm7", 
                       "Gpr158")


p <- plotBrowserTrack(
  ArchRProj = projPAG2, 
  groupBy = "Sample", 
  geneSymbol = browsertrackgenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getCoAccessibility(projPAG2),
  discreteSet = "darjeeling"
)

grid::grid.newpage()
grid::grid.draw(p$Gpr158)

# custom modification of plotting browser tracks with formatting from ggplot
plotBrowserTrack(
  ArchRProj = projPAG2,
  region = NULL,
  groupBy = "Sample",
  useGroups = NULL,
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
  sizes = c(10, 1.5, 3, 4),
  features = getPeakSet(ArchRProj),
  loops = getCoAccessibility(ArchRProj),
  geneSymbol = NULL,
  useMatrix = NULL,
  log2Norm = TRUE,
  upstream = 50000,
  downstream = 50000,
  tileSize = 250,
  minCells = 25,
  normMethod = "ReadsInTSS",
  threads = getArchRThreads(),
  ylim = NULL,
  pal = NULL,
  discreteSet = "darjeeling",
  baseSize = 7,
  scTileSize = 0.5,
  scCellsMax = 100,
  borderWidth = 0.4,
  tickWidth = 0.4,
  facetbaseSize = 7,
  geneAnnotation = getGeneAnnotation(ArchRProj),
  title = "",
  verbose = TRUE,
  logFile = createLogFile("plotBrowserTrack")
)

# custom plot groups function allows to plot any specified groups for visual comparison
plotGroups(
  ArchRProj = NULL,
  groupBy = "Sample",
  colorBy = "colData",
  name = "TSSEnrichment",
  imputeWeights = if (!grepl("coldata", tolower(colorBy[1])))
    getImputeWeights(ArchRProj),
  maxCells = 1000,
  quantCut = c(0.002, 0.998),
  log2Norm = NULL,
  pal = NULL,
  discreteSet = "darjeeling",
  ylim = NULL,
  size = 0.5,
  baseSize = 6,
  ratioYX = NULL,
  ridgeScale = 2,
  plotAs = "ridges")

# chapter 15.4 ArchR
seGroupMotif <- getGroupSE(ArchRProj = projPAG2, useMatrix = "MotifMatrix", groupBy = "Clusters")

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

p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="#C06CAB")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )

p
