# Visualizations from completed ArchR
library(ggplot2)
# setting a custom plotting theme

# qplot.spec.coherency(p1 <- plotFragmentSizes(ArchRProj = projHeme1)
ArchRPalettes[["stallion2"]]

pal_1 <- c("Saline_Veh"="#9983BD", "CFA_Veh"="#FCBF6E", "CFA_ApAP"="#90D5E4", "CFA_3DDA"="#F37B7D")

p1 <- plotFragmentSizes(ArchRProj = projPAG2,
                        pal = pal_1, groupBy = "Sample")
p1 + theme(legend.position = "right", legend.direction = "vertical", 
         legend.text = element_text(size = 8, face = 'bold')) +
  guides(color = guide_legend(override.aes = list(size = 3)))

plotPDF(p1, name = "fragment_plot_custom", ArchRProj = projPAG2, addDOC = FALSE, width = 6, height = 6)

theme_set(theme_grey())


#umaps are first
ArchRPalettes[["ironMan"]]
umap_pal <- c("GABAergic Neurons" = "#C06CAB","Glut Neurons" = "#FCBF6E","Oligodendrocytes"="#90D5E4", "OPCs"="#88CFA4",
                     "Microglia"="#D24B27", "Astrocytes"="#0C727C","Ependymal Cells"="#FF0080")

umap_celltype<- plotEmbedding(
  ArchRProj = projPAG2,
  embedding = "UMAP",
  colorBy = "cellColData",
  name = "Clusters",
  size = 0.1,
  sampleCells = NULL,
  highlightCells = NULL,
  rastr = TRUE,
  quantCut = c(0.01, 0.99),
  pal=umap_pal,
  continuousSet = NULL,
  randomize = TRUE,
  keepAxis = FALSE,
  baseSize = 10,
  plotAs = NULL,
  threads = getArchRThreads(),
  logFile = createLogFile("plotEmbedding")) + theme(legend.position = "right", legend.direction = "vertical", 
                                                    legend.text = element_text(size = 8, face = 'bold')) +
  guides(color = guide_legend(override.aes = list(size = 3)))

umap_celltype

plotPDF(umap_celltype, name = "umap_celltype", ArchRProj = projPAG2, addDOC = FALSE, width = 8, height = 6)


# lets make heatmaps of the genes categorized into TF, Ion channels, GPCRs, and heatmap from a comparison
# we will decide what browser trakcs to plot

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

TF_df <- final_data[final_data$name %in% 
                          c('Egr1','Egr2','Elk4','Ezh2','Hdac10',
                            'Irf3','Irf6','Mzf1','Nfatc1','Nfatc4','Nfkb2',
                            'Notch1','Notch4','Pou6f1','Rfx2','Rreb1',
                            'Tcf4','Zfp36l2','Creb3','Grip1',
                            'Sox11','Zbtb7a','Foxo6','Camta1','Ebf4','Kcnip3',
                            'Ovol2','Pou4f3','Rxrg','Crebl2','Egr1','Atf3','Crem',
                            'Ddit3','Sox11','Bach1','Stat5a','Snai2','Fosl2','Arid5a','Sox5a'), ]

TF_df <- TF_df[,-1]

heatmap_palette <- paletteContinuous(set="whitePurple")

TF_heatmap <- pheatmap::pheatmap(TF_df, color = heatmap_palette, cluster_cols = FALSE, main = 'Transcription Factors', border_color = NA)


plotPDF(TF_heatmap, name = "tf_heatmap", ArchRProj = projPAG2, addDOC = FALSE, width = 4, height = 6)

ion_df <- final_data[final_data$name %in% 
                      c('Tpcn2','Slc9a3','Kcnj9','Trpv3','Best1','Fxyd6','Cacng4','Pkd2','Ubc','Ano4','Tmem206',
                        'Kcnd1','Kcnab2','Clcn5','Cacna1a','Scn11a','Tmtm38a','P2rx3','P2rx5','Anxa6','Asic2','Cacnb4',
                        'Garn2','Cacna2d3','Kcna4',
                        'Kcnq3','Grin1','Grik4','Tmem63c','Htr3b','Cacng5','Scn10a','Trpv2','Lrrc8c','Kcnq5',
                        'Clcn1','Grid2','Kcnf2','Gabrg2','Gabra1','Kcne3','Nalcn','Glrb','Trpm8','Kcnq2','Gria1',
                        'Cacna1b','Kcng1','Kcnma1','Asic3','Gabrb3','Gabrg1','Kcnh7','Kcnh1','Gabra2','Scn8a',
                        'Kcnv1','Scn1a','Scn3a','Ryr2','Hcn1','Scn9a','Cacng8','Kcnh5','Trpc5','Cacna2d1','Gabra5','Kcnk12',
                        'Kcnk13','Kcnq2','Kcnma1','Ano1','Kcnk16','Kcna1','Htr3a','Trpv1','Scn11a','Scn9a','Grik4','Kcnb1',
                        'Chrma6','Kcnb2','Gria4','Scn10a','Kcnv1','Kcng2','Piezo1','Piezo2','Prdm12' ), ]

ion_df <- ion_df[,-1]

heatmap_palette <- paletteContinuous(set="whitePurple")

ion_heatmap <- pheatmap::pheatmap(ion_df, color = heatmap_palette, cluster_cols = FALSE, main = 'Ion Channels', border_color = NA)


plotPDF(ion_heatmap, name = "ion_heatmap", ArchRProj = projPAG2, addDOC = FALSE, width = 4, height = 10)


gpcr_df <- final_data[final_data$name %in% 
                       c('Gabbr2', 'Lpar6', 'Grm7', 'Chrm3', 'Cnr1', 'Chrm2', 
                         'Gpr158', 'Oprm1', 'Gpr19', 'Ntsr2', 'Mrgprx3', 'Omg', 'Gpr149'), ]

gpcr_df <- gpcr_df[,-1]

heatmap_palette <- paletteContinuous(set="whitePurple")

gpcr_heatmap <- pheatmap::pheatmap(gpcr_df, color = heatmap_palette, cluster_cols = FALSE, main = 'GPCRs', border_color = NA)


plotPDF(gpcr_heatmap, name = "gpcr_heatmap", ArchRProj = projPAG2, addDOC = FALSE, width = 4, height = 6)


#custom browser tracks with new color palette and formating

sox11 <- plotBrowserTrack(
  ArchRProj = projPAG2,
  groupBy = "Sample",
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
  sizes = c(8, 1.2, 3, 2),
  geneSymbol = "Sox11",
  upstream = 50000,
  downstream = 50000,
  loops = getCoAccessibility(projPAG2),
  pal = pal_1,
  baseSize = 7,
  ylim = c(0,1),
  scTileSize = 0.5,
  scCellsMax = 100,
  borderWidth = 0.4,
  tickWidth = 0.4,
  facetbaseSize = 7,
  title = "Sox11",
  verbose = TRUE,
  logFile = createLogFile("plotBrowserTrack")
)
grid::grid.newpage()
grid::grid.draw(sox11$Sox11)

atf3 <- plotBrowserTrack(
  ArchRProj = projPAG2,
  groupBy = "Sample",
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
  sizes = c(8, 1.2, 3, 2),
  geneSymbol = "Atf3",
  upstream = 50000,
  downstream = 50000,
  loops = getCoAccessibility(projPAG2),
  pal = pal_1,
  baseSize = 7,
  ylim = c(0,1),
  scTileSize = 0.5,
  scCellsMax = 100,
  borderWidth = 0.4,
  tickWidth = 0.4,
  facetbaseSize = 7,
  title = "Atf3",
  verbose = TRUE,
  logFile = createLogFile("plotBrowserTrack")
)


grid::grid.newpage()
grid::grid.draw(atf3$Atf3)

scn9a <- plotBrowserTrack(
  ArchRProj = projPAG2,
  groupBy = "Sample",
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
  sizes = c(8, 1.2, 3, 2),
  geneSymbol = "Scn9a",
  upstream = 50000,
  downstream = 50000,
  loops = getCoAccessibility(projPAG2),
  pal = pal_1,
  baseSize = 7,
  ylim = c(0,1),
  scTileSize = 0.5,
  scCellsMax = 100,
  borderWidth = 0.4,
  tickWidth = 0.4,
  facetbaseSize = 7,
  title = "Scn9a",
  verbose = TRUE,
  logFile = createLogFile("plotBrowserTrack")
)


grid::grid.newpage()
grid::grid.draw(scn9a$Scn9a)

scn1a <- plotBrowserTrack(
  ArchRProj = projPAG2,
  groupBy = "Sample",
  plotSummary =  c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
  sizes = c(8, 1.2, 3, 2),
  geneSymbol = "Scn1a",
  upstream = 50000,
  downstream = 50000,
  loops = getCoAccessibility(projPAG2),
  pal = pal_1,
  baseSize = 7,
  ylim = c(0,1),
  scTileSize = 0.5,
  scCellsMax = 100,
  borderWidth = 0.4,
  tickWidth = 0.4,
  facetbaseSize = 7,
  title = "Scn1a",
  verbose = TRUE,
  logFile = createLogFile("plotBrowserTrack")
)


grid::grid.newpage()
grid::grid.draw(scn1a$Scn1a)

# Marker peaks heatmaps (Marker features are features that are unique to a specific cell grouping)
# which peaks are unique to an individual cluster or a small group of clusters
markersPeaks_allsamples <- getMarkerFeatures(
  ArchRProj = projPAG2, 
  useMatrix = "PeakMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList_peaks <- getMarkers(markersPeaks_allsamples, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList_peaks$CFA_3DDA

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markersPeaks_allsamples, 
  cutOff = "FDR <= 0.1 & Log2FC >= 3",
  transpose = TRUE,
  pal= paletteContinuous(set="whitePurple"),scaleRows = TRUE
)

draw(heatmapPeaks, heatmap_legend_side = "right", annotation_legend_side = "bot")

# marker features for gene matrix (unbiased identification of marker features for any given cell groupings)
# this can also be done for transcription factor motifs from chromVAR daviations
markers_genesallsamples <- getMarkerFeatures(
  ArchRProj = projPAG2, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon")
  
  
geneheatmap <- plotMarkerHeatmap(
  seMarker = markers_genesallsamples, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.75",
  transpose = TRUE,
  pal= paletteContinuous(set="whitePurple"),scaleRows = TRUE
)

draw(geneheatmap, heatmap_legend_side = "right", annotation_legend_side = "bot")

# plotting chromvar daviations
plotVarDev <- getVarDeviations(projPAG, name = "MotifMatrix", plot = TRUE)

motifs <- c("GATA1", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5")
markerMotifs <- getFeatures(projPAG, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs


# plotting umap for specific gene

scn1a_umap <- plotEmbedding(
  ArchRProj = projPAG2, 
  colorBy = "GeneScoreMatrix", 
  name = "Scn1a", # change this to be a list of genes we want
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL, # if you want imputation weights, replace with: getImputeWeights(projPAG)
  pal = paletteContinuous(set="purpleOrange"),
  rastr = TRUE
)

scn1a_umap


scn9a_umap <- plotEmbedding(
  ArchRProj = projPAG2, 
  colorBy = "GeneScoreMatrix", 
  name = "Scn9a", # change this to be a list of genes we want
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL, # if you want imputation weights, replace with: getImputeWeights(projPAG)
  pal = paletteContinuous(set="purpleOrange"),
  rastr = TRUE
)

scn9a_umap


sox11_umap <- plotEmbedding(
  ArchRProj = projPAG2, 
  colorBy = "GeneScoreMatrix", 
  name = "Sox11", # change this to be a list of genes we want
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL, # if you want imputation weights, replace with: getImputeWeights(projPAG)
  pal = paletteContinuous(set="purpleOrange"),
  rastr = TRUE
)

sox11_umap

atf3_umap <- plotEmbedding(
  ArchRProj = projPAG2, 
  colorBy = "GeneScoreMatrix", 
  name = "Atf3", # change this to be a list of genes we want
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL, # if you want imputation weights, replace with: getImputeWeights(projPAG)
  pal = paletteContinuous(set="purpleOrange"),
  rastr = TRUE
)

atf3_umap

# Motif footprinting for motifs identified in heatmap across all samples

getMatrixFromProject(projPAG2, useMatrix = "MotifMatrix")

markersMotifs_all <- getMarkerFeatures(
  ArchRProj = projPAG2, 
  useMatrix = "MotifMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon",
  useSeqnames="z"
)

heatmapEM <- plotMarkerHeatmap(seMarker = markersMotifs_all, cutOff = "MeanDiff >= 0.1 & Pval <= 0.05",
                               pal = paletteContinuous(set = "coolwarm"))

ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")








