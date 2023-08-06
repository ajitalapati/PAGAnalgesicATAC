setwd("~/Documents/PAGAnalgesicATAC")
library(ArchR)
library(parallel)
addArchRThreads(threads = 18) 
set.seed(1) 

projPAG <- loadArchRProject(path = "./PAG_ATAC_Directory3", force = FALSE, showLogo = TRUE)

projPAG2 <- filterDoublets(projPAG)

projPAG2

p1 <- plotGroups(
  ArchRProj = projPAG, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "ridges"
)
p1

plotPDF(p1, name = "QC-PlotGroups-Ridge.pdf", ArchRProj = projPAG, addDOC = FALSE, width = 5, height = 5)


p2 <- plotFragmentSizes(ArchRProj = projPAG)
p2

p3 <- plotTSSEnrichment(ArchRProj = projPAG)
p3

# to save vectorized versions of these plots
plotPDF(p2,p3, name = "QC-Sample-FragSizes-TSSProfile.pdf", ArchRProj = projPAG, addDOC = FALSE, width = 5, height = 5)

# Iterative LSI will do further batch correction with Harmony through ArchR if necessary
# received the Cairo rasterization warning will load Cairo if it effects anything downstream
# does not seem to change anyting so far
projPAG <- addIterativeLSI(
  ArchRProj = projPAG,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force=TRUE
)

# Now we can begin clustering which uses the Seurat FindCLusters function
projPAG <- addClusters(
  input = projPAG,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,
  force= TRUE
)

# access the clustering information
projPAG$Clusters

# get number of cells present in each cluster
table(projPAG$Clusters)

#creating a cluster confusion matrix to better understand the distribution of samples
# in thcluster

cM <- confusionMatrix(paste0(projPAG2$Clusters), paste0(projPAG2$Sample))

# this confusion matrix can then be plotted as a heatmao using pheatmap
library(pheatmap)
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)
p

plotPDF(p, name = "Sample-cluster-distribution-heatmap.pdf", ArchRProj = projPAG, addDOC = FALSE, width = 4, height = 6)

# Time to make umap since dimesnionality reduction has been performed
projPAG <- addUMAP(
  ArchRProj = projPAG, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)



# plotting the embedding
p1 <- plotEmbedding(ArchRProj = projPAG, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p1

p3 <- plotEmbedding(ArchRProj = projPAG, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3

plotPDF(p1,p2, name = "umap-colored-sample-cluster.pdf", ArchRProj = projPAG, addDOC = FALSE, width = 5, height = 5)


# Now time to compute gene score matrices

devtools::install_github("immunogenomics/presto")

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
markerList$C6$name

sampleSE <- getGroupSE(
  ArchRProj = projPAG,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Sample")

sub <- c("Sox11", "Atf3", "Gal1r", "Gabbr2", "Grm7", "Gpr158", "Chrm2", "Chrm3", "Cnr1", "Oprl1", "Oprm1", "Mrgpre")
subsetSE <- markersGS[which(rowData(markersGS)$name %in% sub),] 

heatmapGS <- plotMarkerHeatmap(
  seMarker = markersGS,
  labelMarkers = sub, #can set to marker genes for labelling once defined
  transpose = TRUE,
  pal = ArchRPalettes$coolwarm
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

sechm(subsetSE, features=sub, assayName = "Log2FC")

saveArchRProject(ArchRProj = projPAG2, outputDirectory = "PAG_ATAC_Directory2", load = FALSE)

p <- plotBrowserTrack(
  ArchRProj = projPAG2, 
  groupBy = "Clusters", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000
)

grid::grid.newpage()
grid::grid.draw(p$CD14)

plotPDF(plotList = p, 
        name = "Plot-Tracks-Marker-Genes.pdf", 
        ArchRProj = projPAG2, 
        addDOC = FALSE, width = 5, height = 5)

#extra plots for Dr. Bazan 8/3/2023
p2 <- plotEmbedding(ArchRProj = projPAG, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2

projPAG <- addTSNE(
  ArchRProj = projPAG, 
  reducedDims = "IterativeLSI", 
  name = "TSNE", 
  perplexity = 30
)
p1 <- plotEmbedding(ArchRProj = projPAG, colorBy = "cellColData", name = "Clusters", embedding = "TSNE")
p1

tf_list <- c(
  "Egr1", "Egr2", "Elk4", "Ezh2", "Hdac10", "Irf3", "Irf6", "Mzf1", "Nfatc1", "Nfatc4", 
  "Nfkb2", "Notch1", "Notch4", "Pou6f1", "Rfx2", "Rreb1", "Tcf4", "Zfp36l2", "Creb3", "Grip1",
  "Sox11", "Zbtb7a", "Camta1", "Ebf4", "Kcnip3", "Ovol2", "Pou4f3", "Rxrg", "Crebl2", 
  "Egr1", "Atf3", "Crem", "Ddit3", "Sox11", "Bach1", "Stat5a", "Snai2", "Fosl2", "Arid5a"
)

markerList$C1

# Step 1: Rank the data frame by the column "Log2FC" in descending order
ranked_df <- markerList$C1[order(-markerList$C1$Log2FC), ]

# Step 2: Get the top 9 rows from the ranked data frame
top_9_rows <- ranked_df[1:9, ]

# Step 3: Get the names (GeneName) of the top 9 rows
top_9_gene_names <- top_9_rows$name

cytokines_list <- c(
  "C1qtnf7", "Il6st", "Irf6", "Ccl2", "Irf2bpl", "Socs5", "Tnfrsf12a", "Ifrd1", "Socs2", 
  "Nfil3", "Il34", "Clcf1", "Csf1"
)

p <- plotEmbedding(
  ArchRProj = projPAG, 
  colorBy = "GeneScoreMatrix", 
  name = cytokines_list[1:9], 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)
p$Tmem176a

p2 <- lapply(p, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

projPAG <- addImputeWeights(projPAG)

p3 <- plotEmbedding(
  ArchRProj = projPAG, 
  colorBy = "GeneScoreMatrix", 
  name = cytokines_list[1:9], 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(projPAG)
)

p3 <- lapply(p3, function(x){
  x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p3))
