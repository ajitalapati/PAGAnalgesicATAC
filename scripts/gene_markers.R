setwd("/work/aalapa/pag")
library(ArchR)
library(parallel)
addArchRThreads(threads = 18) 

projPAG <- loadArchRProject(path = "./PAG_ATAC_Directory", force = FALSE, showLogo = TRUE)

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

# access the clustering information
head(projPAG2$Clusters)

# get number of cells present in each cluster
table(projPAG2$Clusters)

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
projPAG2 <- addUMAP(
  ArchRProj = projPAG2, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

# plotting the embedding
p1 <- plotEmbedding(ArchRProj = projPAG2, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p1

p2 <- plotEmbedding(ArchRProj = projPAG2, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2

plotPDF(p1,p2, name = "umap-colored-sample-cluster.pdf", ArchRProj = projPAG, addDOC = FALSE, width = 5, height = 5)


# Now time to compute gene score matrices

devtools::install_github("immunogenomics/presto")

markersGS <- getMarkerFeatures(
  ArchRProj = projPAG2, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Get markers for each of the clusters
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C6

# To visualize all of the marker features
markerGenes  <- c(
  "CD34", #Early Progenitor
  "GATA1", #Erythroid
  "PAX5", "MS4A1", "EBF1", "MME", #B-Cell Trajectory
  "CD14", "CEBPB", "MPO", #Monocytes
  "IRF8", 
  "CD3D", "CD8A", "TBX21", "IL7R" #TCells
)

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = NULL, #can set to marker genes for labelling once defined
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

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


