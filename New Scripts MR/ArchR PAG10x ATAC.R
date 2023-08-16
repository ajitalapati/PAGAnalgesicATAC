if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

devtools::install_github("GreenleafLab/ArchR", ref="master", 
                         repos = BiocManager::repositories())

#notes from manual 
#Archr arrow files are comparable to Seurat objects which store data in levels
#as metadata withing the object

library(ArchR)
set.seed(1)

addArchRThreads(threads = 2)

#ArchR has a precompiled mm10 genome with blacklist regions and annotations so no need to compile one
#just kidding this is rat data so we will need to compilea custom rn6 or rn7 ArchR genome object guided by the tutorial

#first identify and install and load the relevant BSgenome object.
if (!requireNamespace("BSgenome.Rnorvegicus.UCSC.rn7", quietly = TRUE)){
  BiocManager::install("BSgenome.Rnorvegicus.UCSC.rn7")
}
library(BSgenome.Rnorvegicus.UCSC.rn7)

# create genome annotation
genomeAnnotation.rn7 <- createGenomeAnnotation(genome = BSgenome.Rnorvegicus.UCSC.rn7)

# install transcript and organism databased for annotaiton

if (!requireNamespace("TxDb.Rnorvegicus.UCSC.rn7.refGene", quietly = TRUE)){
  BiocManager::install("TxDb.Rnorvegicus.UCSC.rn7.refGene")
}
if (!requireNamespace("org.Rn.eg.db", quietly = TRUE)){
  BiocManager::install("org.Rn.eg.db")
}
library(TxDb.Rnorvegicus.UCSC.rn7.refGene)
library(org.Rn.eg.db)

geneAnnotation.rn7 <- createGeneAnnotation(TxDb = TxDb.Rnorvegicus.UCSC.rn7.refGene, OrgDb = org.Rn.eg.db)

geneAnnotation.rn7<-createGeneAnnotation(TSS = geneAnnotation.rn7$TSS, exons = geneAnnotation.rn7$exons, 
                                     genes = geneAnnotation.rn7$genes)

# load the fragment files before data aggregation with Cellranger

inputFiles <- c("CFA_Veh_atac_fragments.tsv.gz", 
                "CFA_ApAP_atac_fragments.tsv.gz",
                "CFA_3DDA_atac_fragments.tsv.gz",
                "Saline_Veh_atac_fragments.tsv.gz")

sampleNames <- c("CFA_Veh", "CFA_ApAP", "CFA_3DDA", "Saline_Veh") 

inputFiles

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = sampleNames,
  geneAnnotation = geneAnnotation.rn7,
  genomeAnnotation = genomeAnnotation.rn7,
  minTSS = 1, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

#view arrow file/files to check that it loaded a string of charachters
ArrowFiles

# Doublet inference using ArchR on Arrow Files
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1
)

# Create ArchR project using the loaded ArrowFiles
projPAG <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  geneAnnotation = geneAnnotation.rn7,
  genomeAnnotation = genomeAnnotation.rn7,
  outputDirectory ="PAG_ATAC_Directory",
  copyArrows = TRUE #This is recommended so that if you modify the Arrow files you have an original copy for later usage.
)

saveArchRProject(ArchRProj = projPAG, outputDirectory = "PAG_ATAC_Directory", load = FALSE)

projPAG <- readRDS("Save-ArchR-Project.rds")

# We save this as a new ArchRProject for the purposes of this stepwise tutorial but you can overwrite the ArchRProject instead
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

# To visualize all of the marker features
markerGenes  <- c(
  "CD34", #Early Progenitor
  "GATA1", #Erythroid
  "PAX5", "MS4A1", "EBF1", "MME", #B-Cell Trajectory
  "CD14", "CEBPB", "MPO", #Monocytes
  "IRF8", 
  "CD3D", "CD8A", "TBX21", "IL7R" #TCells
)

library(magick)
library(Cairo)

heatmapGS <- markerHeatmap(
  seMarker = markersGS, 
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25", 
  labelMarkers = markerGenes, #can set to marker genes for labelling once defined
  transpose = TRUE
)

ComplexHeatmap::draw(heatmapGS, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapGS, name = "GeneScores-Marker-Heatmap", width = 8, height = 6, ArchRProj = projPAG2, addDOC = FALSE)

# plotting marker genes in a UMAP embedding will change the marker gene vector to match our data set once we decide cell types and top marker genes
p <- plotEmbedding(
  ArchRProj = projPAG2, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, 
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL
)
# go through ArchR tutorial to do different plots of this type to define markers for each clouster later and perform imputation
# imputations using MAGIC smooths the embedding markeplot.spec.coherency and follow the rest of the sections in chapter 7

# skipping to chapter 10 calling peaks to see if RNA integration is necessary for this analysis not using macs2 instead using TileMatrix
# ArchR developed tile matrix for instances when Macs2 cannot be used but since it may be easier I want to check and see how it works

# this step is required before peaks can be called changed the group by since we didn't create a new cluster category from RNA integration
projPAG2 <- addGroupCoverages(ArchRProj = projPAG2, groupBy = "Clusters")

library(reticulate)
BiocManager::install("Herper")
library(Herper)

install_CondaTools(tools="macs2", env="PeakCalling_analysis", pathToMiniConda="~/Library/r-miniconda")
pathToMacs2 <- "/Users/madiganreid/Library/r-miniconda/envs/PeakCalling_analysis/bin/macs2"


projPAG2 <- addReproduciblePeakSet(
  ArchRProj = projPAG2, 
  groupBy = "Clusters",
  pathToMacs2 = pathToMacs2,
  genomeSize = 2647915728,
  method = "p"
)

getPeakSet(projPAG2)


saveArchRProject(ArchRProj = projPAG2, outputDirectory = "Save-ProjPAG2", load = FALSE)

projPAG2 <- addPeakMatrix(projPAG2)

getAvailableMatrices(projPAG2)

table(projPAG2$Clusters)

markersPeaks <- getMarkerFeatures(
  ArchRProj = projPAG2, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

draw(heatmapPeaks, heatmap_legend_side = "bot", annotation_legend_side = "bot")

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = projPAG2, addDOC = FALSE)

pma <- markerPlot(seMarker = markersPeaks, name = "C1", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "MA")
pma


pv <- markerPlot(seMarker = markersPeaks, name = "C1", cutOff = "FDR <= 0.1 & Log2FC >= 1", plotAs = "Volcano")
pv


plotPDF(pma, pv, name = "C1-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = projPAG2, addDOC = FALSE)


p <- plotBrowserTrack(
  ArchRProj = projPAG2, 
  groupBy = "Clusters", 
  geneSymbol = c("KCNIP3"),
  features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 1", returnGR = TRUE)["C1"],
  upstream = 50000,
  downstream = 50000
)

grid::grid.newpage()
grid::grid.draw(p$Kcnip3)

plotPDF(p, name = "Plot-Tracks-With-Features", width = 5, height = 5, ArchRProj = projPAG2, addDOC = FALSE)

markerTest <- getMarkerFeatures(
  ArchRProj = projPAG2, 
  useMatrix = "PeakMatrix",
  groupBy = "Clusters",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "C1",
  bgdGroups = "C2"
)

pma <- markerPlot(seMarker = markerTest, name = "C1", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "MA")
pma

pv <- markerPlot(seMarker = markerTest, name = "C1", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
pv

plotPDF(pma, pv, name = "C1-vs-C2-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = projPAG2, addDOC = FALSE)


devtools::install_github("GreenleafLab/chromVARmotifs")
BiocManager::install("JASPAR2020")

projPAG2 <- addMotifAnnotations(ArchRProj = projPAG2, motifSet = "encode", name = "Motif", species = getGenome(ArchRProj = projPAG2))

motifsUp <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = projPAG2,
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


motifsDo <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = projPAG2,
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

plotPDF(ggUp, ggDo, name = "C1-vs-C2-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = projPAG2, addDOC = FALSE)

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = projPAG2,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)

ComplexHeatmap::draw(heatmapEM, heatmap_legend_side = "bot", annotation_legend_side = "bot")


plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-Heatmap", width = 8, height = 6, ArchRProj = projPAG2, addDOC = FALSE)


if("Motif" %ni% names(projPAG2@peakAnnotation)){
  projPAG2 <- addMotifAnnotations(ArchRProj = projPAG2, motifSet = "encode", name = "Motif", species = getGenome(ArchRProj = projPAG2))
}

projPAG2 <- addBgdPeaks(projPAG2)


projPAG2 <- addDeviationsMatrix(
  ArchRProj = projPAG2, 
  peakAnnotation = "Motif",
  force = TRUE
)
