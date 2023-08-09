setwd("~/Documents/PAGAnalgesicATAC")
library(ArchR)
library(parallel)
addArchRThreads(threads = 18) 
set.seed(1) 

library(cowplot)

projPAG <- loadArchRProject(path = "./PAG_ATAC_Directory4", force = FALSE, showLogo = FALSE)

#2a __________________________________________
p <- plotFragmentSizes(
  ArchRProj = projPAG,
  groupBy = "Sample",
  chromSizes = getChromSizes(projPAG),
  maxSize = 750,
  pal = NULL,
  returnDF = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("plotFragmentSizes")
)

p3 <- plotEmbedding(ArchRProj = projPAG, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")  + theme(plot.margin = margin(0, 0, 5, 5) , labelMeans =FALSE)
p3

counter_hashmap <- table(projPAG2$Clusters)

# Convert the counter hashmap to a data frame
data_df <- data.frame(values = names(counter_hashmap), counts = as.vector(counter_hashmap))

# Sort the data frame by 'values' in alphanumerical order
data_df$values <- factor(data_df$values, levels = sort(unique(data_df$values)))

# Plot the histogram using ggplot2
count_hist <- ggplot(data_df, aes(x = counts, y = values, fill = values)) +
  geom_bar(stat = "identity", fill = "lightblue", color = "black", position = "dodge") +
  labs(title = "Number of Cells",
       x = "Frequency", y = "") +
  theme_minimal()

widths <- c(3, 1.5)

# Combine the plots side by side using cowplot
final_plot <- plot_grid(p1, count_hist + geom_text(aes(label = "")), ncol = 2, rel_widths = widths)

# Print the final plot
print(final_plot)
#2a end__________________________________________

#2c __________________________________________
markerGenes  <- c(
  "CD34",  #Early Progenitor
  "GATA1", #Erythroid
  "PAX5", "MS4A1", "MME", #B-Cell Trajectory
  "CD14", "MPO", #Monocytes
  "CD3D", "CD8A"#TCells
)

p <- plotEmbedding(
  ArchRProj = projPAG, 
  colorBy = "GeneScoreMatrix", 
  name = markerGenes, # change this to be a list of genes we want
  embedding = "UMAP",
  quantCut = c(0.01, 0.95),
  imputeWeights = NULL # if you want imputation weights, replace with: getImputeWeights(projPAG)
)

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
do.call(cowplot::plot_grid, c(list(ncol = 3),p2)) #adjust num of columns at ncol parameter
#2c end__________________________________________

#2d __________________________________________
markersPeaks <- getMarkerFeatures(
  ArchRProj = projPAG, 
  useMatrix = "PeakMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

heatmapPeaks <- markerHeatmap(
  seMarker = markersPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

#2d end __________________________________________

#2e __________________________________________

#if needed refer to: https://www.archrproject.com/bookdown/motif-deviations.html
plotVarDev <- getVarDeviations(projPAG, name = "MotifMatrix", plot = TRUE)

motifs <- c("GATA1", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5")
markerMotifs <- getFeatures(projPAG, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs

p <- plotEmbedding(
  ArchRProj = projPAG, 
  colorBy = "MotifMatrix", 
  name = sort(markerMotifs), 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(projPAG)
)

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

#2e end__________________________________________

markersPeaks <- getMarkerFeatures(
  ArchRProj = projPAG, 
  useMatrix = "PeakMatrix", 
  groupBy = "Sample",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerTest <- getMarkerFeatures(
  ArchRProj = projPAG, 
  useMatrix = "PeakMatrix",
  groupBy = "Sample",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "CFA_3DDA",
  bgdGroups = "CFA_ApAP"
)

pv <- plotMarkers(seMarker = markerTest, name = "CFA_3DDA", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
pv

motifsUp <- peakAnnoEnrichment(
  seMarker = markerTest,
  ArchRProj = projPAG,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

enrichMotifs <- peakAnnoEnrichment(
  seMarker = markersPeaks,
  ArchRProj = projPAG,
  peakAnnotation = "Motif",
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)

addPeakSet(
  ArchRProj = projPAG,
  peakSet = NULL,
  genomeAnnotation = getGenomeAnnotation(projPAG),
  force = FALSE
)






