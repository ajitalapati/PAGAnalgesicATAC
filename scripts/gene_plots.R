setwd("~/Documents/PAGAnalgesicATAC")
library(ArchR)
library(parallel)
addArchRThreads(threads = 18) 

projPAG <- loadArchRProject(path = "./PAG_ATAC_Directory4", force = FALSE, showLogo = TRUE)

markersGS <- getMarkerFeatures(
  ArchRProj = projPAG, 
  useMatrix = "GeneScoreMatrix", 
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Get markers for each of the clusters
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")
markerList$C6$name

tf_list <- c(
  "Egr1", "Egr2", "Elk4", "Ezh2", "Hdac10", "Irf3", "Irf6", "Mzf1", "Nfatc1", "Nfatc4", 
  "Nfkb2", "Notch1", "Notch4", "Pou6f1", "Rfx2", "Rreb1", "Tcf4", "Zfp36l2", "Creb3", "Grip1",
  "Sox11", "Zbtb7a", "Camta1", "Ebf4", "Kcnip3", "Ovol2", "Pou4f3", "Rxrg", "Crebl2", 
  "Egr1", "Atf3", "Crem", "Ddit3", "Sox11", "Bach1", "Stat5a", "Snai2", "Fosl2", "Arid5a"
)

markerList$C1

# Step 1: Rank the data frame by the column "Log2FC" in descending order
ranked_df <- as.data.frame(markerList$C1[order(-markerList$C1$Log2FC), ])

# Step 2: Get the top 9 rows from the ranked data frame
top_9_rows <- ranked_df[1:9, ]

# Step 3: Get the names (GeneName) of the top 9 rows
top_9_gene_names <- top_9_rows$name
top_9_gene_names

cytokines_list <- c(
  "C1qtnf7", "Il6st", "Irf6", "Ccl2", "Irf2bpl", "Socs5", "Tnfrsf12a", "Ifrd1", "Socs2", 
  "Nfil3", "Il34", "Clcf1", "Csf1"
)

p <- plotEmbedding(
  ArchRProj = projPAG, 
  colorBy = "GeneScoreMatrix", 
  name = top_9_gene_names[1:9], 
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