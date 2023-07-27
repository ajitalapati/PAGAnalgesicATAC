setwd("~/Documents/PAGAnalgesicATAC")
library(ArchR)
library(parallel)
addArchRThreads(threads = 22) 

library(BSgenome.Rnorvegicus.UCSC.rn7)
library(TxDb.Rnorvegicus.UCSC.rn7.refGene)
library(org.Rn.eg.db)

library(Cairo)

projPAG5 <- loadArchRProject(path = "./PAG_ATAC_Directory4", force = FALSE, showLogo = FALSE)

#15.1
projPAG5 <- addCoAccessibility(
  ArchRProj = projPAG5,
  reducedDims = "IterativeLSI"
)

cA <- getCoAccessibility(
  ArchRProj = projPAG5,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = FALSE
)
metadata(cA)[[1]]

p <- plotBrowserTrack(
  ArchRProj = projPAG5, 
  groupBy = "Sample", 
  geneSymbol = c("CFA_3DDA", "CFA_Veh", "Saline_Veh"), 
  upstream = 50000,
  downstream = 50000,
  loops = getCoAccessibility(projPAG4)
)
p

seGroupMotif <- getGroupSE(ArchRProj = projPAG5, useMatrix = "MotifMatrix", groupBy = "Sample")

seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs

corGSM_MM <- correlateMatrices(
  ArchRProj = projPAG5,
  useMatrix1 = "GeneScoreMatrix",
  useMatrix2 = "MotifMatrix",
  reducedDims = "IterativeLSI"
)

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
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )
p

p1 <- plotEmbedding(ArchRProj = projPAG5, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p1
