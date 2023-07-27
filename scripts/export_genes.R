setwd("~/Documents/PAGAnalgesicATAC")
library(ArchR)
library(parallel)
addArchRThreads(threads = 22) 

library(BSgenome.Rnorvegicus.UCSC.rn7)
library(TxDb.Rnorvegicus.UCSC.rn7.refGene)
library(org.Rn.eg.db)

library(Cairo)

projPAG5 <- loadArchRProject(path = "./PAG_ATAC_Directory4", force = FALSE, showLogo = FALSE)

#section 11.1 in manual
markerTest <- getMarkerFeatures(
  ArchRProj = projPAG5, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "Sample",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  bgdGroups = "CFA_ApAP",
  useGroups = "CFA_3DDA",
)

t <- head(getGenes(projPAG5), -3)

df <- data.frame(
                symbol = t,
                 Log2FC = assays(markerTest)$Log2FC, 
                 Pval = assays(markerTest)$Pval, 
                 FDR = assays(markerTest)$FDR)

colnames(df) <- c("chromosome", "start", "end", "width","strand","gene_id","symbol", "Log2FC", "Pval", "FDR")

write.csv(df, "~/Documents/PAGAnalgesicATAC/brian/CFA_ApAP2CFA_3DDA.csv", row.names=FALSE)
