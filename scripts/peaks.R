setwd("/work/aalapa/pag")
library(ArchR)
library(parallel)
addArchRThreads(threads = 18) 

library(BSgenome.Rnorvegicus.UCSC.rn7)
library(TxDb.Rnorvegicus.UCSC.rn7.refGene)
library(org.Rn.eg.db)

projPAG3 <- loadArchRProject(path = "./PAG_ATAC_Directory2", force = FALSE, showLogo = TRUE)
 
projPAG3 <- addGroupCoverages(ArchRProj = projPAG3,
                  groupBy = "Clusters")

pathToMacs2 <- "/home/aalapa/.local/share/r-miniconda/envs/PeakCalling_analysis/bin/macs2"
pathToMacs2 <- "/work/aalapa/"

projPAG3 <- addReproduciblePeakSet(
  ArchRProj = projPAG3, 
  groupBy = "Clusters", 
  genomeSize = 2647915728,
  pathToMacs2 = pathToMacs2
)

getPeakSet(projPAG3)
projPAG3 <- addPeakMatrix(projPAG3)

getAvailableMatrices(projPAG3)
table(projPAG3$Clusters)

saveArchRProject(ArchRProj = projPAG3, outputDirectory = "PAG_ATAC_Directory3", load = FALSE)