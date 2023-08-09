setwd("/work/aalapa/pag")
library(ArchR)
library(parallel)
addArchRThreads(threads = 18) 

library(BSgenome.Rnorvegicus.UCSC.rn7)
library(TxDb.Rnorvegicus.UCSC.rn7.refGene)
library(org.Rn.eg.db)

projPAG3 <- loadArchRProject(path = "./PAG_ATAC_Directory2", force = FALSE, showLogo = TRUE)
 
projPAG <- addGroupCoverages(ArchRProj = projPAG,
                  groupBy = "Clusters", force = TRUE)

pathToMacs2 <- "/home/aalapa/.local/share/r-miniconda/envs/PeakCalling_analysis/bin/macs2"
pathToMacs2 <- "/work/aalapa/"
pathToMacs2 <- findMacs2()

projPAG <- addReproduciblePeakSet(
  ArchRProj = projPAG, 
  groupBy = "Clusters", 
  genomeSize = 2647915728,
  pathToMacs2 = pathToMacs2
)

getPeakSet(projPAG)
projPAG <- addPeakMatrix(projPAG)

getAvailableMatrices(projPAG3)
table(projPAG3$Clusters)

saveArchRProject(ArchRProj = projPAG, outputDirectory = "PAG_ATAC_Directory4", load = FALSE)
