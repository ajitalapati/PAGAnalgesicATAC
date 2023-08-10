setwd("~/Documents/PAGAnalgesicATAC")
library(ArchR)
library(parallel)
addArchRThreads(threads = 22) 

library(BSgenome.Rnorvegicus.UCSC.rn7)
library(TxDb.Rnorvegicus.UCSC.rn7.refGene)
library(org.Rn.eg.db)

library(Cairo)

quartz(width=6, height=4)

projPAG <- loadArchRProject(path = "./PAG_ATAC_Directory4", force = FALSE, showLogo = FALSE)
pal_1 <- c("Saline_Veh"="#9983BD", "CFA_Veh"="#FCBF6E", "CFA_ApAP"="#90D5E4", "CFA_3DDA"="#F37B7D")
magma_pal <- c("Saline_Veh"="#4B2991", "CFA_Veh"="#ED5983", "CFA_ApAP"="#F89078", "CFA_3DDA"="#F1C18E")

gpcrGeneSymbols = c("Omg", "Chrm3", "Gpr149", "Gpr19")
ionGeneSymbols = c(
  "Scn8a",
  "Kcnma1",
  "Cacna1a",
  "Kcnh5",
  "Scn3a",
  "Kcng2",
  "Ano1",
  "Asic2",
  "Gabrb3",
  "Kcnq2",
  "Gabra1",
  "Htr3a",
  "Gabrg2",
  "Cacna1b",
  "Trpv2",
  "Scn9a",
  "Scn1a"
)
projPAG <- addCoAccessibility(
  ArchRProj = projPAG,
  reducedDims = "IterativeLSI"
)

for (gene in gpcrGeneSymbols){
  currGroup <- plotBrowserTrack(
    ArchRProj = projPAG,
    groupBy = "Sample",
    plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
    sizes = c(8, 1.2, 3, 2),
    geneSymbol = gene,
    upstream = 50000,
    downstream = 50000,
    loops = getCoAccessibility(projPAG),
    pal = magma_pal,
    baseSize = 7,
    ylim = c(0,1),
    scTileSize = 0.5,
    scCellsMax = 100,
    borderWidth = 0.4,
    tickWidth = 0.4,
    facetbaseSize = 7,
    title = gene,
    verbose = TRUE,
    logFile = createLogFile("plotBrowserTrack")
  )
  
  setwd("~/Documents/PAGAnalgesicATAC/Browser\ Tracks")
  
  grid::grid.newpage()
  grid::grid.draw(currGroup[[gene]])
  ggsave(paste(gene,".pdf",sep = ''), currGroup[[gene]])
}



for (gene in ionGeneSymbols){
  setwd("~/Documents/PAGAnalgesicATAC/UMAPs")
  currUMAP <- plotEmbedding(
    ArchRProj = projPAG, 
    colorBy = "GeneScoreMatrix", 
    name = gene, # change this to be a list of genes we want
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = NULL, # if you want imputation weights, replace with: getImputeWeights(projPAG)
    pal = paletteContinuous(set="purpleOrange"),
    rastr = TRUE
  )
  plotPDF(currUMAP, name = gene, addDOC = FALSE, width = 6, height = 6)
}

projPAG <- addArchRAnnotations(ArchRProj = projPAG, collection = "ATAC")

projPAG <- addDeviationsMatrix(
  ArchRProj = projPAG, 
  peakAnnotation = "ATAC",
  force = TRUE
)


motifs <- c("Sox2.HMG_261", "Sox15.HMG_260", "Sox10.HMG_259", "Sox6.HMG_264", "Sox9.HMG_265", "Sox4.HMG_263", "Sox3.HMG_262", 
               "SpiB.ETS_268", "Rfx2.HTH_242", "Rfx1.HTH_241", "RORgt.NR_246", "ELF5.ETS_67", "ETS1.ETS_76")
markerMotifs <- getFeatures(projPAG, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

deviationUMAP <- plotEmbedding(
  ArchRProj = projPAG, 
  colorBy = "MotifMatrix", 
  name = sort(markerMotifs), 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(projPAG)
)
deviationUMAP

for (motif in motifs){
  setwd("~/Documents/PAGAnalgesicATAC/DevUMAP")
  currUmap = paste("deviations:", motif, sep='')
  plotPDF(deviationUMAP[currUmap], name = motif, addDOC = FALSE, width = 6, height = 6)
}

deviationUMAP$"deviations:ELF5.ETS_67"
currUmap = paste("deviations:", "ELF5.ETS_67", sep='')
deviationUMAP[currUmap]

motifPositions <- getPositions(projPAG)
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
markerMotifs

seFoot <- getFootprints(
  ArchRProj = projPAG, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Sample"
)

plotFootprints(
  seFoot = seFoot,
  pal = magma_pal,
  ArchRProj = projPAG, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 5
)
