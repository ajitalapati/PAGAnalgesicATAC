### Brian Analysis ###

GenesSaline_Veh <- getMarkerFeatures(
  ArchRProj = projPAG5, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "Sample",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  bgdGroups = c("Saline_Veh"),
  useGroups = c("CFA_Veh", "CFA_ApAP", "CFA_3DDA"))

PeaksSaline_Veh <- getMarkerFeatures(
  ArchRProj = projPAG5, 
  useMatrix = "PeakMatrix",
  groupBy = "Sample",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  bgdGroups = c("Saline_Veh"),
  useGroups = c("CFA_Veh", "CFA_ApAP", "CFA_3DDA"))

GenesCFA_Veh <- getMarkerFeatures(
  ArchRProj = projPAG5, 
  useMatrix = "GeneScoreMatrix",
  groupBy = "Sample",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  bgdGroups = c("CFA_Veh"),
  useGroups = c("Saline_Veh", "CFA_ApAP", "CFA_3DDA"))

PeaksCFA_Veh <- getMarkerFeatures(
  ArchRProj = projPAG5, 
  useMatrix = "PeakMatrix",
  groupBy = "Sample",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  bgdGroups = c("CFA_Veh"),
  useGroups = c("Saline_Veh", "CFA_ApAP", "CFA_3DDA"))

GenesHeatmapSaline_Veh <- plotMarkerHeatmap(GenesSaline_Veh, plotLog2FC = TRUE, clusterCols = FALSE, transpose = TRUE)
GenesHeatmapCFA_Veh <- plotMarkerHeatmap(GenesCFA_Veh, plotLog2FC = TRUE, clusterCols = FALSE,transpose = TRUE, returnMatrix = TRUE)
PeaksHeatmapSaline_Veh <- plotMarkerHeatmap(PeaksSaline_Veh, plotLog2FC = TRUE, clusterCols = FALSE,transpose = TRUE)
PeaksHeatmapCFA_Veh <- plotMarkerHeatmap(PeaksCFA_Veh, plotLog2FC = TRUE, clusterCols = FALSE,transpose = TRUE)


GenesHeatmapSaline_Veh
PeaksHeatmapSaline_Veh
GenesHeatmapCFA_Veh
PeaksHeatmapCFA_Veh


write.csv(GenesHeatmapCFA_Veh, file = "MarkerGenes.csv")

### Browsertrack Plots for DAR Genes ###

# Saline Veh vs CFA Veh DAR

Acap3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Acap3", useMatrix = "GeneScoreMatrix")
Acap3

plotPDF(Acap3, name = "Acap3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Eif3c <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Eif3c", useMatrix = "GeneScoreMatrix")
Eif3c

plotPDF(Eif3c, name = "Eif3c browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Map2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Map2", useMatrix = "GeneScoreMatrix")
Map2

plotPDF(Map2, name = "Map2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Mpp4 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Mpp4", useMatrix = "GeneScoreMatrix")
Mpp4

plotPDF(Mpp4, name = "Mpp4 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Xirp2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Xirp2", useMatrix = "GeneScoreMatrix")
Xirp2

plotPDF(Xirp2, name = "Xirp2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ncoa2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ncoa2", useMatrix = "GeneScoreMatrix")
Ncoa2

plotPDF(Ncoa2, name = "Ncoa2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ube2z <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ube2z", useMatrix = "GeneScoreMatrix")
Ube2z

plotPDF(Ube2z, name = "Ube2z browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Car1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Car1", useMatrix = "GeneScoreMatrix")
Car1

plotPDF(Car1, name = "Car1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Vom2r18 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Vom2r18", useMatrix = "GeneScoreMatrix")
Vom2r18

plotPDF(Vom2r18, name = "Vom2r18 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Or2ak6b <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Or2ak6b", useMatrix = "GeneScoreMatrix")
Or2ak6b

plotPDF(Or2ak6b, name = "Or2ak6b browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Hsfy2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Hsfy2", useMatrix = "GeneScoreMatrix")
Hsfy2

plotPDF(Hsfy2, name = "Hsfy2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Pag1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Pag1", useMatrix = "GeneScoreMatrix")
Pag1

plotPDF(Pag1, name = "Pag1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Lif <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Lif", useMatrix = "GeneScoreMatrix")
Lif

plotPDF(Lif, name = "Lif browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Mir3562 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Mir3562", useMatrix = "GeneScoreMatrix")
Mir3562

plotPDF(Mir3562, name = "Mir3562 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ube2j2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ube2j2", useMatrix = "GeneScoreMatrix")
Ube2j2

plotPDF(Ube2j2, name = "Ube2j2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ambn <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ambn", useMatrix = "GeneScoreMatrix")
Ambn

plotPDF(Ambn, name = "Ambn browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Astn2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Astn2", useMatrix = "GeneScoreMatrix")
Astn2

plotPDF(Astn2, name = "Astn2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Rps4y2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Rps4y2", useMatrix = "GeneScoreMatrix")
Rps4y2

plotPDF(Rps4y2, name = "Rps4y2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ces2j <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ces2j", useMatrix = "GeneScoreMatrix")
Ces2j

plotPDF(Ces2j, name = "Ces2j browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Mfsd5 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Mfsd5", useMatrix = "GeneScoreMatrix")
Mfsd5

plotPDF(Mfsd5, name = "Mfsd5 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

# CFA Veh vs CFA ApAP DAR

Fam111a <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Fam111a", useMatrix = "GeneScoreMatrix")
Fam111a

plotPDF(Fam111a, name = "Fam111a browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Or2j3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Or2j3", useMatrix = "GeneScoreMatrix")
Or2j3

plotPDF(Or2j3, name = "Or2j3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Smok2a <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Smok2a", useMatrix = "GeneScoreMatrix")
Smok2a

plotPDF(Smok2a, name = "Smok2a browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Nfia <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Nfia", useMatrix = "GeneScoreMatrix")
Nfia

plotPDF(Nfia, name = "Nfia browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Trmt1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Trmt1", useMatrix = "GeneScoreMatrix")
Trmt1

plotPDF(Trmt1, name = "Trmt1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Or2t29 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Or2t29", useMatrix = "GeneScoreMatrix")
Or2t29

plotPDF(Or2t29, name = "Or2t29 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Apof <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Apof", useMatrix = "GeneScoreMatrix")
Apof

plotPDF(Apof, name = "Apof browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Phkb <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Phkb", useMatrix = "GeneScoreMatrix")
Phkb

plotPDF(Phkb, name = "Phkb browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Mir21 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Mir21", useMatrix = "GeneScoreMatrix")
Mir21

plotPDF(Mir21, name = "Mir21 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Cyp4f18 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Cyp4f18", useMatrix = "GeneScoreMatrix")
Cyp4f18

plotPDF(Cyp4f18, name = "Cyp4f18 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Nell1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Nell1", useMatrix = "GeneScoreMatrix")
Nell1

plotPDF(Nell1, name = "Nell1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Trim69 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Trim69", useMatrix = "GeneScoreMatrix")
Trim69

plotPDF(Trim69, name = "Trim69 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Vom2r67 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Vom2r67", useMatrix = "GeneScoreMatrix")
Vom2r67

plotPDF(Vom2r67, name = "Vom2r67 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Slc6a2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Slc6a2", useMatrix = "GeneScoreMatrix")
Slc6a2

plotPDF(Slc6a2, name = "Slc6a2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Fam151b <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Fam151b", useMatrix = "GeneScoreMatrix")
Fam151b

plotPDF(Fam151b, name = "Fam151b browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Retnla <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Retnla", useMatrix = "GeneScoreMatrix")
Retnla

plotPDF(Retnla, name = "Retnla browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Slx4ip <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Slx4ip", useMatrix = "GeneScoreMatrix")
Slx4ip

plotPDF(Slx4ip, name = "Slx4ip browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Grap2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Grap2", useMatrix = "GeneScoreMatrix")
Grap2

plotPDF(Grap2, name = "Grap2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Rbis <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Rbis", useMatrix = "GeneScoreMatrix")
Rbis

plotPDF(Rbis, name = "Rbis browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

# CFA Veh vs CFA 3DDA

Lrp10 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Lrp10", useMatrix = "GeneScoreMatrix")
Lrp10

plotPDF(Lrp10, name = "Lrp10 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Rem2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Rem2", useMatrix = "GeneScoreMatrix")
Rem2

plotPDF(Rem2, name = "Rem2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Mettl3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Mettl3", useMatrix = "GeneScoreMatrix")
Mettl3

plotPDF(Mettl3, name = "Mettl3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Or10g1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Or10g1", useMatrix = "GeneScoreMatrix")
Or10g1

plotPDF(Or10g1, name = "Or10g1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Or4c125 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Or4c125", useMatrix = "GeneScoreMatrix")
Or4c125

plotPDF(Or4c125, name = "Or4c125 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Acap3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Acap3", useMatrix = "GeneScoreMatrix")
Acap3

plotPDF(Acap3, name = "Acap3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Tox4 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Tox4", useMatrix = "GeneScoreMatrix")
Tox4

plotPDF(Tox4, name = "Tox4 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Gsk3b <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Gsk3b", useMatrix = "GeneScoreMatrix")
Gsk3b

plotPDF(Gsk3b, name = "Gsk3b browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Tnp2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Tnp2", useMatrix = "GeneScoreMatrix")
Tnp2

plotPDF(Tnp2, name = "Tnp2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Or4a15 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Or4a15", useMatrix = "GeneScoreMatrix")
Or4a15

plotPDF(Or4a15, name = "Or4a15 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Or4e2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Or4e2", useMatrix = "GeneScoreMatrix")
Or4e2

plotPDF(Or4e2, name = "Or4e2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Slc24a4 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Slc24a4", useMatrix = "GeneScoreMatrix")
Slc24a4

plotPDF(Slc24a4, name = "Slc24a4 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Prmt5 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Prmt5", useMatrix = "GeneScoreMatrix")
Prmt5

plotPDF(Prmt5, name = "Prmt5 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Haus4 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Haus4", useMatrix = "GeneScoreMatrix")
Haus4

plotPDF(Haus4, name = "Haus4 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Akna <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Akna", useMatrix = "GeneScoreMatrix")
Akna

plotPDF(Akna, name = "Akna browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Tie1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Tie1", useMatrix = "GeneScoreMatrix")
Tie1

plotPDF(Tie1, name = "Tie1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Inpp4b <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Inpp4b", useMatrix = "GeneScoreMatrix")
Inpp4b

plotPDF(Inpp4b, name = "Inpp4b browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Or52n2b <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Or52n2b", useMatrix = "GeneScoreMatrix")
Or52n2b

plotPDF(Or52n2b, name = "Or52n2b browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

# CFA ApAP vs CFA 3DDA

Fam111a <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Fam111a", useMatrix = "GeneScoreMatrix")
Fam111a

plotPDF(Fam111a, name = "Fam111a browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Car2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Car2", useMatrix = "GeneScoreMatrix")
Car2

plotPDF(Car2, name = "Car2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Cd2ap <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Cd2ap", useMatrix = "GeneScoreMatrix")
Cd2ap

plotPDF(Cd2ap, name = "Cd2ap browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Or2j3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Or2j3", useMatrix = "GeneScoreMatrix")
Or2j3

plotPDF(Or2j3, name = "Or2j3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Scml4 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Scml4", useMatrix = "GeneScoreMatrix")
Scml4

plotPDF(Scml4, name = "Scml4 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ncoa2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ncoa2", useMatrix = "GeneScoreMatrix")
Ncoa2

plotPDF(Ncoa2, name = "Ncoa2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Sox14 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Sox14", useMatrix = "GeneScoreMatrix")
Sox14

plotPDF(Sox14, name = "Sox14 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Eif3c <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Eif3c", useMatrix = "GeneScoreMatrix")
Eif3c

plotPDF(Eif3c, name = "Eif3c browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kdsr <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kdsr", useMatrix = "GeneScoreMatrix")
Kdsr

plotPDF(Kdsr, name = "Kdsr browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ly6e <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ly6e", useMatrix = "GeneScoreMatrix")
Ly6e

plotPDF(Ly6e, name = "Ly6e browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Tmcc3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Tmcc3", useMatrix = "GeneScoreMatrix")
Tmcc3

plotPDF(Tmcc3, name = "Tmcc3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Clip1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Clip1", useMatrix = "GeneScoreMatrix")
Clip1

plotPDF(Clip1, name = "Clip1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Oas1b <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Oas1b", useMatrix = "GeneScoreMatrix")
Oas1b

plotPDF(Oas1b, name = "Oas1b browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Herpud1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Herpud1", useMatrix = "GeneScoreMatrix")
Herpud1

plotPDF(Herpud1, name = "Herpud1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Calr3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Calr3", useMatrix = "GeneScoreMatrix")
Calr3

plotPDF(Calr3, name = "Calr3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Liph <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Liph", useMatrix = "GeneScoreMatrix")
Liph

plotPDF(Liph, name = "Liph browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Tob2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Tob2", useMatrix = "GeneScoreMatrix")
Tob2

plotPDF(Tob2, name = "Tob2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Plin2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Plin2", useMatrix = "GeneScoreMatrix")
Plin2

plotPDF(Plin2, name = "Plin2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)


### PAIN GENES OF INTEREST ###
#Transcription Factors (TF)
#1. Egr1,Egr2,Elk4,Ezh2,Hdac10,Irf3,Irf6,Mzf1,Nfatc1,Nfatc4,Nfkb2,Notch1,Notch4,Pou6f1,Rfx2,Rreb1,Tcf4,Zfp36l2,Creb3,Grip1,

#2. Sox11,Zbtb7a,Foxo6,Camta1,Ebf4,Kcnip3,Ovol2,Pou4f3,Rxrg,Crebl2,Egr1,Atf3,Crem,Ddit3,Sox11,Bach1,Stat5a,Snai2,Fosl2,Arid5a,Sox5a

Egr1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Egr1", useMatrix = "GeneScoreMatrix")
Egr1

plotPDF(Egr1, name = "Egr1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Egr2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Egr2", useMatrix = "GeneScoreMatrix")
Egr2

plotPDF(Egr2, name = "Egr2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Elk4 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Elk4", useMatrix = "GeneScoreMatrix")
Elk4

plotPDF(Elk4, name = "Elk4 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ezh2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ezh2", useMatrix = "GeneScoreMatrix")
Ezh2

plotPDF(Ezh2, name = "Ezh2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Hdac10 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Hdac10", useMatrix = "GeneScoreMatrix")
Hdac10

plotPDF(Hdac10, name = "Hdac10 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Irf3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Irf3", useMatrix = "GeneScoreMatrix")
Irf3

plotPDF(Irf3, name = "Irf3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Irf6 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Irf6", useMatrix = "GeneScoreMatrix")
Irf6

plotPDF(Irf6, name = "Irf6 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Mzf1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Mzf1", useMatrix = "GeneScoreMatrix")
Mzf1

plotPDF(Mzf1, name = "Mzf1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)


Nfatc4 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Nfatc4", useMatrix = "GeneScoreMatrix")
Nfatc4

plotPDF(Nfatc4, name = "Nfatc4 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Nfatc1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Nfatc1", useMatrix = "GeneScoreMatrix")
Nfatc1

plotPDF(Nfatc1, name = "Nfatc1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Nfkb2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Nfkb2", useMatrix = "GeneScoreMatrix")
Nfkb2

plotPDF(Nfkb2, name = "Nfkb2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Notch1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Notch1", useMatrix = "GeneScoreMatrix")
Notch1

plotPDF(Notch1, name = "Notch1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Notch4 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Notch4", useMatrix = "GeneScoreMatrix")
Notch4

plotPDF(Notch4, name = "Notch4 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Pou6f1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Pou6f1", useMatrix = "GeneScoreMatrix")
Pou6f1

plotPDF(Pou6f1, name = "Pou6f1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Rfx2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Rfx2", useMatrix = "GeneScoreMatrix")
Rfx2

plotPDF(Rfx2, name = "Rfx2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Rreb1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Rreb1", useMatrix = "GeneScoreMatrix")
Rreb1

plotPDF(Rreb1, name = "Rreb1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Tcf4 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Tcf4", useMatrix = "GeneScoreMatrix")
Tcf4

plotPDF(Tcf4, name = "Tcf4 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Zfp36l2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Zfp36l2", useMatrix = "GeneScoreMatrix")
Zfp36l2

plotPDF(Zfp36l2, name = "Zfp36l2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Creb3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Creb3", useMatrix = "GeneScoreMatrix")
Creb3

plotPDF(Creb3, name = "Creb3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Grip1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Grip1", useMatrix = "GeneScoreMatrix")
Grip1

plotPDF(Grip1, name = "Grip1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Sox11 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Sox11", useMatrix = "GeneScoreMatrix")
Sox11

plotPDF(Sox11, name = "Sox11 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Zbtb7a <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Zbtb7a", useMatrix = "GeneScoreMatrix")
Zbtb7a

plotPDF(Zbtb7a, name = "Zbtb7a browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Camta1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Camta1", useMatrix = "GeneScoreMatrix")
Camta1

plotPDF(Camta1, name = "Camta1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ebf4 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ebf4", useMatrix = "GeneScoreMatrix")
Ebf4

plotPDF(Ebf4, name = "Ebf4 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcnip3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcnip3", useMatrix = "GeneScoreMatrix")
Kcnip3

plotPDF(Kcnip3, name = "Kcnip3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ovol2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ovol2", useMatrix = "GeneScoreMatrix")
Ovol2

plotPDF(Ovol2, name = "Ovol2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Pou4f3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Pou4f3", useMatrix = "GeneScoreMatrix")
Pou4f3

plotPDF(Pou4f3, name = "Pou4f3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Rxrg <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Rxrg", useMatrix = "GeneScoreMatrix")
Rxrg

plotPDF(Rxrg, name = "Rxrg browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Crebl2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Crebl2", useMatrix = "GeneScoreMatrix")
Crebl2

plotPDF(Crebl2, name = "Crebl2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Egr1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Egr1", useMatrix = "GeneScoreMatrix")
Egr1

plotPDF(Egr1, name = "Egr1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Atf3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Atf3", useMatrix = "GeneScoreMatrix")
Atf3

plotPDF(Atf3, name = "Atf3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Crem <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Crem", useMatrix = "GeneScoreMatrix")
Crem

plotPDF(Crem, name = "Crem browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ddit3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ddit3", useMatrix = "GeneScoreMatrix")
Ddit3

plotPDF(Ddit3, name = "Ddit3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Bach1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Bach1", useMatrix = "GeneScoreMatrix")
Bach1

plotPDF(Bach1, name = "Bach1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Stat5a <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Stat5a", useMatrix = "GeneScoreMatrix")
Stat5a

plotPDF(Stat5a, name = "Stat5a browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Snai2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Snai2", useMatrix = "GeneScoreMatrix")
Snai2

plotPDF(Snai2, name = "Snai2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Fosl2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Fosl2", useMatrix = "GeneScoreMatrix")
Fosl2

plotPDF(Fosl2, name = "Fosl2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Arid5a <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Arid5a", useMatrix = "GeneScoreMatrix")
Arid5a

plotPDF(Arid5a, name = "Arid5a browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

#GPCRs
#1. Tas1r1,Vom2r44,Gpr157,Gpr81,Npr3,C5ar1,Lpar1,Ptger1,Fzd2,Fzd6,Fzd7,Lpar6,Smo,Gpr176,Ptgfr,Gabbr2,Adgrl2,Gpr85,Gpr68,Omg,Mrgprx1,Gpr149,

#2. Mrgpre,Adgrg2,Htr1b,Grm7,Chrm3,F2rl2,Cnr1,Prokr2,Oprl1,Bdkrb2,Oprm1,Htr2a,Ptger4,Chrm2,Galr1,Gpr158,Nyp2r,Htr2b,Olr1746,Cckar,Gpr153,

#3. Gpr158,Adra1b,Ackr1,Adgrd1,Cckbr,Gpr151,Qrfpr,Gpr19,Galr1,Ntsr2,F2rl2,Ptger1,Gpr68,Adra2c,Gpr149,Gpr45,Htr1a,Oprk1,Mrgpr,Mrgprx3,Rgs9, 

Tas1r1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Tas1r1", useMatrix = "GeneScoreMatrix")
Tas1r1

plotPDF(Tas1r1, name = "Tas1r1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Vom2r44 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Vom2r44", useMatrix = "GeneScoreMatrix")
Vom2r44

plotPDF(Vom2r44, name = "Vom2r44 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Gpr157 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Gpr157", useMatrix = "GeneScoreMatrix")
Gpr157

plotPDF(Gpr157, name = "Gpr157 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Npr3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Npr3", useMatrix = "GeneScoreMatrix")
Npr3

plotPDF(Npr3, name = "Npr3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

C5ar1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "C5ar1", useMatrix = "GeneScoreMatrix")
C5ar1

plotPDF(C5ar1, name = "C5ar1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Lpar1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Lpar1", useMatrix = "GeneScoreMatrix")
Lpar1

plotPDF(Lpar1, name = "Lpar1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ptger1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ptger1", useMatrix = "GeneScoreMatrix")
Ptger1

plotPDF(Ptger1, name = "Ptger1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Fzd2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Fzd2", useMatrix = "GeneScoreMatrix")
Fzd2

plotPDF(Fzd2, name = "Fzd2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Fzd6 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Fzd6", useMatrix = "GeneScoreMatrix")
Fzd6

plotPDF(Fzd6, name = "Fzd6 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Fzd7 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Fzd7", useMatrix = "GeneScoreMatrix")
Fzd7

plotPDF(Fzd7, name = "Fzd7 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Lpar6 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Lpar6", useMatrix = "GeneScoreMatrix")
Lpar6

plotPDF(Lpar6, name = "Lpar6 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Smo <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Smo", useMatrix = "GeneScoreMatrix")
Smo

plotPDF(Smo, name = "Smo browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Gpr176 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Gpr176", useMatrix = "GeneScoreMatrix")
Gpr176

plotPDF(Gpr176, name = "Gpr176 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ptgfr <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ptgfr", useMatrix = "GeneScoreMatrix")
Ptgfr

plotPDF(Ptgfr, name = "Ptgfr browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Gabbr2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Gabbr2", useMatrix = "GeneScoreMatrix")
Gabbr2

plotPDF(Gabbr2, name = "Gabbr2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Adgrl2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Adgrl2", useMatrix = "GeneScoreMatrix")
Adgrl2

plotPDF(Adgrl2, name = "Adgrl2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Gpr85 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Gpr85", useMatrix = "GeneScoreMatrix")
Gpr85

plotPDF(Gpr85, name = "Gpr85 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Gpr68 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Gpr68", useMatrix = "GeneScoreMatrix")
Gpr68

plotPDF(Gpr68, name = "Gpr68 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Omg <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Omg", useMatrix = "GeneScoreMatrix")
Omg

plotPDF(Omg, name = "Omg browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Mrgprx1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Mrgprx1", useMatrix = "GeneScoreMatrix")
Mrgprx1

plotPDF(Mrgprx1, name = "Mrgprx1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Gpr149 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Gpr149", useMatrix = "GeneScoreMatrix")
Gpr149

plotPDF(Gpr149, name = "Gpr149 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Mrgpre <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Mrgpre", useMatrix = "GeneScoreMatrix")
Mrgpre

plotPDF(Mrgpre, name = "Mrgpre browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Adgrg2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Adgrg2", useMatrix = "GeneScoreMatrix")
Adgrg2

plotPDF(Adgrg2, name = "Adgrg2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Htr1b <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Htr1b", useMatrix = "GeneScoreMatrix")
Htr1b

plotPDF(Htr1b, name = "Htr1b browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Grm7 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Grm7", useMatrix = "GeneScoreMatrix")
Grm7

plotPDF(Grm7, name = "Grm7 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Chrm3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Chrm3", useMatrix = "GeneScoreMatrix")
Chrm3

plotPDF(Chrm3, name = "Chrm3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

F2rl2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "F2rl2", useMatrix = "GeneScoreMatrix")
F2rl2

plotPDF(F2rl2, name = "F2rl2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Cnr1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Cnr1", useMatrix = "GeneScoreMatrix")
Cnr1

plotPDF(Cnr1, name = "Cnr1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Prokr2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Prokr2", useMatrix = "GeneScoreMatrix")
Prokr2

plotPDF(Prokr2, name = "Prokr2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Oprl1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Oprl1", useMatrix = "GeneScoreMatrix")
Oprl1

plotPDF(Oprl1, name = "Oprl1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Bdkrb2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Bdkrb2", useMatrix = "GeneScoreMatrix")
Bdkrb2

plotPDF(Bdkrb2, name = "Bdkrb2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Oprm1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Oprm1", useMatrix = "GeneScoreMatrix")
Oprm1

plotPDF(Oprm1, name = "Oprm1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Htr2a <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Htr2a", useMatrix = "GeneScoreMatrix")
Htr2a

plotPDF(Htr2a, name = "Htr2a browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ptger4 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ptger4", useMatrix = "GeneScoreMatrix")
Ptger4

plotPDF(Ptger4, name = "Ptger4 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Chrm2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Chrm2", useMatrix = "GeneScoreMatrix")
Chrm2

plotPDF(Chrm2, name = "Chrm2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Galr1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Galr1", useMatrix = "GeneScoreMatrix")
Galr1

plotPDF(Galr1, name = "Galr1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Gpr158 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Gpr158", useMatrix = "GeneScoreMatrix")
Gpr158

plotPDF(Gpr158, name = "Gpr158 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Htr2b <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Htr2b", useMatrix = "GeneScoreMatrix")
Htr2b

plotPDF(Htr2b, name = "Htr2b browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Cckar <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Cckar", useMatrix = "GeneScoreMatrix")
Cckar

plotPDF(Cckar, name = "Cckar browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Gpr153 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Gpr153", useMatrix = "GeneScoreMatrix")
Gpr153

plotPDF(Gpr153, name = "Gpr153 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Adra1b <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Adra1b", useMatrix = "GeneScoreMatrix")
Adra1b

plotPDF(Adra1b, name = "Adra1b browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Cckbr <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Cckbr", useMatrix = "GeneScoreMatrix")
Cckbr

plotPDF(Cckbr, name = "Cckbr browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Gpr151 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Gpr151", useMatrix = "GeneScoreMatrix")
Gpr151

plotPDF(Gpr151, name = "Gpr151 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Qrfpr <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Qrfpr", useMatrix = "GeneScoreMatrix")
Qrfpr

plotPDF(Qrfpr, name = "Qrfpr browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Gpr19 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Gpr19", useMatrix = "GeneScoreMatrix")
Gpr19

plotPDF(Gpr19, name = "Gpr19 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ntsr2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ntsr2", useMatrix = "GeneScoreMatrix")
Ntsr2

plotPDF(Ntsr2, name = "Ntsr2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

F2rl2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "F2rl2", useMatrix = "GeneScoreMatrix")
F2rl2

plotPDF(F2rl2, name = "F2rl2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Adra2c <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Adra2c", useMatrix = "GeneScoreMatrix")
Adra2c

plotPDF(Adra2c, name = "Adra2c browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Gpr45 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Gpr45", useMatrix = "GeneScoreMatrix")
Gpr45

plotPDF(Gpr45, name = "Gpr45 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Htr1a <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Htr1a", useMatrix = "GeneScoreMatrix")
Htr1a

plotPDF(Htr1a, name = "Htr1a browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Oprk1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Oprk1", useMatrix = "GeneScoreMatrix")
Oprk1

plotPDF(Oprk1, name = "Oprk1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Mrgprx3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Mrgprx3", useMatrix = "GeneScoreMatrix")
Mrgprx3

plotPDF(Mrgprx3, name = "Mrgprx3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Rgs9 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Rgs9", useMatrix = "GeneScoreMatrix")
Rgs9

plotPDF(Rgs9, name = "Rgs9 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

#Growth Factors and receptors
#1. Fgf3,Mydgf,Tgfb1,Fgf11,Igf2r,Fgfr1,Tgfa,Hgf,Pdgfa,Fgf13,Fgf1,Vegfa,Ngfr,Vegfb,Igf1r,Hrh1

Fgf3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Fgf3", useMatrix = "GeneScoreMatrix")
Fgf3

plotPDF(Fgf3, name = "Fgf3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Tgfb1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Tgfb1", useMatrix = "GeneScoreMatrix")
Tgfb1

plotPDF(Tgfb1, name = "Tgfb1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Fgf11 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Fgf11", useMatrix = "GeneScoreMatrix")
Fgf11

plotPDF(Fgf11, name = "Fgf11 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Igf2r <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Igf2r", useMatrix = "GeneScoreMatrix")
Igf2r

plotPDF(Igf2r, name = "Igf2r browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Fgfr1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Fgfr1", useMatrix = "GeneScoreMatrix")
Fgfr1

plotPDF(Fgfr1, name = "Fgfr1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Tgfa <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Tgfa", useMatrix = "GeneScoreMatrix")
Tgfa

plotPDF(Tgfa, name = "Tgfa browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Hgf <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Hgf", useMatrix = "GeneScoreMatrix")
Hgf

plotPDF(Hgf, name = "Hgf browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Pdgfa <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Pdgfa", useMatrix = "GeneScoreMatrix")
Pdgfa

plotPDF(Pdgfa, name = "Pdgfa browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Fgf13 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Fgf13", useMatrix = "GeneScoreMatrix")
Fgf13

plotPDF(Fgf13, name = "Fgf13 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Fgf1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Fgf1", useMatrix = "GeneScoreMatrix")
Fgf1

plotPDF(Fgf1, name = "Fgf1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Vegfa <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Vegfa", useMatrix = "GeneScoreMatrix")
Vegfa

plotPDF(Vegfa, name = "Vegfa browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ngfr <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ngfr", useMatrix = "GeneScoreMatrix")
Ngfr

plotPDF(Ngfr, name = "Ngfr browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Vegfb <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Vegfb", useMatrix = "GeneScoreMatrix")
Vegfb

plotPDF(Vegfb, name = "Vegfb browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Igf1r <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Igf1r", useMatrix = "GeneScoreMatrix")
Igf1r

plotPDF(Igf1r, name = "Igf1r browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Hrh1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Hrh1", useMatrix = "GeneScoreMatrix")
Hrh1

plotPDF(Hrh1, name = "Hrh1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)


#Cytokines and related molecules
#2. C1qtnf7,Il6st,Irf2bp2,Irf6,Ccl2,Irf2bpl,Socs5,Tnfrsf12a,Ifrd1,Socs2,Nfil3,Il34,Fam19a4,Clcf1,Ifitm10,Csf1,Ackr1

C1qtnf7 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "C1qtnf7", useMatrix = "GeneScoreMatrix")
C1qtnf7

plotPDF(C1qtnf7, name = "C1qtnf7 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Il6st <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Il6st", useMatrix = "GeneScoreMatrix")
Il6st

plotPDF(Il6st, name = "Il6st browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Irf6 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Irf6", useMatrix = "GeneScoreMatrix")
Irf6

plotPDF(Irf6, name = "Irf6 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ccl2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ccl2", useMatrix = "GeneScoreMatrix")
Ccl2

plotPDF(Ccl2, name = "Ccl2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Irf2bpl <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Irf2bpl", useMatrix = "GeneScoreMatrix")
Irf2bpl

plotPDF(Irf2bpl, name = "Irf2bpl browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Socs5 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Socs5", useMatrix = "GeneScoreMatrix")
Socs5

plotPDF(Socs5, name = "Socs5 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Tnfrsf12a <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Tnfrsf12a", useMatrix = "GeneScoreMatrix")
Tnfrsf12a

plotPDF(Tnfrsf12a, name = "Tnfrsf12a browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ifrd1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ifrd1", useMatrix = "GeneScoreMatrix")
Ifrd1

plotPDF(Ifrd1, name = "Ifrd1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Socs2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Socs2", useMatrix = "GeneScoreMatrix")
Socs2

plotPDF(Socs2, name = "Socs2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Nfil3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Nfil3", useMatrix = "GeneScoreMatrix")
Nfil3

plotPDF(Nfil3, name = "Nfil3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Il34 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Il34", useMatrix = "GeneScoreMatrix")
Il34

plotPDF(Il34, name = "Il34 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Clcf1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Clcf1", useMatrix = "GeneScoreMatrix")
Clcf1

plotPDF(Clcf1, name = "Clcf1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Csf1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Csf1", useMatrix = "GeneScoreMatrix")
Csf1

plotPDF(Csf1, name = "Csf1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

#Ion channels
#1. Tpcn2,Slc9a3,Kcnj9,Trpv3,Best1,Fxyd6,Cacng4,Pkd2,Ubc,Ano4,Tmem206,Kcnd1,Kcnab2,Clcn5,Cacna1a,Scn11a,Tmtm38a,P2rx3,P2rx5,Anxa6,Asic2,

#2. Cacnb4,Garn2,Cacna2d3,Kcna4,Kcnq3,Grin1,Grik4,Tmem63c,Htr3b,Cacng5,Scn10a,Trpv2,Lrrc8c,Kcnq5,Clcn1,Grid2,Kcnf2,Gabrg2,Gabra1,Kcne3,

#3. Nalcn,Glrb,Trpm8,Kcnq2,Gria1,Cacna1b,Kcng1,Kcnma1,Asic3,Gabrb3,Gabrg1,Kcnh7,Kcnh1,Gabra2,Scn8a,Kcnv1,Scn1a,Scn3a,Ryr2,Hcn1,Scn9a,Cacng8,

#4.Kcnh5,Trpc5,Cacna2d1,Gabra5,Kcnk12,Kcnk13,Kcnq2,Kcnma1,Ano1,Kcnk16,Kcna1,Htr3a,Trpv1,Scn11a,Scn9a,Grik4,Kcnb1,Chrma6,Kcnb2,Gria4,Scn10a,Kcnv1,Kcng2,Piezo1,Piezo2,Prdm12 

Tpcn2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Tpcn2", useMatrix = "GeneScoreMatrix")
Tpcn2

plotPDF(Tpcn2, name = "Tpcn2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Slc9a3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Slc9a3", useMatrix = "GeneScoreMatrix")
Slc9a3

plotPDF(Slc9a3, name = "Slc9a3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcnj9 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcnj9", useMatrix = "GeneScoreMatrix")
Kcnj9

plotPDF(Kcnj9, name = "Kcnj9 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Trpv3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Trpv3", useMatrix = "GeneScoreMatrix")
Trpv3

plotPDF(Trpv3, name = "Trpv3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Best1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Best1", useMatrix = "GeneScoreMatrix")
Best1

plotPDF(Best1, name = "Best1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Fxyd6 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Fxyd6", useMatrix = "GeneScoreMatrix")
Fxyd6

plotPDF(Fxyd6, name = "Fxyd6 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Cacng4 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Cacng4", useMatrix = "GeneScoreMatrix")
Cacng4

plotPDF(Cacng4, name = "Cacng4 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Pkd2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Pkd2", useMatrix = "GeneScoreMatrix")
Pkd2

plotPDF(Pkd2, name = "Pkd2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ubc <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ubc", useMatrix = "GeneScoreMatrix")
Ubc

plotPDF(Ubc, name = "Ubc browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ano4 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ano4", useMatrix = "GeneScoreMatrix")
Ano4

plotPDF(Ano4, name = "Ano4 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcnd1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcnd1", useMatrix = "GeneScoreMatrix")
Kcnd1

plotPDF(Kcnd1, name = "Kcnd1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcnab2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcnab2", useMatrix = "GeneScoreMatrix")
Kcnab2

plotPDF(Kcnab2, name = "Kcnab2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Clcn5 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Clcn5", useMatrix = "GeneScoreMatrix")
Clcn5

plotPDF(Clcn5, name = "Clcn5 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Cacna1a <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Cacna1a", useMatrix = "GeneScoreMatrix")
Cacna1a

plotPDF(Cacna1a, name = "Cacna1a browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Scn11a <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Scn11a", useMatrix = "GeneScoreMatrix")
Scn11a

plotPDF(Scn11a, name = "Scn11a browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

P2rx3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "P2rx3", useMatrix = "GeneScoreMatrix")
P2rx3

plotPDF(P2rx3, name = "P2rx3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

P2rx5 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "P2rx5", useMatrix = "GeneScoreMatrix")
P2rx5

plotPDF(P2rx5, name = "P2rx5 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Anxa6 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Anxa6", useMatrix = "GeneScoreMatrix")
Anxa6

plotPDF(Anxa6, name = "Anxa6 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Asic2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Asic2", useMatrix = "GeneScoreMatrix")
Asic2

plotPDF(Asic2, name = "Asic2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Cacnb4 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Cacnb4", useMatrix = "GeneScoreMatrix")
Cacnb4

plotPDF(Cacnb4, name = "Cacnb4 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Cacna2d3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Cacna2d3", useMatrix = "GeneScoreMatrix")
Cacna2d3

plotPDF(Cacna2d3, name = "Cacna2d3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcna4 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcna4", useMatrix = "GeneScoreMatrix")
Kcna4

plotPDF(Kcna4, name = "Kcna4 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcnq3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcnq3", useMatrix = "GeneScoreMatrix")
Kcnq3

plotPDF(Kcnq3, name = "Kcnq3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Grin1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Grin1", useMatrix = "GeneScoreMatrix")
Grin1

plotPDF(Grin1, name = "Grin1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Grik4 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Grik4", useMatrix = "GeneScoreMatrix")
Grik4

plotPDF(Grik4, name = "Grik4 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Tmem63c <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Tmem63c", useMatrix = "GeneScoreMatrix")
Tmem63c

plotPDF(Tmem63c, name = "Tmem63c browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Htr3b <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Htr3b", useMatrix = "GeneScoreMatrix")
Htr3b

plotPDF(Htr3b, name = "Htr3b browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Cacng5 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Cacng5", useMatrix = "GeneScoreMatrix")
Cacng5

plotPDF(Cacng5, name = "Cacng5 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Scn10a <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Scn10a", useMatrix = "GeneScoreMatrix")
Scn10a

plotPDF(Scn10a, name = "Scn10a browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Trpv2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Trpv2", useMatrix = "GeneScoreMatrix")
Trpv2

plotPDF(Trpv2, name = "Trpv2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Lrrc8c <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Lrrc8c", useMatrix = "GeneScoreMatrix")
Lrrc8c

plotPDF(Lrrc8c, name = "Lrrc8c browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcnq5 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcnq5", useMatrix = "GeneScoreMatrix")
Kcnq5

plotPDF(Kcnq5, name = "Kcnq5 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Clcn1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Clcn1", useMatrix = "GeneScoreMatrix")
Clcn1

plotPDF(Clcn1, name = "Clcn1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Grid2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Grid2", useMatrix = "GeneScoreMatrix")
Grid2

plotPDF(Grid2, name = "Grid2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Gabrg2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Gabrg2", useMatrix = "GeneScoreMatrix")
Gabrg2

plotPDF(Gabrg2, name = "Gabrg2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Gabra1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Gabra1", useMatrix = "GeneScoreMatrix")
Gabra1

plotPDF(Gabra1, name = "Gabra1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcne3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcne3", useMatrix = "GeneScoreMatrix")
Kcne3

plotPDF(Kcne3, name = "Kcne3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Nalcn <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Nalcn", useMatrix = "GeneScoreMatrix")
Nalcn

plotPDF(Nalcn, name = "Nalcn browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Glrb <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Glrb", useMatrix = "GeneScoreMatrix")
Glrb

plotPDF(Glrb, name = "Glrb browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Trpm8 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Trpm8", useMatrix = "GeneScoreMatrix")
Trpm8

plotPDF(Trpm8, name = "Trpm8 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcnq2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcnq2", useMatrix = "GeneScoreMatrix")
Kcnq2

plotPDF(Kcnq2, name = "Kcnq2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Gria1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Gria1", useMatrix = "GeneScoreMatrix")
Gria1

plotPDF(Gria1, name = "Gria1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Cacna1b <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Cacna1b", useMatrix = "GeneScoreMatrix")
Cacna1b

plotPDF(Cacna1b, name = "Cacna1b browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcng1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcng1", useMatrix = "GeneScoreMatrix")
Kcng1

plotPDF(Kcng1, name = "Kcng1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcnma1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcnma1", useMatrix = "GeneScoreMatrix")
Kcnma1

plotPDF(Kcnma1, name = "Kcnma1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Asic3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Asic3", useMatrix = "GeneScoreMatrix")
Asic3

plotPDF(Asic3, name = "Asic3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Gabrb3 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Gabrb3", useMatrix = "GeneScoreMatrix")
Gabrb3

plotPDF(Gabrb3, name = "Gabrb3 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Gabrg1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Gabrg1", useMatrix = "GeneScoreMatrix")
Gabrg1

plotPDF(Gabrg1, name = "Gabrg1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcnh7 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcnh7", useMatrix = "GeneScoreMatrix")
Kcnh7

plotPDF(Kcnh7, name = "Kcnh7 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcnh1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcnh1", useMatrix = "GeneScoreMatrix")
Kcnh1

plotPDF(Kcnh1, name = "Kcnh1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Gabra2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Gabra2", useMatrix = "GeneScoreMatrix")
Gabra2

plotPDF(Gabra2, name = "Gabra2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Scn8a <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Scn8a", useMatrix = "GeneScoreMatrix")
Scn8a

plotPDF(Scn8a, name = "Scn8a browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcnv1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcnv1", useMatrix = "GeneScoreMatrix")
Kcnv1

plotPDF(Kcnv1, name = "Kcnv1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Scn1a <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Scn1a", useMatrix = "GeneScoreMatrix")
Scn1a

plotPDF(Scn1a, name = "Scn1a browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Scn3a <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Scn3a", useMatrix = "GeneScoreMatrix")
Scn3a

plotPDF(Scn3a, name = "Scn3a browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ryr2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ryr2", useMatrix = "GeneScoreMatrix")
Ryr2

plotPDF(Ryr2, name = "Ryr2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Hcn1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Hcn1", useMatrix = "GeneScoreMatrix")
Hcn1

plotPDF(Hcn1, name = "Hcn1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Scn9a <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Scn9a", useMatrix = "GeneScoreMatrix")
Scn9a

plotPDF(Scn9a, name = "Scn9a browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Cacng8 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Cacng8", useMatrix = "GeneScoreMatrix")
Cacng8

plotPDF(Cacng8, name = "Cacng8 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcnh5 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcnh5", useMatrix = "GeneScoreMatrix")
Kcnh5

plotPDF(Kcnh5, name = "Kcnh5 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Trpc5 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Trpc5", useMatrix = "GeneScoreMatrix")
Trpc5

plotPDF(Trpc5, name = "Trpc5 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Cacna2d1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Cacna2d1", useMatrix = "GeneScoreMatrix")
Cacna2d1

plotPDF(Cacna2d1, name = "Cacna2d1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Gabra5 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Gabra5", useMatrix = "GeneScoreMatrix")
Gabra5

plotPDF(Gabra5, name = "Gabra5 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcnk12 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcnk12", useMatrix = "GeneScoreMatrix")
Kcnk12

plotPDF(Kcnk12, name = "Kcnk12 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcnk13 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcnk13", useMatrix = "GeneScoreMatrix")
Kcnk13

plotPDF(Kcnk13, name = "Kcnk13 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcnq2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcnq2", useMatrix = "GeneScoreMatrix")
Kcnq2

plotPDF(Kcnq2, name = "Kcnq2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Ano1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Ano1", useMatrix = "GeneScoreMatrix")
Ano1

plotPDF(Ano1, name = "Ano1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcnk16 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcnk16", useMatrix = "GeneScoreMatrix")
Kcnk16

plotPDF(Kcnk16, name = "Kcnk16 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcna1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcna1", useMatrix = "GeneScoreMatrix")
Kcna1

plotPDF(Kcna1, name = "Kcna1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Htr3a <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Htr3a", useMatrix = "GeneScoreMatrix")
Htr3a

plotPDF(Htr3a, name = "Htr3a browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Trpv1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Trpv1", useMatrix = "GeneScoreMatrix")
Trpv1

plotPDF(Trpv1, name = "Trpv1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Scn11a <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Scn11a", useMatrix = "GeneScoreMatrix")
Scn11a

plotPDF(Scn11a, name = "Scn11a browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcnb1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcnb1", useMatrix = "GeneScoreMatrix")
Kcnb1

plotPDF(Kcnb1, name = "Kcnb1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcnb2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcnb2", useMatrix = "GeneScoreMatrix")
Kcnb2

plotPDF(Kcnb2, name = "Kcnb2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Gria4 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Gria4", useMatrix = "GeneScoreMatrix")
Gria4

plotPDF(Gria4, name = "Gria4 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Scn10a <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Scn10a", useMatrix = "GeneScoreMatrix")
Scn10a

plotPDF(Scn10a, name = "Scn10a browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcnv1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcnv1", useMatrix = "GeneScoreMatrix")
Kcnv1

plotPDF(Kcnv1, name = "Kcnv1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Kcng2 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Kcng2", useMatrix = "GeneScoreMatrix")
Kcng2

plotPDF(Kcng2, name = "Kcng2 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)

Piezo1 <- plotBrowserTrack(ArchRProj = projPAG5, groupBy = "Sample", useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), geneSymbol = "Piezo1", useMatrix = "GeneScoreMatrix")
Piezo1

plotPDF(Piezo1, name = "Piezo1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)



### PAG6 ###
# add GeneIntegration Matrix
# add Peak2GeneLinks
# plot Peak2GeneHeatmap by Sample
projPAG6 <- projPAG5


#~5 minutes
projPAG6 <- addGeneIntegrationMatrix(
  ArchRProj = projPAG6, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  addToArrow = TRUE,
  force= TRUE,
  groupList = groupList,
  groupRNA = "BioClassification",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore"
)



projPAG6 <- addPeak2GeneLinks(
  ArchRProj = projPAG6,
  reducedDims = "IterativeLSI"
)


p2g <- getPeak2GeneLinks(
  ArchRProj = projPAG5,
  corCutOff = 0.45,
  resolution = 1000,
  returnLoops = TRUE
)

p2g[[1]]

p <- plotPeak2GeneHeatmap(ArchRProj = projPAG5, groupBy = "Sample")


### Upset Plot ###

install.packages("UpSetR")
library(UpSetR)

Upset <- getMatrixFromProject(ArchRProj = projPAG5, useMatrix = "PeakMatrix")
UpsetPlot <- plotMarkerHeatmap(Upset, )



UpsetPlot <- plotMarkerHeatmap(Upset, clusterCols = FALSE)


PeaksHeatmapSaline_Veh <- plotMarkerHeatmap(PeaksSaline_Veh, plotLog2FC = TRUE, clusterCols = FALSE)
PeaksHeatmapCFA_Veh <- plotMarkerHeatmap(PeaksCFA_Veh, plotLog2FC = TRUE, clusterCols = FALSE, returnMatrix = TRUE)

peaksCFA_Veh <- as.data.frame(PeaksHeatmapCFA_Veh)

write.csv(peaksCFA_Veh, "~/Documents/PAGAnalgesicATAC/UpSetPlot.csv", row.names=TRUE)
UpsetPlot <- read.csv("~/Documents/PAGAnalgesicATAC/UpSetPlot.csv")



p <- upset(UpsetPlot, order.by = "freq")
p







projPAG5 <- addMotifAnnotations(ArchRProj = projPAG5, motifSet = "encode", name = "Motif", force = TRUE)

projPAG5 <- addBgdPeaks(projPAG5)

projPAG5 <- addDeviationsMatrix(
  ArchRProj = projPAG4, 
  peakAnnotation = "Motif",
)

plotVarDev <- getVarDeviations(projPAG5, name = "MotifMatrix", plot = TRUE)




Cacna1a <- plotBrowserTrack(ArchRProj = projPAG5, 
                            groupBy = "Sample", 
                            useGroups = c("CFA_3DDA", "CFA_ApAP", "CFA_Veh", "Saline_Veh"), 
                            geneSymbol = "Cacna1a", 
                            useMatrix = "GeneScoreMatrix", 
                            loops = getCoAccessibility(projPAG5)
                            )


Cacna1a


plotPDF(Piezo1, name = "Piezo1 browser track", ArchRProj = projPAG5, addDOC = FALSE, height = 10, width = 12)




coaccessibility <- getCoAccessibility(
  ArchRProj = projPAG5,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = TRUE
)















