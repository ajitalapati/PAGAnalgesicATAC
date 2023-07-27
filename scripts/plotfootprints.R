setwd("~/Documents/PAGAnalgesicATAC")
library(ArchR)
library(parallel)
addArchRThreads(threads = 22) 

library(BSgenome.Rnorvegicus.UCSC.rn7)
library(TxDb.Rnorvegicus.UCSC.rn7.refGene)
library(org.Rn.eg.db)

library(Cairo)

projPAG5 <- loadArchRProject(path = "./PAG_ATAC_Directory4", force = FALSE, showLogo = FALSE)

motifPositions <- getPositions(projPAG5)

# Transcription Factors (TF)
tf_list <- lapply(list(
  "Egr1", "Egr2", "Elk4", "Ezh2", "Hdac10", "Irf3", "Irf6", "Mzf1", "Nfatc1", "Nfatc4", 
  "Nfkb2", "Notch1", "Notch4", "Pou6f1", "Rfx2", "Rreb1", "Tcf4", "Zfp36l2", "Creb3", "Grip1",
  "Sox11", "Zbtb7a", "Foxo6", "Camta1", "Ebf4", "Kcnip3", "Ovol2", "Pou4f3", "Rxrg", "Crebl2", 
  "Egr1", "Atf3", "Crem", "Ddit3", "Sox11", "Bach1", "Stat5a", "Snai2", "Fosl2", "Arid5a", "Sox5a"
), toupper)

# GPCRs
gpcr_list <- lapply(list(
  "Tas1r1", "Vom2r44", "Gpr157", "Gpr81", "Npr3", "C5ar1", "Lpar1", "Ptger1", "Fzd2", "Fzd6", 
  "Fzd7", "Lpar6", "Smo", "Gpr176", "Ptgfr", "Gabbr2", "Adgrl2", "Gpr85", "Gpr68", "Omg", "Mrgprx1", "Gpr149",
  "Mrgpre", "Adgrg2", "Htr1b", "Grm7", "Chrm3", "F2rl2", "Cnr1", "Prokr2", "Oprl1", "Bdkrb2", "Oprm1", 
  "Htr2a", "Ptger4", "Chrm2", "Galr1", "Gpr158", "Nyp2r", "Htr2b", "Olr1746", "Cckar", "Gpr153",
  "Gpr158", "Adra1b", "Ackr1", "Adgrd1", "Cckbr", "Gpr151", "Qrfpr", "Gpr19", "Galr1", "Ntsr2", "F2rl2", 
  "Ptger1", "Gpr68", "Adra2c", "Gpr149", "Gpr45", "Htr1a", "Oprk1", "Mrgpr", "Mrgprx3", "Rgs9"
), toupper)

# Ion channels
ion_channels_list <- lapply(list(
  "Tpcn2", "Slc9a3", "Kcnj9", "Trpv3", "Best1", "Fxyd6", "Cacng4", "Pkd2", "Ubc", "Ano4", "Tmem206", 
  "Kcnd1", "Kcnab2", "Clcn5", "Cacna1a", "Scn11a", "Tmtm38a", "P2rx3", "P2rx5", "Anxa6", "Asic2",
  "Cacnb4", "Garn2", "Cacna2d3", "Kcna4", "Kcnq3", "Grin1", "Grik4", "Tmem63c", "Htr3b", "Cacng5", "Scn10a", 
  "Trpv2", "Lrrc8c", "Kcnq5", "Clcn1", "Grid2", "Kcnf2", "Gabrg2", "Gabra1", "Kcne3",
  "Nalcn", "Glrb", "Trpm8", "Kcnq2", "Gria1", "Cacna1b", "Kcng1", "Kcnma1", "Asic3", "Gabrb3", "Gabrg1", 
  "Kcnh7", "Kcnh1", "Gabra2", "Scn8a", "Kcnv1", "Scn1a", "Scn3a", "Ryr2", "Hcn1", "Scn9a", "Cacng8",
  "Kcnh5", "Trpc5", "Cacna2d1", "Gabra5", "Kcnk12", "Kcnk13", "Kcnq2", "Kcnma1", "Ano1", "Kcnk16", "Kcna1", 
  "Htr3a", "Trpv1", "Scn11a", "Scn9a", "Grik4", "Kcnb1", "Chrma6", "Kcnb2", "Gria4", "Scn10a", "Kcnv1", "Kcng2", "Piezo1", "Piezo2", "Prdm12"
), toupper)

# Growth Factors and receptors
growth_factors_list <- lapply(list(
  "Fgf3", "Mydgf", "Tgfb1", "Fgf11", "Igf2r", "Fgfr1", "Tgfa", "Hgf", "Pdgfa", "Fgf13", "Fgf1", 
  "Vegfa", "Ngfr", "Vegfb", "Igf1r", "Hrh1"
), toupper)

# Cytokines and related molecules
cytokines_list <- lapply(list(
  "C1qtnf7", "Il6st", "Irf2bp2", "Irf6", "Ccl2", "Irf2bpl", "Socs5", "Tnfrsf12a", "Ifrd1", "Socs2", 
  "Nfil3", "Il34", "Fam19a4", "Clcf1", "Ifitm10", "Csf1", "Ackr1"
), toupper)


markerMotifs <- unlist(lapply(ion_channels_list, function(x) grep(x, names(motifPositions), value = TRUE)))

seFoot <- getFootprints(
  ArchRProj = projPAG5, 
  positions = motifPositions[markerMotifs], 
  groupBy = "Sample"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = projPAG5, 
  normMethod = "Subtract",
  plotName = "ion-channels-Footprints-Subtract-Bias-By-Sample",
  addDOC = FALSE,
  smoothWindow = 5
)
