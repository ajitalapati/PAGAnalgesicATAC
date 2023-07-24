setwd("/work/aalapa/pag")

library(dplyr)
library(Seurat)
library(patchwork)

# Load the PBMC dataset
pag.data <- Read10X(data.dir = "raw_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
PAG <- CreateSeuratObject(counts = pag.data$'Gene Expression', min.cells = 3, min.features = 200, names.delim = "-", names.field = 2)

PAG[["percent.mt"]] <- PercentageFeatureSet(PAG, pattern = "^Mt-")

VlnPlot(PAG, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


plot1 <- FeatureScatter(PAG, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(PAG, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

PAG <- subset(PAG, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

PAG <- subset(PAG, idents = c("1", "2"))

PAG <- FindVariableFeatures(PAG, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(PAG), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(PAG)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

all.genes <- rownames(PAG)
PAG <- ScaleData(PAG, features = all.genes)

PAG <- RunPCA(PAG, features = VariableFeatures(object = PAG))

VizDimLoadings(PAG, dims = 1:2, reduction = "pca")

DimPlot(PAG, reduction = "pca")

DimHeatmap(PAG, dims = 1, cells = 500, balanced = TRUE)

DimHeatmap(PAG, dims = 1:15, cells = 500, balanced = TRUE)

ElbowPlot(PAG)

PAG <- FindNeighbors(PAG, dims = 1:15)
PAG <- FindClusters(PAG, resolution = 0.5)

PAG <- RunUMAP(PAG, dims = 1:15)



# load libraries and functions
lapply(c("dplyr","Seurat","HGNChelper","openxlsx"), library, character.only = T)
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R"); source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# get cell-type-specific gene sets from our in-built database (DB)
gs_list = gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Immune system") # e.g. Immune system, Liver, Pancreas, Kidney, Eye, Brain

# assign cell types
scRNAseqData = readRDS(gzcon(url('https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/exampleData.RDS'))); #load example scRNA-seq matrix
es.max = sctype_score(scRNAseqData = scRNAseqData, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# View results, cell-type by cell matrix. See the complete example below
View(es.max)
# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Brain" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)

# get cell-type by cell matrix
es.max = sctype_score(scRNAseqData = PAG[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix. 
# In case Seurat is used, it is either pbmc[["RNA"]]@scale.data (default), pbmc[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or pbmc[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

# merge by cluster
cL_resutls = do.call("rbind", lapply(unique(PAG@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(PAG@meta.data[PAG@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(PAG@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
print(sctype_scores[,1:3])



PAG@meta.data$customclassif = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  PAG@meta.data$customclassif[PAG@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(PAG, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'customclassif')        

BiocManager::install("DESeq2")

markers <- FindMarkers(PAG, ident.1 = "1", ident.2 = "2")