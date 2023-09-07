if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
install.packages('Rcpp')
install.packages('dplyr')
install.packages('Seurat')
install.packages('batchbench')
install.packages('devtools')
devtools::install_github('xzhoulab/iDEA')
install.packages("lattice")
BiocManager::install("clusterProfiler", force =T)
BiocManager::install("org.Hs.eg.db")
source("https://bioconductor.org/biocLite.R")
BiocManager::install("KEGG.db")
BiocManager::install("GO.db")
BiocManager::install("DO.db")
BiocManager::install("GenomeInfoDbData")
install.packages("scCustomize")


require(devtools)
devtools::install_version('flexmix', '2.3-13')
devtools::install_github('hms-dbmi/scde', build_vignettes = FALSE)

if(!requireNamespace("remotes", quietly = TRUE))
  install.packages("remotes")
remotes::install_github("YuLab-SMU/clusterProfiler")
library(Seurat)
library(dplyr)
library(cowplot)
library(reticulate)
library(ggplot2)
library(tictoc)
library(future)
library(ggpubr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(KEGG.db)
library(iDEA)
library(vegan)
library(SingleR)
library(tidyverse)
library(RColorBrewer)
library(scCustomize)
set.seed(42)
setwd("~/Documents/scRNA SEQ/mice/micedata")

#read in data
NM23.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/Ailiang_scSeq/NM2-NM3/outs/raw_feature_bc_matrix")
BM2.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/Ailiang_scSeq/BM2/outs/raw_feature_bc_matrix")
TM2.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/Ailiang_scSeq/TM2/outs/raw_feature_bc_matrix")
NF3.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/Ailiang_scSeq/NF3/outs/raw_feature_bc_matrix")
BF3.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/Ailiang_scSeq/BF3/outs/raw_feature_bc_matrix")
TF3.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/Ailiang_scSeq/TF3/outs/raw_feature_bc_matrix")
BM3.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/Ailiang_scSeq/BM3/outs/raw_feature_bc_matrix")
TM3.data <- Read10X(data.dir ="/rsrch3/scratch/genomic_med/dzamler/Ailiang_scSeq/TM3/outs/raw_feature_bc_matrix")
NF4.data <- Read10X(data.dir ="NF4")
BF4.data <- Read10X(data.dir ="BF4")
TF4.data <- Read10X(data.dir ="TF4")
NF6.data <- Read10X(data.dir ="NF6")
BF6.data <- Read10X(data.dir ="BF6")
TF6.data <- Read10X(data.dir ="TF6")

#create Suerat objects
NM23 <- CreateSeuratObject(counts = NM23.data, project = "NM2&3", min.cells = 3)
NM23$num <- "2&3"
NM23$BTN <- "Periphery"
NM23$MF <- "Male"
NM23[["percent.mt"]] <- PercentageFeatureSet(NM23, pattern = "^mt-")
NM23 <- subset(NM23, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
VlnPlot(NM23, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

BM2 <- CreateSeuratObject(counts = BM2.data, project = "BM2", min.cells = 3)
BM2$num <- "2"
BM2$BTN <- "Border"
BM2$MF <- "Male"
BM2[["percent.mt"]] <- PercentageFeatureSet(BM2, pattern = "^mt-")
BM2 <- subset(BM2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
VlnPlot(BM2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

TM2 <- CreateSeuratObject(counts = TM2.data, project = "TM2", min.cells = 3)
TM2$num <- "2"
TM2$BTN <- "Core"
TM2$MF <- "Male"
TM2[["percent.mt"]] <- PercentageFeatureSet(TM2, pattern = "^mt-")
TM2 <- subset(TM2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
VlnPlot(TM2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

NF3 <- CreateSeuratObject(counts = NF3.data, project = "NF3", min.cells = 3)
NF3$num <- "3"
NF3$BTN <- "Periphery"
NF3$MF <- "Female"
NF3[["percent.mt"]] <- PercentageFeatureSet(object = NF3, pattern = "^mt-")
NF3 <- subset(NF3, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
VlnPlot(NF3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

BF3 <- CreateSeuratObject(counts = BF3.data, project = "BF3", min.cells = 3)
BF3$num <- "3"
BF3$BTN <- "Border"
BF3$MF <- "Female"
BF3[["percent.mt"]] <- PercentageFeatureSet(BF3, pattern = "^mt-")
BF3 <- subset(BF3, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
VlnPlot(BF3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

TF3 <- CreateSeuratObject(counts = TF3.data, project = "TF3", min.cells = 3)
TF3$num <- "3"
TF3$BTN <- "Core"
TF3$MF <- "Female"
TF3[["percent.mt"]] <- PercentageFeatureSet(TF3, pattern = "^mt-")
TF3 <- subset(TF3, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
VlnPlot(TF3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

BM3 <- CreateSeuratObject(counts = BM3.data, project = "BM3", min.cells = 3)
BM3$num <- "3"
BM3$BTN <- "Border"
BM3$MF <- "Male"
BM3[["percent.mt"]] <- PercentageFeatureSet(BM3, pattern = "^mt-")
BM3 <- subset(BM3, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
VlnPlot(BM3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

TM3 <- CreateSeuratObject(counts = TM3.data, project = "TM3", min.cells = 3)
TM3$num <- "3"
TM3$BTN <- "Core"
TM3$MF <- "Male"
TM3[["percent.mt"]] <- PercentageFeatureSet(TM3, pattern = "^mt-")
TM3 <- subset(TM3, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
VlnPlot(TM3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


NF4 <- CreateSeuratObject(counts = NF4.data, project = "NF4", min.cells = 3)
NF4$num <- "4"
NF4$BTN <- "Periphery"
NF4$MF <- "Female"
NF4[["percent.mt"]] <- PercentageFeatureSet(NF4, pattern = "^mt-")
FeatureScatter(NF4,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
NF4 <- subset(NF4, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
VlnPlot(NF4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


BF4 <- CreateSeuratObject(counts = BF4.data, project = "BF4", min.cells = 3)
BF4$num <- "4"
BF4$BTN <- "Border"
BF4$MF <- "Female"
BF4[["percent.mt"]] <- PercentageFeatureSet(BF4, pattern = "^mt-")
FeatureScatter(BF4,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
BF4 <- subset(BF4, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
VlnPlot(BF4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


TF4 <- CreateSeuratObject(counts = TF4.data, project = "TF4", min.cells = 3)
TF4$num <- "4"
TF4$BTN <- "Core"
TF4$MF <- "Female"
TF4[["percent.mt"]] <- PercentageFeatureSet(TF4, pattern = "^mt-")
FeatureScatter(TF4,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
TF4 <- subset(TF4, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
VlnPlot(TF4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


NF6 <- CreateSeuratObject(counts = NF6.data, project = "NF6", min.cells = 3)
NF6$num <- "6"
NF6$BTN <- "Periphery"
NF6$MF <- "Female"
NF6[["percent.mt"]] <- PercentageFeatureSet(NF6, pattern = "^mt-")
FeatureScatter(NF6,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
NF6 <- subset(NF6, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
VlnPlot(NF6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


BF6 <- CreateSeuratObject(counts = BF6.data, project = "BF6", min.cells = 3)
BF6$num <- "6"
BF6$BTN <- "Border"
BF6$MF <- "Female"
BF6[["percent.mt"]] <- PercentageFeatureSet(BF6, pattern = "^mt-")
FeatureScatter(BF6,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
BF6 <- subset(BF6, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
VlnPlot(BF6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


TF6 <- CreateSeuratObject(counts = TF6.data, project = "TF6", min.cells = 3)
TF6$num <- "6"
TF6$BTN <- "Core"
TF6$MF <- "Female"
TF6[["percent.mt"]] <- PercentageFeatureSet(TF6, pattern = "^mt-")
FeatureScatter(TF6,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
TF6 <- subset(TF6, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
VlnPlot(TF6, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#merge dataset
datasets <- list(NM23, NF3, NF4,NF6,BM2,BM3,BF3,BF4,BF6,TM2,TM3,TF3,TF4,TF6)
datasets=lapply( datasets, function(x) SCTransform(x,vars.to.regress = "percent.mt", verbose = TRUE))
immune.anchors <- FindIntegrationAnchors(object.list = datasets, dims = 1:20)
immune.integrated <- IntegrateData(anchorset = immune.anchors, dims = 1:20, k.weight = 85)
immune.integrated <- ScaleData(immune.integrated, features = rownames(immune.integrated))
immune.integrated <- RunPCA(immune.integrated, npcs = 30, verbose = FALSE)
immune.integrated <- RunUMAP(immune.integrated, reduction = "pca", dims = 1:30)
immune.integrated <- FindNeighbors(immune.integrated, reduction = "pca", dims = 1:30)
DefaultAssay(immune.integrated) <- "integrated"

#Change resolution to modify number or clusters
immune.integrated<- FindClusters(immune.integrated, resolution = 0.50)
immune.integrated$BTN <- factor(immune.integrated$BTN, levels = c('Periphery', 'Border',"Core"))
DimPlot(immune.integrated, reduction = "umap", label = TRUE, pt.size = 0.1,cols=DiscretePalette_scCustomize(num_colors = 30, palette = "varibow",shuffle_pal = TRUE,seed =16 ),split.by = "BTN") + NoLegend()
DimPlot_scCustom(immune.integrated, reduction = "umap", label = TRUE, pt.size = 0.1,split.by ="BTN" ,colors_use= "polychrome")+NoLegend()

#Find markers for every cluster compared to all remaining cells, report only the positive ones
immune.markers <- FindAllMarkers(immune.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Predicting cell types
immune.integrated_for_SingleR <- GetAssayData(immune.integrated, slot="scale.data")
clusters=immune.integrated@meta.data$seurat_clusters
pred.mouseImmu <- SingleR(test = immune.integrated_for_SingleR, ref = mouseImmu, labels = mouseImmu$label.fine,
                          method = "cluster", clusters = clusters,
                          assay.type.test = "logcounts", assay.type.ref = "logcounts")
pred.mouseRNA <- SingleR(test = immune.integrated_for_SingleR, ref = mouseRNA, labels = mouseRNA$label.fine ,
                         method = "cluster", clusters = clusters,
                         assay.type.test = "logcounts", assay.type.ref = "logcounts")
cellType=data.frame(ClusterID=levels(immune.integrated@meta.data$seurat_clusters),
                    mouseImmu=pred.mouseImmu$labels,
                    mouseRNA=pred.mouseRNA$labels )
cellType


#######Validate cell types
#t cell
FeaturePlot(immune.integrated, features = c("Cd3d", "Cd3e", "Cd3g"), min.cutoff = "q9")
#####nk
FeaturePlot(immune.integrated, features = c("Klrd1", "Nkg7", "Nktr"), min.cutoff = "q9")
VlnPlot(immune.integrated, features = c("Klrd1", "Nkg7", "Nktr"), pt.size = 0)
#CD4, CD8 T
VlnPlot(immune.integrated, features = c("Cd4", "Cd8a", "Cd8b"), pt.size = 0)
###Gamma-delta T
FeaturePlot(immune.integrated, features = c("Tcrg-C1"), min.cutoff = "q9")
#####cyotxic nk
VlnPlot(immune.integrated, features = c("Prf1","Cx3cr1","Fcgr3","Ncam1","Cxcr1","Itgb2","Klrc2","Klrg1"), pt.size = 0)
############reg nk
VlnPlot(immune.integrated, features = c("Ncam1","Ccr7","CSF2","Cxcr3","Ifng","Il2rb","Il7r","Kit","Klcr1","Klrd1","Ncr1","Sell"), pt.size = 0)
#B cell
VlnPlot(immune.integrated, features = c("Cd79a", "Igkc", "Ebf1"))
FeaturePlot(immune.integrated, reduction = "tsne", features = c("Cd79a", "Igkc", "Ebf1"))
#microglia
VlnPlot(immune.integrated, features = c("Cx3cr1", "P2ry12", "Tmem119" ,"Hexb", "Cst3"), pt.size = 0)
FeaturePlot(immune.integrated, features = c("Cx3cr1", "P2ry12", "Tmem119" ,"Hexb", "Cst3"))
#Neutrophil
VlnPlot(immune.integrated, features = c("Cd24a", "S100a8", "S100a9"))
FeaturePlot(immune.integrated, features = c("Cd24a", "S100a8", "S100a9"))

###############RENAME THE CLUSTER
new.cluster.ids <- c("TIM2","ERM", "ILC","CD8 T","TIM1",  "cDC",  "Regulatory NK", "Cytotoxic NK", "mo-DC","CD4 T","Treg","B cell","Monocyte","Microglia","pDC",
                     "Gamma-delta T","Neutrophil","Proli. CD8 T")
names(new.cluster.ids) <- levels(immune.integrated)
immune.integrated <- RenameIdents(immune.integrated, new.cluster.ids) 

#Dimplot visulization
DimPlot(immune.integrated, reduction = "umap", label = TRUE, pt.size = 0.1,cols=DiscretePalette_scCustomize(num_colors = 30, palette = "varibow",shuffle_pal = TRUE,seed =16 ),split.by = "BTN") + NoLegend()
DimPlot_scCustom(immune.integrated, reduction = "umap", label = TRUE, pt.size = 0.1,split.by ="BTN" ,colors_use= "polychrome")+NoLegend()

######################################
###Gene Ontology Part#################
######################################
imp.neut1 <- subset(immune.integrated, idents = c("TIM2","ERM"))

agg=imp.neut1

DefaultAssay(agg)="RNA"
immune.markers <- FindAllMarkers(agg, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Mm.eg.db")
library(org.Mm.eg.db)

mMarkers <- immune.markers
# Filter
mMarkers <- mMarkers[ mMarkers$p_val_adj<0.01, ]
# Sort
mMarkers <- mMarkers[ order(mMarkers$cluster, -mMarkers$avg_log2FC), ]

View(mMarkers)

# Set number of upregulated genes
NTop <-50

#Set your number of clusters
agg@active.ident
mGroups <- c("TIM2","ERM")

#Entrez4Enrichment- convert genes to EntrezID for GO conversion
for (group in mGroups)
{
  GeneSymbol4Enrichment<-mMarkers[mMarkers$cluster==group,] %>% group_by(cluster) %>% top_n(n = NTop, wt = avg_log2FC)
  GeneSymbol4Enrichment.group<-GeneSymbol4Enrichment$gene
  converted <- bitr(GeneSymbol4Enrichment.group, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
  Entrez4Enrichment <- converted$ENTREZID
  assign(paste0("mEntrez4Enrichment", group), converted$ENTREZID)
}

# background gene selection
NTopBkgGenes<-2000
mExpressionVector<-rowMeans(as.matrix(GetAssayData(agg, slot = "data")))
mExpressionVector<-sort(mExpressionVector, decreasing = T)
mGeneSymbol4Bkg<-names(mExpressionVector[1:NTopBkgGenes])
mconverted <- bitr(mGeneSymbol4Bkg, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
mEntrez4Bkg <- mconverted$ENTREZID

#choose the gene ontologies you want to check
Onts <- c("CC","MF","BP")

for (group in mGroups)
{
  Entrez_counter <- paste("mEntrez4Enrichment", group, sep = "")
  Entrez_counter <- eval(parse(text = Entrez_counter))
  for (ont in Onts)
  {
    print("Working on")
    print(group)
    print(ont)
    enrich.res<-enrichGO(gene= Entrez_counter, universe = mEntrez4Bkg,
                         ont=ont,OrgDb='org.Mm.eg.db' ,readable=T,
                         minGSSize = 20, maxGSSize = 500, qvalueCutoff = 0.05)
    if(is.null(enrich.res) == T){
      
    }
    else{
      print("computed")
      assign(paste0("enrich.res", group, ont), enrich.res)
    }
  }
}

p1=dotplot(enrich.resERMBP, showCategory=20) +ggtitle("Significantly Enriched Pathways in ERM  Biological Process")
p2=dotplot(enrich.resTIMBP, showCategory=20) +ggtitle("Significantly Enriched Pathways in TIM2 Biological Process")
p3=dotplot(enrich.resERMCC, showCategory=20) +ggtitle("Significantly Enriched Pathways in ERM Cellular Component")
p4=dotplot(enrich.resTIMCC, showCategory=20) +ggtitle("Significantly Enriched Pathways in TIM2 Cellular Component")
p5=dotplot(enrich.resERMMF, showCategory=20) +ggtitle("Significantly Enriched Pathways in ERM Molecular Function")
p6=dotplot(enrich.resTIMMF, showCategory=20) +ggtitle("Significantly Enriched Pathways in TIM2 Molecular Function")

plot_grid(p1, p2)
plot_grid(p3, p4)
plot_grid(p5, p6)


######################################
###CellChat###########################
######################################
install.packages("remotes")
remotes::install_github("sqjin/CellChat")

library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

cellchat <- createCellChat(object = immune.integrated, group.by = "ident", assay = "RNA")
CellChatDB <- CellChatDB.mouse

# use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
# use a subset of CellChatDB for cell-cell communication analysis
# use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB 
# simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) 
# This step is necessary even if using the whole database
future::plan("multiprocess", workers = 4) 
# do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

#############OPTIANL I DIDNOT USE THIS TIME
cellchat <- projectData(cellchat, PPI.human)

options(future.globals.maxSize = 8000 * 1024^2)

cellchat <- computeCommunProb(cellchat,population.size = TRUE,raw.use = FALSE)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- aggregateNet(cellchat)

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways

# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
for (i in 1:length(pathways.show.all)) {
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
}

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
cellchat@net$count
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

write.csv(mat,"cellchatweight.csv")

######################################
###Gene Set Scoring###################
######################################
install.packages("ggpubr")
install.packages("ggsignif")
library(ggsignif)
library(ggpubr)

geneset.list<-list(
  "Adhesion"=c("Ada","Alcam","Amica1","Angpt2","Cd164","Cd22","Cd33","Cd6","Cd63","Cd84","Cd9","Cd96","Cd97","Cdh5","Col3a1","Csf3r","Cyfip2","Fn1","Glycam1","Icam2","Icam4","Irf2","Itga1","Itga2","Itga2b","Itga4","Itga5","Itga6","Itgae","Itgal","Itgax","Itgb1","Itgb4","Jam3","Kdr","Klra4","Klra5","Klra6","Klra7","Ly9","Lyve1","Map2k1","Mcam","Mertk","Mfge8","Mmp9","Msln","Ncam1","Nrp1","Pecam1","Plau","Pnma1","Pvrl2","Saa1","Sell","Selplg","Siglec1","Spink5","Spn","Spp1","Tek","Thy1","Tnfrsf12a","Vcam1","Vwf"), 
  "Adaptive_immunity"=c("C1qbp","C3ar1","Camp","Ccl1","Ccl11","Ccl12","Ccl24","Ccl25","Ccl26","Ccl3","Ccl4","Ccl5","Ccl7","Ccl8","Ccr1","Ccr2","Ccr4","Ccr5","Ccr6","Ccr7","Ccr8","Ccrl2","Cd28","Cd4","Cd40","Cd40lg","Cd80","Cd86","Cd8a","Cd97","Cklf","Cma1","Creb5","Crp","Csf2","Cxcl1","Cxcl10","Cxcl13","Cxcl14","Cxcl2","Cxcl3","Cxcl9","Cxcr1","Cxcr2","Cxcr3","Fasl","Fcer1a","Fcer2a","Foxp3","Fpr2","Gata3","H2-Q10","Hmgb1","Icam1","Ido1","Ifna2","Ifnar1","Ifnb1","Ifng","Ifngr1","Il10","Il13","Il17a","Il18","Il1a","Il1b","Il1r1","Il2","Il22","Il23a","Il24","Il25","Il2ra","Il4","Il5","Il6","Il6st","Il9","Irf3","Irf7","Itgam","Itk","Jak2","Jam3","Klrb1","Lta","Mapk8","Mbl2","Ms4a2","Mx1","Nfatc4","Nfkb1","Nod2","Nos2","Pnma1","Rag1","Rorc","S100a8","Sele","Slc11a1","Stat1","Stat3","Stat4","Stat6","Tbx21","Thbs1","Tlr4","Tlr6","Tnf","Tnfrsf11a","Txk","Xcr1"),
  "Antigen_processing"=c("Ccr7","Cd1d1","Cd1d2","Cd74","Fcer1g","Fcgr1","Fcgr2b","Fcgr3","H2-Aa","H2-Ab1","H2-D1","H2-DMa","H2-DMb1","H2-DMb2","H2-Ea-ps","H2-Eb1","H2-K1","H2-M3","H2-Ob","H2-Q1","H2-Q10","H2-Q2","H2-T23","Icam1","Mr1","Nod1","Nod2","Psmb8","Psmb9","Relb","Slc11a1","Tap1","Tap2","Tapbp"),
  "Apoptosis"=c("Abl1","Adora2a","Angpt1","Apoe","App","Atg5","Atg7","Atm","Bax","Bcl2l1","Bid","Birc5","Btk","C6","C9","Casp3","Casp8","Cd38","Cd3g","Cd5","Cd59b","Cdk1","Clec5a","Clu","Ctsh","Cyfip2","Cyld","Dusp6","Egfr","Egr3","Ep300","Ets1","Fadd","Fn1","Gpi1","Gzma","Gzmb","Hif1a","Hmgb1","Ifih1","Igf2r","Ikbke","Il19","Il24","Il3","Inpp5d","Itga1","Itga6","Jun","Kdr","Lck","Lcn2","Litaf","Lrp1","Ltbr","Ltk","Map2k4","Map3k1","Map3k5","Map3k7","Mapk1","Mapk3","Mapk8","Mef2c","Mertk","Mfge8","Mmp9","Muc1","Myc","Nefl","Nlrp3","Nos2","Osm","Pdcd1","Pik3cg","Plaur","Pml","Prkcd","Psen1","Ptgs2","Pycard","Rps6","Runx3","S100b","Sell","Spn","Spp1","Tcf7","Tek","Tgfb2","Tgfb3","Tgfbr1","Tgfbr2","Tmem173","Tnfaip3","Tnfrsf10b","Tnfrsf11b","Tnfrsf12a","Tnfrsf18","Tnfrsf8","Tnfsf10","Tnfsf12","Tnfsf14","Tnfsf15","Traf2","Traf3","Trp53","Twist1","Txnip","Vegfa","Xaf1"),
  "Autophagy"=c("Atg10","Atg12","Atg16l1","Atg5","Atg7","Lamp1"),
  "B_activation"=c("Cd28","Cd4","Cd69","Cd70","Cd83","Cd86","Cr2","Cxcr5","Dpp4","Fcer2a","Icosl","Il2ra","Il6","Ms4a1","Tgfb1"),
  "B_differentiation"=c("Ada","Cd79a","Flt3","Il10","Il11","Ptprc","Rag1"),
  "B_function"=c("Ada","Atm","Bcl10","Bcl2","Bcl6","Blk","Blnk","Bmi1","Btla","Card11","Casp3","Cd19","Cd200r1","Cd22","Cd27","Cd28","Cd37","Cd38","Cd4","Cd40","Cd40lg","Cd69","Cd70","Cd74","Cd79a","Cd79b","Cd81","Cd83","Cd86","Cdkn1a","Cr2","Ctla4","Cxcl13","Cxcr5","Dpp4","Fas","Fcer2a","Fcgr2b","Flt3","Foxp3","Gpr183","H2-Q10","Icosl","Ikbkb","Ikbkg","Ikzf1","Il10","Il11","Il13","Il13ra1","Il1r2","Il21","Il2ra","Il2rg","Il4","Il5","Il6","Il7","Il7r","Inpp5d","Itga2","Jak3","Lck","Lyn","Mapk1","Mef2c","Mif","Ms4a1","Nt5e","Prdm1","Prkcd","Ptprc","Rag1","Sh2b2","Stat5b","Syk","Tgfb1","Ticam1","Tirap","Tnfaip3","Tnfrsf13b","Tnfrsf13c","Tnfrsf4","Tnfsf13b","Tnfsf4"),
  "B_proliferation"=c("Bcl2","Cd38","Cd81","Cdkn1a","Ctla4","Il7","Prkcd","Tnfrsf13b","Tnfrsf13c","Tnfsf13b"),
  "B_regulation"=c("Btla","Cd27","Cd40","Cd40lg","Foxp3","Il4"),
  "Bacterial_response"=c("Ccl2","Cd14","Fos","Hmgb1","Il10","Il1b","Irak1","Jun","Lta","Ly86","Ly96","Nfkbia","Ptgs2","Ripk2","Tlr4","Tlr6","Tnfrsf1a"),
  "Basic_cell_function"=c("Arg1","Arg2","Chil3","Chit1","Cmpk2","Ctsg","Ctsl","Ddx60","Dock9","Dusp4","Emr1","Epsti1","Ewsr1","F12","F13a1","Fpr2","Gbp2b","Hamp","Hcst","Herc6","Hsd11b1","Isg15","Isg20","Map2k2","Mapk11","Mpo","Mpped1","Ncf4","Notch1","Oas2","Oas3","Oasl1","Pdgfc","Pla2g1b","Pla2g6","Pmch","Pou2af1","Prg2","Psmb7","Raet1c","Reps1","Rrad","Rsad2","S100a8","St6gal1","Stat2","Tab1","Tank","Timd4","Trem2","Ubc","Usp18","Usp9y","Ythdf2","Zfp13"),
  "Cancer_progression"=c("Akt3","Angpt1","Angpt2","Apoe","C3","C3ar1","Camp","Casp8","Ccl11","Ccl5","Ccl7","Ccl8","Ccr2","Ccr3","Cd163","Cd34","Cd36","Cd44","Cd46","Cdh1","Cdkn1a","Ceacam1","Cfp","Clu","Cma1","Col1a1","Col3a1","Col4a1","Crebbp","Csf2rb","Cspg4","Dll4","Egfr","Erbb2","Fap","Hif1a","Hspb2","Kdr","Mmp9","Msln","Psma2","Sele","Smad2","Smad3","Smad4","Smn1","Snai1","Tdo2","Tek","Tgfbr1","Tgfbr2","Tie1","Twist1","Vegfc","Vhl","Vim","Vwf"),
  "CD4_T_function"=c("Cd27","Cd4","Crebbp","Ctla4","Il15","Il7","Jak2","Mapk8","Ptprc","Socs3","Tgfb3","Tnfrsf4","Tnfrsf8","Tnfsf4","Tyk2","Yy1"),
  "CD_molecules"=c("Abcb1a","Alcam","Bst1","Bst2","Btla","C1qbp","Ccr1","Ccr2","Ccr3","Ccr4","Ccr5","Ccr6","Ccr7","Ccr9","Cd14","Cd160","Cd163","Cd164","Cd180","Cd19","Cd1d1","Cd1d2","Cd2","Cd200","Cd207","Cd209e","Cd22","Cd244","Cd247","Cd27","Cd274","Cd276","Cd28","Cd33","Cd34","Cd36","Cd37","Cd38","Cd3d","Cd3e","Cd3eap","Cd3g","Cd4","Cd40","Cd40lg","Cd46","Cd47","Cd48","Cd5","Cd53","Cd55","Cd6","Cd63","Cd68","Cd7","Cd70","Cd74","Cd79a","Cd79b","Cd80","Cd81","Cd83","Cd84","Cd86","Cd8a","Cd8b1","Cd9","Cd96","Cd97","Cd99","Cdh1","Cdh5","Ceacam1","Cr2","Csf1r","Csf2rb","Csf3r","Ctla4","Ctsw","Cxcr1","Cxcr2","Cxcr3","Cxcr4","Cxcr5","Cxcr6","Dpp4","Eng","Entpd1","Epcam","Fas","Fcer2a","Fcgr1","Fcgr2b","Fcgr3","Fcgr4","Flt3","Icam1","Icam2","Icam4","Icos","Icosl","Ifitm1","Ifngr1","Igf1r","Igf2r","Igll1","Il10ra","Il12rb1","Il13ra1","Il13ra2","Il15ra","Il17ra","Il18r1","Il18rap","Il1r1","Il1r2","Il21r","Il2ra","Il2rb","Il2rg","Il3ra","Il4ra","Il5ra","Il6ra","Il6st","Il7r","Itga1","Itga2","Itga2b","Itga4","Itga5","Itga6","Itgae","Itgal","Itgam","Itgax","Itgb1","Itgb2","Itgb3","Itgb4","Kit","Klrb1","Klrc1","Klrc2","Klrk1","Lag3","Lamp1","Lamp2","Lamp3","Lilra5","Lrp1","Lrrn3","Ly9","Mcam","Mme","Mrc1","Ms4a1","Msr1","Mst1r","Muc1","Ncam1","Ncr1","Nrp1","Nt5e","Pdcd1","Pdcd1lg2","Pdgfrb","Pecam1","Plaur","Psmd7","Ptgdr2","Ptprc","Pvr","Sele","Sell","Selplg","Siglec1","Slamf1","Slamf6","Slamf7","Spn","Tfrc","Thbd","Thy1","Tlr1","Tlr2","Tlr3","Tlr4","Tlr6","Tlr8","Tlr9","Tnfrsf10b","Tnfrsf11a","Tnfrsf12a","Tnfrsf13b","Tnfrsf13c","Tnfrsf14","Tnfrsf17","Tnfrsf1a","Tnfrsf1b","Tnfrsf4","Tnfrsf8","Tnfrsf9","Tnfsf10","Tnfsf11","Tnfsf13","Tnfsf13b","Tnfsf14","Tnfsf4","Tnfsf8","Trem1","Vcam1"),
  "Cell_cycle"=c("Abcb1a","Abl1","Anp32b","Anxa1","App","Atm","Bid","Birc5","Casp3","Ccnd3","Cdk1","Cdkn1a","Cxcl15","Cyld","Ets1","Il12a","Il12b","Itgb1","Map2k1","Mapk3","Muc1","Myc","Nfatc1","Pin1","Pml","Prkce","Rps6","Runx3","Smpd3","Stat5b","Tal1","Tgfb1","Tgfb2","Thbs1","Trp53","Txnip"),
  "Chemokine_receptors"=c("Ccr1","Ccr2","Ccr3","Ccr4","Ccr5","Ccr6","Ccr7","Ccr8","Ccr9","Ccrl2","Cx3cr1","Cxcr1","Cxcr2","Cxcr3","Cxcr4","Cxcr5","Cxcr6","Xcr1"),
  "Chemokines"=c("Ccl1","Ccl11","Ccl12","Ccl17","Ccl19","Ccl2","Ccl20","Ccl21a","Ccl22","Ccl24","Ccl25","Ccl26","Ccl28","Ccl3","Ccl4","Ccl5","Ccl6","Ccl7","Ccl8","Ccl9","Cklf","Csf1r","Cx3cl1","Cxcl1","Cxcl10","Cxcl11","Cxcl12","Cxcl13","Cxcl14","Cxcl15","Cxcl16","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Elane","Hc","Il18","Il1rl1","Il22ra1","Il4ra","Il6ra","Ilf3","Itch","Jak1","Lbp","Myd88","Sigirr","Ticam1","Tirap","Tlr3","Tlr4","Tlr7","Tlr9","Tnfsf4","Xcl1"),
  "Chemokines_chemokine_receptors"=c("C5ar1","Ccl1","Ccl11","Ccl12","Ccl17","Ccl19","Ccl2","Ccl20","Ccl21a","Ccl22","Ccl24","Ccl25","Ccl26","Ccl28","Ccl3","Ccl4","Ccl5","Ccl6","Ccl7","Ccl8","Ccl9","Ccr1","Ccr2","Ccr3","Ccr4","Ccr5","Ccr6","Ccr7","Ccr8","Ccr9","Ccrl2","Cklf","Cmklr1","Csf1r","Cx3cl1","Cx3cr1","Cxcl1","Cxcl10","Cxcl11","Cxcl12","Cxcl13","Cxcl14","Cxcl15","Cxcl16","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Cxcr1","Cxcr2","Cxcr3","Cxcr4","Cxcr5","Cxcr6","Elane","Hc","Ifng","Il16","Il18","Il1b","Il1rl1","Il22ra1","Il4","Il4ra","Il6","Il6ra","Ilf3","Itch","Jak1","Lbp","Myd88","Ppbp","Sigirr","Tgfb1","Ticam1","Tirap","Tlr3","Tlr4","Tlr7","Tlr9","Tnf","Tnfsf4","Xcl1","Xcr1"),
  "Chemotaxis"=c("C5ar1","Cmklr1","Ifng","Il16","Il1b","Il4","Il6","Ppbp","Tgfb1","Tnf"),
  "Complement_pathway"=c("A2m","C1qa","C1qb","C1qbp","C1ra","C1s1","C2","C3","C3ar1","C4b","C5ar1","C6","C7","C8a","C8b","C8g","C9","Cd46","Cd55","Cd59b","Cfb","Cfd","Cfh","Cfi","Cfp","Cr2","Crp","Hc","Masp1","Masp2","Mbl2","Serping1"),
  "Cytokine_receptors"=c("Ccr1","Ccr2","Ccr3","Ccr5","Ccr9","Csf1r","Cxcr1","Cxcr4","Erbb2","Flt3","Lyn"),
  "Cytokines"=c("Bcl2l1","Birc5","Card11","Casp1","Ccl1","Ccl11","Ccl12","Ccl17","Ccl19","Ccl2","Ccl20","Ccl22","Ccl24","Ccl25","Ccl26","Ccl27a","Ccl28","Ccl3","Ccl4","Ccl6","Ccl7","Ccl8","Ccl9","Cd14","Cd40lg","Cd44","Cd70","Cklf","Clec4n","Clec5a","Col3a1","Csf1","Csf2","Csf2rb","Csf3","Csf3r","Cx3cl1","Cx3cr1","Cxcl1","Cxcl10","Cxcl11","Cxcl12","Cxcl13","Cxcl14","Cxcl15","Cxcl16","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Cxcr2","Cxcr3","Ebi3","Elane","F2rl1","Fasl","Flt3l","Foxp3","Gpi1","Hc","Ido1","Ifna1","Ifna2","Ifna4","Ifnar2","Ifnb1","Ifng","Ifngr1","Igf1r","Il10","Il11","Il11ra1","Il12a","Il12b","Il12rb1","Il12rb2","Il13","Il13ra1","Il13ra2","Il15","Il16","Il17a","Il17b","Il17f","Il17ra","Il17rb","Il18","Il18rap","Il19","Il1a","Il1r1","Il1r2","Il1rap","Il1rapl2","Il1rl1","Il1rl2","Il1rn","Il2","Il21","Il21r","Il22","Il22ra2","Il23a","Il23r","Il24","Il25","Il27","Il2rb","Il2rg","Il34","Il3ra","Il4","Il4ra","Il5","Il5ra","Il6","Il6ra","Il6st","Il7","Il7r","Il9","Irak1","Irak2","Irak3","Irak4","Irf3","Itk","Jak1","Jak2","Jak3","Kit","Lif","Lta","Ltb","Mif","Mme","Ms4a2","Myd88","Nfatc2","Nfkb1","Nod2","Osm","Pin1","Pparg","Prkce","Ptprc","Rel","S100b","Saa1","Sh2b2","Sigirr","Socs1","Spp1","Stat1","Stat4","Stat6","Tgfb1","Tgfb2","Ticam2","Tnf","Tnfrsf11a","Tnfrsf1a","Tnfrsf4","Tnfsf10","Tnfsf11","Tnfsf12","Tnfsf13","Tnfsf13b","Tnfsf14","Tnfsf15","Tnfsf18","Tnfsf4","Tnfsf8","Traf2","Traf3","Traf6","Txk","Vegfa","Xcl1","Xcr1"),
  "Cytokines_cytokine_receptors"=c("Bcl2l1","Birc5","Card11","Casp1","Ccl1","Ccl11","Ccl12","Ccl17","Ccl19","Ccl2","Ccl20","Ccl22","Ccl24","Ccl25","Ccl26","Ccl27a","Ccl28","Ccl3","Ccl4","Ccl6","Ccl7","Ccl8","Ccl9","Ccr1","Ccr2","Ccr3","Ccr5","Ccr9","Cd14","Cd40lg","Cd44","Cd70","Cklf","Clec4n","Clec5a","Col3a1","Csf1","Csf1r","Csf2","Csf2rb","Csf3","Csf3r","Cx3cl1","Cx3cr1","Cxcl1","Cxcl10","Cxcl11","Cxcl12","Cxcl13","Cxcl14","Cxcl15","Cxcl16","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Cxcr1","Cxcr2","Cxcr3","Cxcr4","Ebi3","Elane","Erbb2","F2rl1","Fasl","Flt3","Flt3l","Foxp3","Gpi1","Hc","Ido1","Ifna1","Ifna2","Ifna4","Ifnar2","Ifnb1","Ifng","Ifngr1","Igf1r","Il10","Il11","Il11ra1","Il12a","Il12b","Il12rb1","Il12rb2","Il13","Il13ra1","Il13ra2","Il15","Il16","Il17a","Il17b","Il17f","Il17ra","Il17rb","Il18","Il18rap","Il19","Il1a","Il1r1","Il1r2","Il1rap","Il1rapl2","Il1rl1","Il1rl2","Il1rn","Il2","Il21","Il21r","Il22","Il22ra2","Il23a","Il23r","Il24","Il25","Il27","Il2rb","Il2rg","Il34","Il3ra","Il4","Il4ra","Il5","Il5ra","Il6","Il6ra","Il6st","Il7","Il7r","Il9","Irak1","Irak2","Irak3","Irak4","Irf3","Itk","Jak1","Jak2","Jak3","Kit","Lif","Lta","Ltb","Lyn","Mif","Mme","Ms4a2","Myd88","Nfatc2","Nfkb1","Nod2","Osm","Pin1","Pparg","Prkce","Ptprc","Rel","S100b","Saa1","Sh2b2","Sigirr","Socs1","Spp1","Stat1","Stat4","Stat6","Tgfb1","Tgfb2","Ticam2","Tnf","Tnfrsf11a","Tnfrsf1a","Tnfrsf4","Tnfsf10","Tnfsf11","Tnfsf12","Tnfsf13","Tnfsf13b","Tnfsf14","Tnfsf15","Tnfsf18","Tnfsf4","Tnfsf8","Traf2","Traf3","Traf6","Txk","Vegfa","Xcl1","Xcr1"),
  "Cytotoxicity"=c("Cd1d2","Cd8a","Ctsh","Fcgr1","Fcgr3","Gzmb","Gzmk","Gzmm","H2-D1","H2-K1","H2-M3","H2-T23","H60a","Il12a","Il21","Il23a","Il7r","Klrb1c","Klrk1","Lag3","Prf1","Ptprc","Pvr","Pvrl2","Sh2d1a","Sh2d1b1","Stat5b","Tap1","Ulbp1","Xcl1"),
  "DC_function"=c("Ccl19","Ccl5","Ccr1","Ccr2","Ccr5","Cd40","Cd40lg","Cd83","Cd86","Cr2","Cxcr1","Cxcr4","Il10","Lyn","Rag1","Relb","Tgfb1"),
  "Humoral_immune_rersponse"=c("Aire","Blnk","Bmi1","Bst1","Bst2","C4b","C5ar1","Ccl12","Ccl2","Ccl22","Ccl3","Ccl7","Ccr2","Ccr6","Ccr7","Cd28","Cd37","Cd40","Cd83","Crp","Cxcl13","Ebi3","Fcer2a","Foxj1","Gpi1","Gpr183","Ifnb1","Ifng","Il10","Il1b","Il6","Il7","Itgb2","Lta","Ltf","Ly86","Ly96","Mbl2","Mef2c","Mnx1","Ms4a1","Ms4a2","Nfkb1","Nod2","Pax5","Pdcd1","Pou2af1","Pou2f2","Psmb10","Sh2d1a","St6gal1","Tfe3","Tfeb","Tnf","Trem1","Trem2","Ythdf2"),
  "Inflammation"=c("Adora2a","Anxa1","Axl","Bcl6","C1qbp","C3","C3ar1","C4b","Camp","Ccl1","Ccl11","Ccl12","Ccl17","Ccl19","Ccl2","Ccl20","Ccl22","Ccl24","Ccl25","Ccl26","Ccl3","Ccl4","Ccl5","Ccl7","Ccl8","Ccr1","Ccr2","Ccr3","Ccr4","Ccr7","Ccrl2","Cd14","Cd163","Cd180","Cd276","Cd28","Cd40","Cd40lg","Cd47","Cd97","Cebpb","Cklf","Clec7a","Cma1","Creb5","Crp","Csf1","Cspg4","Cxcl1","Cxcl10","Cxcl11","Cxcl13","Cxcl14","Cxcl15","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Cxcr1","Cxcr2","Cxcr4","Elane","F2rl1","Fas","Fasl","Fcer1a","Fcer2a","Fcgr2b","Fos","Foxp3","Fpr2","Hc","Hck","Hmgb1","Ido1","Il10","Il13","Il17a","Il17b","Il17f","Il18","Il1a","Il1b","Il1r1","Il1rap","Il1rl1","Il1rn","Il22","Il23r","Il24","Il25","Il27","Il2ra","Il34","Il4","Il4ra","Il5ra","Il6","Il6st","Il9","Itgb2","Jam3","Klrb1","Lbp","Lta","Ly86","Ly96","Lyn","Mapkapk2","Mefv","Mif","Ms4a2","Myd88","Nfatc4","Nfkb1","Nlrp3","Nos2","Nt5e","Pik3cd","Pik3cg","Pnma1","Pparg","Ptgs2","Pycard","Ripk2","S100a8","Sbno2","Sele","Stat3","Tgfb1","Thbs1","Ticam1","Ticam2","Tirap","Tlr1","Tlr2","Tlr3","Tlr4","Tlr5","Tlr6","Tlr7","Tlr8","Tlr9","Tnf","Tnfaip3","Tnfrsf1a","Tnfrsf1b","Tnfrsf4","Tnfsf4","Tollip","Xcl1","Xcr1"),
  "Innate_immune_response"=c("A2m","Abca1","Abcg1","Abl1","Aire","Angpt1","App","Atf1","Atf2","Atg12","Atg5","Axl","Bcl10","Bcl2","Bcl2l1","Bid","Bst2","Btk","C1qa","C1qb","C1ra","C1s1","C2","C3","C4b","C6","C7","C8a","C8b","C8g","C9","Camp","Card9","Casp1","Casp8","Ccl17","Ccl2","Ccl5","Ccr1","Ccr3","Ccr6","Ccr9","Cd14","Cd180","Cd1d1","Cd1d2","Cd36","Cd4","Cd46","Cd55","Cd74","Cd8a","Cd97","Cdk1","Cebpb","Cfb","Cfd","Cfh","Cfi","Cfp","Chuk","Clec4a2","Clec4n","Clec5a","Clec7a","Clu","Colec12","Cr2","Creb1","Crebbp","Crp","Csf1","Csf1r","Csf2","Ctss","Cxcl10","Cxcl11","Cxcl16","Cxcl2","Cxcl9","Cxcr2","Cxcr3","Cxcr6","Cybb","Cyfip2","Cyld","Ddx58","Defb1","Dmbt1","Dusp4","Dusp6","Ecsit","Elk1","Ep300","F12","F2rl1","Fadd","Fcgr1","Fos","Gzmk","Gzmm","Hamp","Hc","Hck","Hmgb1","Ifih1","Ifitm1","Ifitm2","Ifna1","Ifna2","Ifnb1","Ifngr1","Ikbkb","Ikbke","Ikbkg","Il18r1","Il18rap","Il1a","Il1b","Il1r1","Il1rap","Il1rapl2","Il1rl1","Il1rl2","Il23a","Il23r","Il27","Il34","Il4","Irak1","Irak2","Irak3","Irak4","Irf3","Irf7","Irgm2","Isg15","Isg20","Itch","Itga5","Itgam","Itgax","Jak1","Jak2","Jak3","Klrg1","Lbp","Lcn2","Lgals3","Lilra5","Ly86","Ly96","Lyn","Map2k1","Map2k2","Map2k4","Map3k1","Map3k5","Map3k7","Map4k2","Mapk1","Mapk11","Mapk14","Mapk3","Mapk8","Mapkapk2","Marco","Masp1","Masp2","Mavs","Mbl2","Mefv","Mif","Mst1r","Mx1","Mx2","Myd88","Ncf4","Nfkb1","Nfkb2","Nfkbia","Nlrc5","Nlrp3","Nod1","Nod2","Pik3cd","Pik3cg","Pin1","Pparg","Pycard","Rela","Ripk2","S100b","Saa1","Serping1","Sigirr","Slamf7","Stat1","Syk","Tab1","Tank","Tbk1","Ticam1","Ticam2","Tirap","Tlr1","Tlr2","Tlr3","Tlr4","Tlr5","Tlr6","Tlr7","Tlr8","Tlr9","Tmem173","Tnfaip3","Tollip","Traf2","Traf3","Traf6","Txnip","Tyk2","Ubc","Xcl1","Zbp1"),
  "Interferon_response"=c("Ccr7","Cd3e","Ciita","Cxcl16","Ddx58","Eomes","Fadd","Gbp5","H2-Aa","H2-Ab1","H60a","Ifi27","Ifi35","Ifi44","Ifi44l","Ifih1","Ifit1","Ifit2","Ifit3","Ifitm1","Ifitm2","Ifna1","Ifna2","Ifna4","Ifnar1","Ifnar2","Ifngr1","Ifnl2","Il12rb2","Irf7","Irf8","Irgm2","Mavs","Nlrc5","Nos2","Runx3","Sh2d1b1","Tbk1","Tmem173","Ulbp1"),
  "Interleukins"=c("A2m","Bcl10","Card11","Card9","Casp1","Ccl2","Ccr2","Ccr7","Cd1d1","Cd1d2","Cd276","Cd28","Cd34","Cd3e","Cd40lg","Cd83","Cma1","Cmklr1","Csf2","Cxcl15","Cxcr1","Cxcr2","Ebi3","Egr1","Elane","Fadd","Fcer1a","Fcer1g","Fcgr2b","Foxj1","Foxp3","Gfi1","H2-Eb1","H2-Q1","Icosl","Ido1","Ifng","Il10","Il11","Il11ra1","Il12a","Il12b","Il12rb1","Il13","Il13ra1","Il15","Il16","Il17a","Il17f","Il17ra","Il18","Il18rap","Il19","Il1a","Il1b","Il1r1","Il1r2","Il1rap","Il1rapl2","Il1rl1","Il1rl2","Il1rn","Il2","Il21","Il22","Il22ra2","Il23a","Il23r","Il24","Il25","Il27","Il2rb","Il2rg","Il3","Il4","Il5","Il5ra","Il6","Il6ra","Il6st","Il7","Il9","Inpp5d","Irak2","Irak3","Irak4","Irf1","Irf4","Irf8","Itk","Jak2","Jak3","Lag3","Ltb","Mapk3","Mapkapk2","Mavs","Myd88","Nfkb1","Nlrp3","Nod1","Nod2","Pycard","Rel","Rela","Sele","Stat5b","Syk","Tcf7","Ticam1","Ticam2","Tigit","Tirap","Tlr1","Tlr3","Tlr4","Tlr5","Tlr6","Tlr7","Tlr8","Tlr9","Tnfrsf11a","Tnfrsf1a","Tnfsf4","Tollip","Traf2","Traf6","Txk"),
  "Leukocyte_function"=c("Ccl19","Ccl25","Ccl4","Ccr1","Ccr7","Cd34","Cklf","Clec7a","Cxcl10","Cxcl12","Cxcl2","Cxcl5","Elane","F2rl1","Foxj1","Fut7","Hc","Icam1","Il16","Il23r","Il3","Itga2","Itga2b","Itga4","Itga5","Itga6","Itgal","Itgam","Itgb1","Itgb2","Itgb3","Jam3","Lbp","Pecam1","Psen1","Psen2","Sele","Selplg","Syk","Tlr2","Vcam1"),
  "Macrophage_function"=c("C3ar1","Casp8","Ccl2","Ccl5","Ccr2","Ccr5","Ccr7","Cd1d1","Cd69","Cklf","Cmklr1","Crp","Csf1","Csf1r","Csf2","Cx3cl1","Cx3cr1","Eng","Fcer1a","Fcer2a","H60a","Hc","Il13","Il17f","Il18","Il1b","Il1rl1","Il23a","Il34","Il4","Il4ra","Itgb3","Lbp","Msr1","Nfkbia","Pparg","Prkce","Rora","Saa1","Sbno2","Syk","Thbs1","Tlr1","Ulbp1","Vegfa"),
  "Mast_cell_function"=c("C5ar1","Cd48","Fcer1a","Tpsab1"),
  "Mature_B_function"=c("Atm","Bcl10","Bcl6","Blk","Blnk","Bmi1","Card11","Casp3","Cd19","Cd200r1","Cd22","Cd37","Cd74","Cd79b","Cxcl13","Fas","Fcgr2b","Gpr183","H2-Q10","Ikbkb","Ikbkg","Ikzf1","Il13","Il13ra1","Il1r2","Il21","Il2rg","Il5","Il7r","Inpp5d","Itga2","Jak3","Lck","Lyn","Mapk1","Mef2c","Mif","Nt5e","Prdm1","Sh2b2","Stat5b","Syk","Ticam1","Tirap","Tnfaip3","Tnfrsf4","Tnfsf4"),
  "Mature_T_function"=c("Anxa1","Bcl10","Card11","Ccl19","Ccl2","Ccl3","Ccr2","Ccr7","Cd1d2","Cd247","Cd274","Cd276","Cd48","Cd5","Cd59b","Cd80","Cd83","Cd86","Ctsh","Cxcl13","Eomes","Fasl","Fcgr4","Foxj1","Fut7","H2-Aa","H2-Ab1","H2-D1","H2-DMa","H2-K1","H2-M3","H2-Q1","H2-T23","Icam1","Ido1","Ifnar1","Ifnb1","Ikzf1","Il12a","Il12rb1","Il23a","Il2rg","Il6st","Il7r","Itgal","Itgam","Itgax","Itgb2","Itk","Lcp1","Mapk1","Mill2","Psen1","Psen2","Psmb10","Pvr","Pvrl2","Rag1","Rorc","Rps6","Spn","Stat5b","Syk","Tap1","Tcf7","Tgfb1","Tgfb2","Tigit","Tnfsf18","Tnfsf8","Traf2","Txk","Vcam1","Xcl1","Zap70"),
  "MHC1_MHC2"=c("Cd160","Cd1d1","Cd40lg","Cd74","Ciita","Ctsh","Fcer1g","Fcgr1","Fcgr2b","Fcgr3","H2-Aa","H2-Ab1","H2-D1","H2-DMa","H2-DMb1","H2-DMb2","H2-Eb1","H2-K1","H2-M3","H2-Ob","H2-Q1","H2-Q10","H2-Q2","H2-T23","Il10","Klrk1","Lag3","Mr1","Nlrc5","Pml","Tap1","Tap2","Tapbp"),
  "Microglial_function"=c("Casp1","Cx3cr1","Nod2","Tlr6","Tlr7"),
  "NK_function"=c("Axl","Ccl2","Ccl3","Ccl4","Ccl5","Ccl7","Cd2","Cd244","Cd247","Cd7","Cd96","Flt3l","H2-M3","H60a","Ifi27","Ikzf1","Il11ra1","Il12a","Il12b","Il12rb1","Il15","Il15ra","Il21","Il23a","Il2rb","Itgb2","Klra1","Klra15","Klra17","Klra2","Klra20","Klra21","Klra27","Klra3","Klra4","Klra5","Klra6","Klra7","Klrb1c","Klrc1","Klrd1","Klrk1","Lag3","Mertk","Mill2","Ncam1","Pvr","Pvrl2","Sh2d1a","Sh2d1b1","Slamf7","Stat5b","Ulbp1"),
  "Pathogen_response"=c("Ccl2","Cd14","Fos","Hmgb1","Ifnb1","Ifng","Il10","Il12a","Il1b","Il6","Irak1","Irf3","Jun","Lta","Ly86","Ly96","Nfkbia","Ptgs2","Rela","Ripk2","Tbk1","Ticam1","Tlr3","Tlr4","Tlr6","Tlr7","Tlr8","Tnf","Tnfrsf1a"),
  "Phagocytosis"=c("Anxa1","C3","Cd14","Cd36","Cd44","Cd47","Clec7a","Colec12","Crp","Csf1","Csf2","Fas","Fcer1g","Fcgr1","Fcgr2b","Fcgr3","Ifng","Il1rl1","Itgam","Itgb2","Marco","Mbl2","Mfge8","Mif","Myd88","Nod1","Pecam1","Siglec1","Ticam1","Tlr3","Tlr9","Tnf","Tnfsf11"),
  "Regulation_inflammatory_response"=c("Bcl6","C3","C3ar1","C4b","Ccl1","Ccl12","Ccl17","Ccl19","Ccl2","Ccl20","Ccl22","Ccl24","Ccl25","Ccl3","Ccl4","Ccl7","Ccl8","Ccr2","Ccr3","Ccr4","Ccr7","Cd14","Cd40","Cd40lg","Cebpb","Crp","Csf1","Cxcl1","Cxcl10","Cxcl11","Cxcl2","Cxcl3","Cxcl5","Cxcl9","Cxcr1","Cxcr2","Cxcr4","Fos","Il10","Il1r1","Il1rap","Il22","Il9","Itgb2","Ly96","Myd88","Nfkb1","Nos2","Ripk2","Sele","Tirap","Tlr1","Tlr2","Tlr3","Tlr4","Tlr5","Tlr6","Tlr7","Tlr9","Tollip"),
  "Treg_function"=c("Ccr6","Foxp3","Ikzf2","Il9","Irf4","Irf8","Pou2f2","Rel","Tnfsf11"),
  "Senescence"=c("Abl1","Atm","Bmi1","Cdkn1a","Egr1","Ets1","Hras","Ifng","Igf1r","Irf3","Irf5","Irf7","Map2k1","Mapk14","Myc","Nfkb1","Plau","Prkcd","Serpinb2","Tgfb1","Trp53","Twist1"),
  "T_cell_function"=c("Ada","Anxa1","Bcl10","Bcl2","Bcl6","Btla","Card11","Casp3","Ccl11","Ccl19","Ccl2","Ccl3","Ccl5","Ccl7","Ccnd3","Ccr2","Ccr3","Ccr4","Ccr5","Ccr6","Ccr7","Cd1d1","Cd1d2","Cd2","Cd247","Cd27","Cd274","Cd276","Cd28","Cd3d","Cd3e","Cd3g","Cd4","Cd40","Cd40lg","Cd47","Cd48","Cd5","Cd59b","Cd74","Cd80","Cd83","Cd86","Cd8a","Cd8b1","Cebpb","Cma1","Crebbp","Csf2","Ctla4","Ctsh","Cxcl12","Cxcl13","Cxcr3","Cxcr4","Dpp4","Egr1","Eomes","Fas","Fasl","Fcgr4","Flt3","Foxj1","Foxp3","Fut7","Gata3","Gfi1","Gpr44","Gzmb","H2-Aa","H2-Ab1","H2-D1","H2-DMa","H2-K1","H2-M3","H2-Q1","H2-T23","Havcr2","Icam1","Icos","Icosl","Ido1","Ifnar1","Ifnb1","Ifng","Ikzf1","Ikzf2","Il10","Il12a","Il12b","Il12rb1","Il12rb2","Il13","Il13ra1","Il15","Il17a","Il18","Il18r1","Il18rap","Il1b","Il1r1","Il2","Il21","Il23a","Il25","Il27","Il2ra","Il2rg","Il4","Il4ra","Il5","Il6","Il6st","Il7","Il7r","Il9","Irf1","Irf4","Irf8","Itch","Itga1","Itgal","Itgam","Itgax","Itgb2","Itk","Jak1","Jak2","Jak3","Lag3","Lck","Lcp1","Lgals3","Maf","Map3k7","Mapk1","Mapk8","Mill2","Nfatc1","Nfatc2","Nfkb1","Nos2","Pdcd1","Pdcd1lg2","Pou2f2","Psen1","Psen2","Psmb10","Ptprc","Pvr","Pvrl2","Rag1","Rel","Relb","Ripk2","Rora","Rorc","Rps6","Sell","Socs1","Socs3","Spn","Spp1","Stat1","Stat4","Stat5b","Stat6","Syk","Tap1","Tbx21","Tcf7","Tgfb1","Tgfb2","Tgfb3","Thy1","Tigit","Tlr4","Tlr6","Tmed1","Tnf","Tnfrsf13c","Tnfrsf14","Tnfrsf4","Tnfrsf8","Tnfsf11","Tnfsf13b","Tnfsf14","Tnfsf18","Tnfsf4","Tnfsf8","Traf2","Traf6","Trp53","Txk","Tyk2","Vcam1","Vegfa","Xcl1","Yy1","Zap70"),
  "T_cell_anergy"=c("Cma1","Ctla4","Gzmb","Icos","Itch","Itga1","Jak1","Lgals3","Pdcd1","Sell"),
  "T_cell_differentiation"=c("Ada","Bcl2","Cd1d1","Cd3d","Cd4","Cd74","Egr1","Flt3","Il27","Il7","Irf4","Nos2","Socs1"),
  "T_cell_proliferation"=c("Casp3","Ccnd3","Cd3e","Cxcl12","Cxcr4","Icosl","Il15","Il1b","Pdcd1lg2","Ptprc","Ripk2","Spp1","Tnfrsf13c","Tnfsf13b","Tnfsf14","Traf6","Trp53"),
  "T_cell_regulation"=c("Btla","Cd27","Cd3g","Cd47","Cd8a","Cd8b1","Dpp4","Fas","Foxp3","Icam1","Il2","Lag3","Lck","Map3k7","Tgfb1","Thy1","Tnfrsf14"),
  "Th1_Th2_differentiation"=c("Cd2","Cd28","Cd40","Cd40lg","Ifng","Il12b","Il18","Il2ra","Il4","Il6","Relb"),
  "Th1_function"=c("Ccr5","Csf2","Cxcr3","Havcr2","Il12rb2","Il18r1","Il18rap","Il2","Il27","Irf1","Nfkb1","Socs1","Stat1","Stat4","Tbx21","Tlr4","Tlr6","Tnf","Vegfa"),
  "Th2_function"=c("Bcl6","Ccl11","Ccl5","Ccl7","Ccr3","Ccr4","Cebpb","Gata3","Gfi1","Gpr44","Icos","Il10","Il13","Il13ra1","Il1r1","Il25","Il4ra","Il5","Jak1","Jak3","Maf","Nfatc1","Nfatc2","Stat6","Tmed1"),
  "Th17_function"=c("Il17a","Il21","Rora","Rorc"),
  "TLR"=c("Cd86","Gfi1","Irak1","Irak2","Irf3","Irf4","Map3k7","Mapkapk2","Myd88","Nfkbia","Prkce","Tbk1","Ticam1","Ticam2","Tirap","Tlr1","Tlr2","Tlr3","Tlr4","Tlr5","Tlr6","Tlr7","Tlr8","Tlr9","Traf3","Traf6"),
  "TNF_superfamily"=c("Cd27","Cd40","Cd40lg","Cd70","Fas","Fasl","Lta","Ltb","Ltbr","Tnf","Tnfrsf10b","Tnfrsf11a","Tnfrsf11b","Tnfrsf12a","Tnfrsf13b","Tnfrsf13c","Tnfrsf14","Tnfrsf17","Tnfrsf18","Tnfrsf1a","Tnfrsf1b","Tnfrsf4","Tnfrsf8","Tnfrsf9","Tnfsf10","Tnfsf11","Tnfsf12","Tnfsf13","Tnfsf13b","Tnfsf14","Tnfsf15","Tnfsf18","Tnfsf4","Tnfsf8"),
  "TNF_superfamily_members"=c("Cd40lg","Cd70","Fasl","Lta","Ltb","Tnf","Tnfsf10","Tnfsf11","Tnfsf12","Tnfsf13","Tnfsf13b","Tnfsf14","Tnfsf15","Tnfsf18","Tnfsf4","Tnfsf8"),
  "TNF_superfamily_receptors"=c("Cd27","Cd40","Fas","Ltbr","Tnfrsf10b","Tnfrsf11a","Tnfrsf11b","Tnfrsf12a","Tnfrsf13b","Tnfrsf13c","Tnfrsf14","Tnfrsf17","Tnfrsf18","Tnfrsf1a","Tnfrsf1b","Tnfrsf4","Tnfrsf8","Tnfrsf9"),
  "Transcription_factors"=c("Atf1","Atf2","Batf","Bcl10","Bcl2","Bcl6","Card11","Card9","Cd34","Cd36","Cd40","Cdh1","Cebpb","Ciita","Clu","Cmklr1","Creb1","Creb5","Crebbp","Cyld","Ddx58","Ecsit","Egr1","Egr2","Elk1","Eomes","Ep300","Ets1","Fos","Foxj1","Foxp3","Gata3","Gfi1","Gtf3c1","Hck","Hmgb1","Icam1","Ikbkb","Ikbkg","Ikzf1","Il10","Il4","Il5","Irak1","Irak2","Irak3","Irf1","Irf2","Irf3","Irf4","Irf5","Irf7","Irf8","Itch","Itgb2","Jun","Maf","Mapk1","Mapk3","Mavs","Mef2c","Mnx1","Myc","Nfatc2","Nfatc3","Nfatc4","Nfkb1","Nfkb2","Nfkbia","Nlrc5","Nlrp3","Nod1","Nod2","Pax5","Pou2f2","Pparg","Pycard","Rel","Rela","Relb","Rora","Rorc","Runx1","Runx3","Sigirr","Smad4","Snai1","Stat1","Stat3","Stat4","Stat6","Tal1","Tbx21","Tcf7","Tfe3","Tfeb","Ticam1","Tlr3","Tlr9","Tmem173","Tnfrsf11a","Tnfrsf4","Tnfrsf8","Tnfsf18","Traf2","Traf3","Traf6","Twist1","Xbp1","Yy1"),
  "Transmembrane_transporter"=c("Abca1","Abcb1a","Abcg1","Akt3","Ambp","Amica1","Apoe","App","Atg10","Atg16l1","Atg7","Axl","Bax","C8g","Ccl3","Ccl4","Ccl5","Ccr1","Ccr5","Cd36","Cd3g","Cmah","Col1a1","Csf2","Cxcl1","Cxcl12","Cxcr4","Cybb","Dmbt1","Fez1","Fyn","Icam1","Ifit1","Igf2r","Il13","Il1b","Il4","Itgb3","Lbp","Lcn2","Lcp1","Lrp1","Ltf","Lyn","Lyve1","Lyz2","Map2k1","Mapk14","Mertk","Msr1","Nefl","Nfatc1","Nup107","Pparg","Prkcd","Prkce","Psen1","Slc11a1","Slc7a11","Syk","Syt17","Tap1","Tap2","Tfrc","Tmed1"),
  "Transporter_function"=c("Abca1","Abcb1a","Abcg1","Akt3","Ambp","Amica1","Anxa1","Apoe","App","Atg10","Atg12","Atg16l1","Atg5","Atg7","Axl","Bax","C3","C8g","Ccl3","Ccl4","Ccl5","Ccr1","Ccr5","Cd14","Cd36","Cd3g","Cd44","Cd47","Clec7a","Cmah","Col1a1","Colec12","Crp","Csf1","Csf2","Cxcl1","Cxcl12","Cxcr4","Cybb","Dmbt1","Fas","Fcer1g","Fcgr1","Fcgr2b","Fcgr3","Fez1","Fyn","Icam1","Ifit1","Ifng","Igf2r","Il13","Il1b","Il1rl1","Il4","Itgam","Itgb2","Itgb3","Lamp1","Lbp","Lcn2","Lcp1","Lrp1","Ltf","Lyn","Lyve1","Lyz2","Map2k1","Mapk14","Marco","Mbl2","Mertk","Mfge8","Mif","Msr1","Myd88","Nefl","Nfatc1","Nod1","Nup107","Pecam1","Pparg","Prkcd","Prkce","Psen1","Siglec1","Slc11a1","Slc7a11","Syk","Syt17","Tap1","Tap2","Tfrc","Ticam1","Tlr3","Tlr9","Tmed1","Tnf","Tnfsf11")
)

immune.integrated<-AddModuleScore(immune.integrated, 
                    features =geneset.list, 
                    name=c("Adhesion","Adaptive_immunity","Antigen_processing",
                           "Apoptosis","Autophagy","B_activation","B_differentiation",
                           "B_function","B_proliferation","B_regulation",
                           "Bacterial_response","Basic_cell_function",
                           "Cancer_progression","CD4_T_function","CD_molecules",
                           "Cell_cycle","Chemokine_receptors","Chemokines",
                           "Chemokines_chemokine_receptors","Chemotaxis",
                           "Complement_pathway","Cytokine_receptors","Cytokines",
                           "Cytokines_cytokine_receptors","Cytotoxicity","DC_function",
                           "Humoral_immune_rersponse","Inflammation",
                           "Innate_immune_response","Interferon_response","Interleukins",
                           "Leukocyte_function","Macrophage_function","Mast_cell_function",
                           "Mature_B_function","Mature_T_function","MHC1_MHC2",
                           "Microglial_function","NK_function","Pathogen_response",
                           "Phagocytosis","Regulation_inflammatory_response",
                           "Treg_function","Senescence","T_cell_function",
                           "T_cell_anergy","T_cell_differentiation","T_cell_proliferation",
                           "T_cell_regulation","Th1_Th2_differentiation","Th1_function",
                           "Th2_function","Th17_function","TLR","TNF_superfamily",
                           "TNF_superfamily_members","TNF_superfamily_receptors",
                           "Transcription_factors","Transmembrane_transporter",
                           "Transporter_function"),ctrl=80)
                    
VlnPlot(immune.integrated, features = c("Adhesion1","Adaptive_immunity2","Antigen_processing3",
                           "Apoptosis4","Autophagy5","B_activation6","B_differentiation7",
                           "B_function8","B_proliferation9","B_regulation10",
                           "Bacterial_response11","Basic_cell_function12",
                           "Cancer_progression13","CD4_T_function14","CD_molecules15",
                           "Cell_cycle16","Chemokine_receptors17","Chemokines18",
                           "Chemokines_chemokine_receptors19","Chemotaxis20",
                           "Complement_pathway21","Cytokine_receptors22","Cytokines23",
                           "Cytokines_cytokine_receptors24","Cytotoxicity25","DC_function26",
                           "Humoral_immune_rersponse27","Inflammation28",
                           "Innate_immune_response29","Interferon_response30","Interleukins31",
                           "Leukocyte_function32","Macrophage_function33","Mast_cell_function34",
                           "Mature_B_function35","Mature_T_function36","MHC1_MHC237",
                           "Microglial_function38","NK_function39","Pathogen_response40",
                           "Phagocytosis41","Regulation_inflammatory_response42",
                           "Treg_function43","Senescence44","T_cell_function45",
                           "T_cell_anergy46","T_cell_differentiation47","T_cell_proliferation48",
                           "T_cell_regulation49","Th1_Th2_differentiation50","Th1_function51",
                           "Th2_function52","Th17_function53","TLR54","TNF_superfamily55",
                           "TNF_superfamily_members56","TNF_superfamily_receptors57",
                           "Transcription_factors58","Transmembrane_transporter59",
                           "Transporter_function60"),split.by = "BTN",pt.size = 0)



######################################
###monocle3 pseudotime################
######################################
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)

DefaultAssay(immune.intergated)="integrated"
immune.momacro=subset(immune.intergated,idents=c("TIM2","ERM","TIM1","Monocyte"))
DimPlot(immune.momacro, reduction = "umap",label = T) 

#Re-clustering for higher resolution of Monocytes and Mo-TAMs
agg=immune.momacro
agg <- FindNeighbors(agg, reduction = "pca", dims = 1:20)
agg <- FindClusters(agg,resolution = 0.50)
immune.markers <- FindAllMarkers(agg, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.15)
write.csv(immune.markers, file="momacrodegres0.5.csv")
new.cluster.ids <- c("Bsg+","Pro-inflammatory (Cxcl9+)","Translation initiating", "Transcription active", "Antigen presenting","Monocytes", "Translation active")
names(new.cluster.ids) <- levels(agg)
agg <- RenameIdents(agg, new.cluster.ids)


# ...1 Convert to cell_data_set object ------------------------
gene_annotation <- as.data.frame(rownames(agg@reductions[["pca"]]@feature.loadings),
                                 row.names = rownames(agg@reductions[["pca"]]@feature.loadings))
colnames(gene_annotation) <- "gene_short_name"
cell_metadata=agg[[]] 
cell_metadata=data.frame("cell_id"=rownames(cell_metadata),cell_metadata)
expression_matrix=GetAssayData(object = agg, slot = "scale.data")
cds_from_seurat <- new_cell_data_set(expression_matrix,
                                     cell_metadata = cell_metadata, 
                                     gene_metadata = gene_annotation)
recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
names(recreate.partition) <- cds_from_seurat@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition
list_cluster <- agg@active.ident
names(list_cluster) <- agg@assays[["integrated"]]@data@Dimnames[[2]]
cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster
cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
cds_from_seurat@int_colData@listData$reducedDims@listData[["UMAP"]] <-agg@reductions[["umap"]]@cell.embeddings
cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = F)

# ...2. Plot the pseudotime ------------------------
plot_cells(cds_from_seurat, 
           color_cells_by = 'cluster',
           label_groups_by_cluster=TRUE,
           label_leaves=FALSE,
           label_branch_points=F,
           label_roots = TRUE,
           show_trajectory_graph=TRUE,
           graph_label_size=3, cell_size = 0.45, trajectory_graph_segment_size = 0.75)

cluster.before.trajectory <- plot_cells(cds_from_seurat,
                                        color_cells_by = 'cluster',
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +theme(legend.position = "right")

cluster.names <- plot_cells(cds_from_seurat,
                            color_cells_by = "cluster",
                            label_groups_by_cluster = FALSE,
                            group_label_size = 5) +
  scale_color_manual(values = c('red', 'blue', 'green', 'maroon', 'yellow', 'grey', 'cyan')) +
  theme(legend.position = "right")


# ...3. Learn trajectory graph ------------------------
cds_from_seurat <- learn_graph(cds_from_seurat, use_partition = FALSE)

plot_cells(cds_from_seurat,
           color_cells_by = 'cluster',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE,
           group_label_size = 5,cell_size = 0.65, trajectory_graph_segment_size = 0.75)


# ...4. Order the cells in pseudotime -------------------

cds_from_seurat <- order_cells(cds_from_seurat, reduction_method = 'UMAP', 
                               root_cells = colnames(cds_from_seurat[,clusters(cds_from_seurat) == "Monocyte"]))


View(cds_from_seurat)

plot_cells(cds_from_seurat,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = T,
           label_branch_points = FALSE,
           label_roots = TRUE,
           label_leaves = F,
           show_trajectory_graph=TRUE,
           genes=list("QK"=c("Qk")),
           trajectory_graph_color="Red")

plot_cells(cds_from_seurat,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = T,
           label_branch_points = FALSE,
           label_roots = T,
           label_leaves = F,
           show_trajectory_graph=TRUE,trajectory_graph_color="Red")


