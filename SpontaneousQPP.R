#read in data
####ctrl & early
Ctrl1.data <- Read10X(data.dir ="57")
Ctrl2.data <- Read10X(data.dir ="58")
Ctrl3.data <- Read10X(data.dir ="61")
Ear1.data <- Read10X(data.dir ="59")
Ear2.data <- Read10X(data.dir ="60")
Ear3.data <- Read10X(data.dir ="62")

Ctrl1<- CreateSeuratObject(counts = Ctrl1.data, project = "Ctrl1", min.cells = 3)
Ctrl1$stim <- "Ctrl"
Ctrl1$gemgroup <- "C1"
Ctrl1[["percent.mt"]] <- PercentageFeatureSet(Ctrl1, pattern = "^mt-")
FeatureScatter(Ctrl1,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
Ctrl1 <- subset(Ctrl1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
VlnPlot(Ctrl1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


Ctrl2 <- CreateSeuratObject(counts = Ctrl2.data, project = "Ctrl2", min.cells = 3)
Ctrl2$stim <- "Ctrl"
Ctrl2$gemgroup <- "C2"
Ctrl2[["percent.mt"]] <- PercentageFeatureSet(Ctrl2, pattern = "^mt-")
FeatureScatter(Ctrl2,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
Ctrl2 <- subset(Ctrl2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
VlnPlot(Ctrl2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Ctrl3 <- CreateSeuratObject(counts = Ctrl3.data, project = "Ctrl3", min.cells = 3)
Ctrl3$stim <- "Ctrl"
Ctrl3$gemgroup <- "C3"
Ctrl3[["percent.mt"]] <- PercentageFeatureSet(Ctrl3, pattern = "^mt-")
FeatureScatter(Ctrl3,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
Ctrl3 <- subset(Ctrl3, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
VlnPlot(Ctrl3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



Ear1<- CreateSeuratObject(counts = Ear1.data, project = "Ear1", min.cells = 3)
Ear1$stim <- "Ear"
Ear1$gemgroup <- "E1"
Ear1[["percent.mt"]] <- PercentageFeatureSet(Ear1, pattern = "^mt-")
FeatureScatter(Ear1,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
Ear1 <- subset(Ear1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
VlnPlot(Ear1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


Ear2 <- CreateSeuratObject(counts = Ear2.data, project = "Ear2", min.cells = 3)
Ear2$stim <- "Ear"
Ear2$gemgroup <- "E2"
Ear2[["percent.mt"]] <- PercentageFeatureSet(Ear2, pattern = "^mt-")
FeatureScatter(Ear2,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
Ear2 <- subset(Ear2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
VlnPlot(Ear2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

Ear3 <- CreateSeuratObject(counts = Ear3.data, project = "Ear3", min.cells = 3)
Ear3$stim <- "Ear"
Ear3$gemgroup <- "E3"
Ear3[["percent.mt"]] <- PercentageFeatureSet(Ear3, pattern = "^mt-")
FeatureScatter(Ear3,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
Ear3 <- subset(Ear3, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
VlnPlot(Ear3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

####late
spoagg.data <- Read10X(data.dir ="SPO")
spo <- CreateSeuratObject(counts = spoagg.data, project = "SPONT_CELL", min.cells = 3)
spo[["percent.mt"]] <- PercentageFeatureSet(spo, pattern = "^mt-")
FeatureScatter(spo,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
VlnPlot(spo, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
spo <- subset(spo, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)
sgemgroup <- sapply(strsplit(rownames(spo@meta.data), split="-"), "[[", 2) 
spo <- AddMetaData(object=spo, metadata=data.frame(gemgroup=sgemgroup, row.names=rownames(spo@meta.data)))
spo$stim <- "Late"

#merge dataset
datasets <- list(spo,Ctrl1,Ctrl2,Ctrl3,Ear1,Ear2,Ear3)
datasets=lapply( datasets, function(x) SCTransform(x,vars.to.regress = "percent.mt", verbose = TRUE))

immune.anchors <- FindIntegrationAnchors(object.list = datasets, dims = 1:20)
immune.integrated <- IntegrateData(anchorset = immune.anchors, dims = 1:20, k.weight = 100)
immune.integrated <- ScaleData(immune.integrated, features = rownames(immune.integrated))
immune.integrated <- RunPCA(immune.integrated, npcs = 30, verbose = FALSE)
immune.integrated <- RunUMAP(immune.integrated, reduction = "pca", dims = 1:30)
immune.integrated <- FindNeighbors(immune.integrated, reduction = "pca", dims = 1:30)
DefaultAssay(immune.integrated)="integrated"
immune.integrated<- FindClusters(immune.integrated, resolution = 0.85)

#Find markers for every cluster compared to all remaining cells, report only the positive ones
immune.markers <- FindAllMarkers(immune.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

################
################For subsequent analysis, please refer to 5MICEPAIRED.R file.
