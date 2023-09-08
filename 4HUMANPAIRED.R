
C1.data <- Read10X(data.dir ="C1")
C2.data <- Read10X(data.dir ="C2")
C3.data <- Read10X(data.dir ="C3")
C4.data <- Read10X(data.dir ="C4")


M1.data <- Read10X(data.dir ="M1")
M2.data <- Read10X(data.dir ="M2")
M3.data <- Read10X(data.dir ="M3")
M4.data <- Read10X(data.dir ="M4")


P1.data <- Read10X(data.dir ="P1")
P2.data <- Read10X(data.dir ="P2")
P3.data <- Read10X(data.dir ="P3")
P4.data <- Read10X(data.dir ="P4")



C1 <- CreateSeuratObject(counts = C1.data, project = "Core1", min.cells = 3)
C1$num <- "1"
C1$BTN <- "Core"
C1[["percent.mt"]] <- PercentageFeatureSet(C1, pattern = "^MT-")
FeatureScatter(C1,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
VlnPlot(C1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C1 <- subset(C1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

C2 <- CreateSeuratObject(counts = C2.data, project = "Core2", min.cells = 3)
C2$num <- "2"
C2$BTN <- "Core"
C2[["percent.mt"]] <- PercentageFeatureSet(C2, pattern = "^MT-")
FeatureScatter(C2,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
VlnPlot(C2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C2 <- subset(C2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

C3 <- CreateSeuratObject(counts = C3.data, project = "Core3", min.cells = 3)
C3$num <- "3"
C3$BTN <- "Core"
C3[["percent.mt"]] <- PercentageFeatureSet(C3, pattern = "^MT-")
FeatureScatter(C3,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
VlnPlot(C3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C3 <- subset(C3, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

C4 <- CreateSeuratObject(counts = C4.data, project = "Core4", min.cells = 3)
C4$num <- "4"
C4$BTN <- "Core"
C4[["percent.mt"]] <- PercentageFeatureSet(C4, pattern = "^MT-")
FeatureScatter(C4,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
VlnPlot(C4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
C4 <- subset(C4, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)



M1 <- CreateSeuratObject(counts = M1.data, project = "Middle1", min.cells = 3)
M1$num <- "1"
M1$BTN <- "Middle"
M1[["percent.mt"]] <- PercentageFeatureSet(M1, pattern = "^MT-")
FeatureScatter(M1,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
VlnPlot(M1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
M1 <- subset(M1, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

M2 <- CreateSeuratObject(counts = M2.data, project = "Middle2", min.cells = 3)
M2$num <- "2"
M2$BTN <- "Middle"
M2[["percent.mt"]] <- PercentageFeatureSet(M2, pattern = "^MT-")
FeatureScatter(M2,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
VlnPlot(M2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
M2 <- subset(M2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

M3 <- CreateSeuratObject(counts = M3.data, project = "Middle3", min.cells = 3)
M3$num <- "3"
M3$BTN <- "Middle"
M3[["percent.mt"]] <- PercentageFeatureSet(M3, pattern = "^MT-")
FeatureScatter(M3,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
VlnPlot(M3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
M3 <- subset(M3, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

M4 <- CreateSeuratObject(counts = M4.data, project = "Middle4", min.cells = 3)
M4$num <- "4"
M4$BTN <- "Middle"
M4[["percent.mt"]] <- PercentageFeatureSet(M4, pattern = "^MT-")
FeatureScatter(M4,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
VlnPlot(M4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
M4 <- subset(M4, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

P1 <- CreateSeuratObject(counts = P1.data, project = "Periphery1", min.cells = 3)
P1$num <- "1"
P1$BTN <- "Periphery"
P1[["percent.mt"]] <- PercentageFeatureSet(P1, pattern = "^MT-")
FeatureScatter(P1,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
VlnPlot(P1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
P1 <- subset(P1, ssubset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

P2 <- CreateSeuratObject(counts = P2.data, project = "Periphery2", min.cells = 3)
P2$num <- "2"
P2$BTN <- "Periphery"
P2[["percent.mt"]] <- PercentageFeatureSet(P2, pattern = "^MT-")
FeatureScatter(P2,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
VlnPlot(P2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
P2 <- subset(P2, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

P3 <- CreateSeuratObject(counts = P3.data, project = "Periphery3", min.cells = 3)
P3$num <- "3"
P3$BTN <- "Periphery"
P3[["percent.mt"]] <- PercentageFeatureSet(P3, pattern = "^MT-")
FeatureScatter(P3,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
VlnPlot(P3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
P3 <- subset(P3,subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

P4 <- CreateSeuratObject(counts = P4.data, project = "Periphery4", min.cells = 3)
P4$num <- "4"
P4$BTN <- "Periphery"
P4[["percent.mt"]] <- PercentageFeatureSet(P4, pattern = "^MT-")
FeatureScatter(P4,feature1 = "nCount_RNA",feature2 ="nFeature_RNA" ) + geom_smooth(method = "lm")
VlnPlot(P4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
P4 <- subset(P4, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 25)

###############INTERGTARE LARGE DATASET AS SEURAT##########

datasets <- list(C1,C2,C3,C4,M1,M2,M3,M4,P1,P2,P3,P4)
saveRDS(datasets, file = "HU4.rds")


HU4 <- lapply(X = datasets, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = HU4)
HU4 <- lapply(X = HU4, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

################Periphery as reference######################
anchors <- FindIntegrationAnchors(object.list = HU4, reference = c(9,10,11,12), reduction = "rpca",
                                  dims = 1:50)

HU4.integrated <- IntegrateData(anchorset = anchors, dims = 1:50,k.weight =75)

HU4.integrated <- ScaleData(HU4.integrated, verbose = FALSE)

#Determine appropriate number of PCA
HU4.integrated <- RunPCA(HU4.integrated, npcs = 30, verbose = T)
print(immune.combined[["pca"]], dims = 1:20, nfeatures = 5)
ElbowPlot(HU4.integrated)
pca.mouse.combined <- immune.combined[["pca"]]
View(pca.mouse.combined)
pca.mouse.combined@stdev
HU4.integrated <- RunUMAP(HU4.integrated, dims = 1:50)

#############Change resolution to modify number or clusters############
HU4.integrated <- FindNeighbors(HU4.integrated, reduction = "pca", dims = 1:20)
HU4.integrated <- FindClusters(HU4.integrated,resolution = 0.15)

#############Visulizing Dimplot##############
DimPlot(HU.IMMUNE, reduction = "umap",label = T) 

################
################For subsequent analysis, please refer to 5MICEPAIRED.R file.
