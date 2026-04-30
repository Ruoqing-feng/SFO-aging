# Ruoqing Feng

library(Seurat)
library(ggplot2)

old1 = read_monkey_data("/Volumes/SFO/The_SFO_project/data/raw/monkey/sc_old1/","old1","/Volumes/SFO/The_SFO_project/data/raw/monkey/sc_old1_meta.csv")
old2 = read_monkey_data("/Volumes/SFO/The_SFO_project/data/raw/monkey/sc_old2/","old2","/Volumes/SFO/The_SFO_project/data/raw/monkey/sc_old2_meta.csv")
young1 = read_monkey_data("/Volumes/SFO/The_SFO_project/data/raw/monkey/sc_young1/","young1","/Volumes/SFO/The_SFO_project/data/raw/monkey/sc_young1_meta.csv")
young2 = read_monkey_data("/Volumes/SFO/The_SFO_project/data/raw/monkey/sc_young2/","young2","/Volumes/SFO/The_SFO_project/data/raw/monkey/sc_young2_meta.csv")

old1 = subset(old1,nCount_RNA <= 1500 & nFeature_RNA <= 1500 & nFeature_RNA >= 100 & area <= 2000)
VlnPlot(old1,c("nCount_RNA","nFeature_RNA","area"),pt.size = 0)

old2 = subset(old2,nCount_RNA <= 2000 & nFeature_RNA <= 1500 & nFeature_RNA >= 100 & area <= 2000)
VlnPlot(old2,c("nCount_RNA","nFeature_RNA","area"),pt.size = 0)

young1 = subset(young1,nCount_RNA <= 1000 & nFeature_RNA <= 1000 & nFeature_RNA >= 100 & area <= 2000)
VlnPlot(young1,c("nCount_RNA","nFeature_RNA","area"),pt.size = 0)

young2 = subset(young2,nCount_RNA <= 3000 & nFeature_RNA <= 2500 & nFeature_RNA >= 100 & area <= 2500)
VlnPlot(young2,c("nCount_RNA","nFeature_RNA","area"),pt.size = 0)

all = merge(old1,c(young1,old2,young2))
all = SCTransform(all,assay = "RNA")
all <- RunPCA(all, features = VariableFeatures(all), npcs = 50)
ElbowPlot(all,ndims = 50)

all <- FindNeighbors(all, dims = 1:30)
all <- RunUMAP(all, dims = 1:30)
DimPlot(all,label = T,split.by = "orig.ident",ncol = 2)

all <- FindClusters(all, resolution = 0.1, algorithm = 4)
all <- FindClusters(all, resolution = 0.2, algorithm = 4)
all <- FindClusters(all, resolution = 0.3, algorithm = 4)

DimPlot(all,label = T)
all = PrepSCTFindMarkers(all)
markers = FindAllMarkers(all,min.pct = 0.1,only.pos = T)

saveRDS(all,file = "/Volumes/SFO/The_SFO_project/data/processed/monkey_all4_cellbin_clean_clustered.rds")
FeaturePlot(all,"CTSB")
VlnPlot(all,"CX3CR1",pt.size = 0)
intersect(VariableFeatures(all),c("TYROBP"))

#######
## BIN 50

data = Read10X("/Volumes/SFO/monkey/round1/old/analysis/out_10x/")
colnames(data) = paste0("old1_",colnames(data))
meta = read.csv("/Volumes/SFO/monkey/round1/old/analysis/old_50_meta.csv",row.names = 1)
rownames(meta) = paste0("old1_",rownames(meta))

old = CreateSeuratObject(counts = data,meta.data = meta,)
old@meta.data$orig.ident = "old1"

old = subset(old,nCount_RNA <= 3500 & nFeature_RNA >= 100 & nFeature_RNA <= 2000)
VlnPlot(old,c("nCount_RNA","nFeature_RNA"),pt.size = 0)

data = Read10X("/Volumes/SFO/monkey/round1/young/analysis/out_10x/")
colnames(data) = paste0("young1_",colnames(data))
meta = read.csv("/Volumes/SFO/monkey/round1/young/analysis/young_50_meta.csv",row.names = 1)
rownames(meta) = paste0("young1_",rownames(meta))

young = CreateSeuratObject(counts = data,meta.data = meta)
young@meta.data$orig.ident = "young1"

young = subset(young,nCount_RNA <= 2300 & nFeature_RNA >= 100 & nFeature_RNA <= 1500)
VlnPlot(young,c("nCount_RNA","nFeature_RNA"),pt.size = 0)

data = Read10X("/Volumes/SFO/monkey/round2/old/analysis/out_10x/")
colnames(data) = paste0("old2_",colnames(data))
meta = read.csv("/Volumes/SFO/monkey/round2/old/analysis/old_50_meta.csv",row.names = 1)
rownames(meta) = paste0("old2_",rownames(meta))

old2 = CreateSeuratObject(counts = data,meta.data = meta,min.features = 30)
old2@meta.data$orig.ident = "old2"

old2 = subset(old2,nCount_RNA <= 4000 & nFeature_RNA >= 300 & nFeature_RNA <= 2300)
VlnPlot(old2,c("nCount_RNA","nFeature_RNA"),pt.size = 0)

data = Read10X("/Volumes/SFO/monkey/round2/young/analysis/out_10x/")
colnames(data) = paste0("young2_",colnames(data))
meta = read.csv("/Volumes/SFO/monkey/round2/young/analysis/young_50_meta.csv",row.names = 1)
rownames(meta) = paste0("young2_",rownames(meta))

young2 = CreateSeuratObject(counts = data,meta.data = meta)
young2@meta.data$orig.ident = "young2"

young2 = subset(young2,nCount_RNA <= 6000 & nFeature_RNA >= 300 & nFeature_RNA <= 3500)
VlnPlot(young2,c("nCount_RNA","nFeature_RNA"),pt.size = 0)

all = merge(old,c(young,old2,young2))
Idents(all) = "orig.ident"
VlnPlot(all,c("nCount_RNA","nFeature_RNA"),pt.size = 0)

all = SCTransform(all,assay = "RNA")
all <- RunPCA(all, features = VariableFeatures(all), npcs = 50)
ElbowPlot(all,ndims = 50)

all <- FindNeighbors(all, dims = 1:30)
all <- RunUMAP(all, dims = 1:30)

all <- FindClusters(all, resolution = 0.3, algorithm = 4)
DimPlot(all,label = T)

saveRDS(all,file = "/Volumes/SFO/The_SFO_project/data/processed/monkey_all4_bin50_clean_clustered.rds")

keep.genes <- rowSums(old2@assays$RNA$counts > 0) >= 50
old2 = old2[names(keep.genes),]
old2 = SCTransform(old2,assay = "RNA")
old2 <- RunPCA(old2, features = VariableFeatures(old2), npcs = 50)
old2 <- FindNeighbors(old2, dims = 1:30)
old2 <- RunUMAP(old2, dims = 1:30)

old2 <- FindClusters(old2, resolution = 0.5, algorithm = 4)

DimPlot(old2,label = T)

ggplot(old2@meta.data,aes(x = x, y = y ,color = SCT_snn_res.0.5))+
  geom_point()+
  theme_classic()

keep.genes <- rowSums(old@assays$RNA$counts > 0) >= 50
old = old[names(keep.genes),]
old = SCTransform(old,assay = "RNA")
old <- RunPCA(old, features = VariableFeatures(old), npcs = 50)
old <- FindNeighbors(old, dims = 1:30)
old <- RunUMAP(old, dims = 1:30)

old <- FindClusters(old, resolution = 0.5, algorithm = 4)

DimPlot(old,label = T)

ggplot(old@meta.data,aes(x = x, y = y ,color = SCT_snn_res.0.5))+
  geom_point()+
  theme_classic()




keep.genes <- rowSums(young2@assays$RNA$counts > 0) >= 50
young2 = young2[names(keep.genes),]
young2 = SCTransform(young2,assay = "RNA")
young2 <- RunPCA(young2, features = VariableFeatures(young2), npcs = 50)
young2 <- FindNeighbors(young2, dims = 1:30)
young2 <- RunUMAP(young2, dims = 1:30)

young2 <- FindClusters(young2, resolution = 0.3, algorithm = 4)

DimPlot(young2,label = T)

ggplot(young2@meta.data,aes(x = x, y = y ,color = SCT_snn_res.0.3))+
  geom_point()+
  theme_classic()

keep.genes <- rowSums(young@assays$RNA$counts > 0) >= 50
young = young[names(keep.genes),]
young = SCTransform(young,assay = "RNA")
young <- RunPCA(young, features = VariableFeatures(young), npcs = 50)
young <- FindNeighbors(young, dims = 1:30)
young <- RunUMAP(young, dims = 1:30)

young <- FindClusters(young, resolution = 0.5, algorithm = 4)

DimPlot(young,label = T)

ggplot(young@meta.data,aes(x = x, y = y ,color = SCT_snn_res.0.5))+
  geom_point()+
  theme_classic()

saveRDS(young,"/Volumes/SFO/The_SFO_project/data/processed/monkey_young1_bin50.rds")
saveRDS(young2,"/Volumes/SFO/The_SFO_project/data/processed/monkey_young2_bin50.rds")
saveRDS(old,"/Volumes/SFO/The_SFO_project/data/processed/monkey_old1_bin50.rds")
saveRDS(old2,"/Volumes/SFO/The_SFO_project/data/processed/monkey_old2_bin50.rds")


young1.clean <- subset(young1, subset = !is.na(x) & !is.na(y))
young2.clean <- subset(young2, subset = !is.na(x) & !is.na(y))
old2.clean <- subset(old2, subset = !is.na(x) & !is.na(y))
old1.clean <- subset(old1, subset = !is.na(x) & !is.na(y))

library(SpatialExperiment)
library(Banksy)

gcm <- tmp@assays$SCT$data[VariableFeatures(tmp),]
locs <- as.matrix(tmp@meta.data[,c("x","y")])
mdata = tmp@meta.data

tmp <- SpatialExperiment(assay = list(counts = gcm),spatialCoords = locs,colData = mdata,sample_id = "old2")

lambda <- c(0.2,0.5,0.8)
k_geom <- 8
npcs <- 20
aname <- "counts"
tmp <- Banksy::computeBanksy(tmp, assay_name = aname, k_geom = k_geom,compute_agf = T)

set.seed(1000)
tmp <- Banksy::runBanksyPCA(tmp, lambda = lambda, npcs = npcs,use_agf = T)

set.seed(1000)
tmp <- Banksy::clusterBanksy(tmp, lambda = lambda, npcs = npcs, resolution = c(0.1,0.3,0.5,0.8), use_agf = TRUE)

ggplot(colData(tmp),aes(x = x, y = y,colour = clust_M1_lam0.2_k50_res0.8))+
  geom_point()+
  theme_classic()

old2.clean.banksy = run_stereoeq_banksy(scp_clean = old2.clean,sample_name = "old2",lam = c(0.2,0.5,0.8),res = c(0.1,0.3,0.5,0.8,1),npc = 20,k_geom = 8,aname = "counts")
young2.clean.banksy = run_stereoeq_banksy(scp_clean = young2.clean,sample_name = "young2",lam = c(0.2,0.5,0.8),res = c(0.1,0.3,0.5,0.8,1),npc = 20,k_geom = 8,aname = "counts")
old1.clean.banksy = run_stereoeq_banksy(scp_clean = old1.clean,sample_name = "old1",lam = c(0.2,0.5,0.8),res = c(0.1,0.3,0.5,0.8,1),npc = 20,k_geom = 8,aname = "counts")
young1.clean.banksy = run_stereoeq_banksy(scp_clean = young1.clean,sample_name = "young1",lam = c(0.2,0.5,0.8),res = c(0.1,0.3,0.5,0.8,1),npc = 20,k_geom = 8,aname = "counts")


colData(young2.clean.banksy)$region = as.numeric(as.character(colData(young2.clean.banksy)$clust_M1_lam0.2_k50_res0.5))
colData(young2.clean.banksy)$region[colData(young2.clean.banksy)$region == 5] = "Ependyma"
colData(young2.clean.banksy)$region[colData(young2.clean.banksy)$region == 10] = "SFO"
colData(young2.clean.banksy)$region[colData(young2.clean.banksy)$region %in% c(2,6,7,9)] = "Thalamus"
colData(young2.clean.banksy)$region[colData(young2.clean.banksy)$region == 8] = "ChP"
colData(young2.clean.banksy)$region[colData(young2.clean.banksy)$region %in% c(3,4)] = "Fimbra"
colData(young2.clean.banksy)$region[colData(young2.clean.banksy)$region %in% c(1,11)] = "FWM"

p1 = ggplot(colData(young2.clean.banksy)[colData(young2.clean.banksy)$region == "Fimbra",],aes(x = x, y = y,colour = region))+
  geom_point()+
  theme_classic()

colData(young2.clean.banksy)$cell = rownames(colData(young2.clean.banksy))
colData(young2.clean.banksy)$region[colData(young2.clean.banksy)$cell %in% CellSelector(p1)] = "Thalamus"

p1 = ggplot(colData(young2.clean.banksy),aes(x = x, y = y,colour = region))+
  geom_point()+
  theme_classic()

s = colData(young2.clean.banksy)[colData(young2.clean.banksy)$cell %in% CellSelector(p1),]
group <- as.factor(s$region)
cols <- rainbow(length(levels(group)))
plot(s$x, s$y, pch=16, cex=0.3, asp=1,col=cols[group])
poly <- locator(type = "l")
inside <- sp::point.in.polygon(s$x, s$y, poly$x, poly$y) > 0
selected_cells <- s$cell[inside]
colData(young2.clean.banksy)$region[colData(young2.clean.banksy)$cell %in% selected_cells] = "SFO"

saveRDS(young2.clean.banksy,"/Volumes/SFO/The_SFO_project/data/processed/monkey_young2_banksy_anno.rds")

colData(old2.clean.banksy)$region = as.numeric(as.character(colData(old2.clean.banksy)$clust_M1_lam0.2_k50_res1))
colData(old2.clean.banksy)$region[colData(old2.clean.banksy)$region %in% c(4,9,12)] = "ChP"
colData(old2.clean.banksy)$region[colData(old2.clean.banksy)$region %in% c(6,11)] = "Ependyma"
colData(old2.clean.banksy)$region[colData(old2.clean.banksy)$region %in% c(8)] = "Thalamus"

colData(old2.clean.banksy)$cell = rownames(colData(old2.clean.banksy))
p2 = ggplot(colData(old2.clean.banksy)[colData(old2.clean.banksy)$clust_M1_lam0.2_k50_res1 %in% c(5,7,10,1,3,2),],aes(x = x, y = y,colour = clust_M1_lam0.2_k50_res1))+
  geom_point()+
  theme_classic()
colData(old2.clean.banksy)$region[colData(old2.clean.banksy)$cell %in% CellSelector(p2)] = "Fimbra"
p2 = ggplot(colData(old2.clean.banksy)[colData(old2.clean.banksy)$region %in% c("ChP"),],aes(x = x, y = y,colour = region))+
  geom_point()+
  theme_classic()
colData(old2.clean.banksy)$region[colData(old2.clean.banksy)$cell %in% CellSelector(p2)] = "Fimbra"
p2 = ggplot(colData(old2.clean.banksy)[colData(old2.clean.banksy)$region %in% c("Ependyma"),],aes(x = x, y = y,colour = region))+
  geom_point()+
  theme_classic()
colData(old2.clean.banksy)$region[colData(old2.clean.banksy)$cell %in% CellSelector(p2)] = "Fimbra"

s = colData(old2.clean.banksy)[colData(old2.clean.banksy)$region %in% c("Ependyma"),]
plot(s$x, s$y, pch=16, cex=0.3, asp=1)
poly <- locator(type = "l")
inside <- sp::point.in.polygon(s$x, s$y, poly$x, poly$y) > 0
selected_cells <- s$cell[inside]
colData(old2.clean.banksy)$region[colData(old2.clean.banksy)$cell %in% selected_cells] = "Fimbra"

s = colData(old2.clean.banksy)[colData(old2.clean.banksy)$region %in% c("Fimbra"),]
plot(s$x, s$y, pch=16, cex=0.3, asp=1)
poly <- locator(type = "l")
inside <- sp::point.in.polygon(s$x, s$y, poly$x, poly$y) > 0
selected_cells <- s$cell[inside]
colData(old2.clean.banksy)$region[colData(old2.clean.banksy)$cell %in% selected_cells] = "FWM"

p1 = ggplot(colData(old2.clean.banksy),aes(x = x, y = y,colour = region))+
  geom_point()+
  theme_classic()
s = colData(old2.clean.banksy)[colData(old2.clean.banksy)$cell %in% CellSelector(p1),]
group <- as.factor(s$region)
cols <- rainbow(length(levels(group)))
plot(s$x, s$y, pch=16, cex=0.3, asp=1,col=cols[group])
poly <- locator(type = "l")
inside <- sp::point.in.polygon(s$x, s$y, poly$x, poly$y) > 0
selected_cells <- s$cell[inside]
colData(old2.clean.banksy)$region[colData(old2.clean.banksy)$cell %in% selected_cells] = "SFO"

colData(old2.clean.banksy)$region[colData(old2.clean.banksy)$region %in% c(10,3)] = "Thalamus"

p1 = ggplot(colData(old2.clean.banksy),aes(x = x, y = y,colour = region))+
  geom_point()+
  theme_classic()

s = colData(old2.clean.banksy)[colData(old2.clean.banksy)$cell %in% CellSelector(p1),]
group <- as.factor(s$region)
cols <- rainbow(length(levels(group)))
plot(s$x, s$y, pch=16, cex=0.3, asp=1,col=cols[group])
poly <- locator(type = "l")
inside <- sp::point.in.polygon(s$x, s$y, poly$x, poly$y) > 0
selected_cells <- s$cell[inside]
colData(old2.clean.banksy)$region[colData(old2.clean.banksy)$cell %in% selected_cells] = "ChP"

saveRDS(old2.clean.banksy,"/Volumes/SFO/The_SFO_project/data/processed/monkey_old2_banksy_anno.rds")

colData(young1.clean.banksy)$region = as.numeric(as.character(colData(young1.clean.banksy)$clust_M1_lam0.2_k50_res0.3))
colData(young1.clean.banksy)$region[colData(young1.clean.banksy)$region %in% c(3,9)] = "Vantricular"
colData(young1.clean.banksy)$cell = rownames(colData(young1.clean.banksy))

s = colData(young1.clean.banksy)[colData(young1.clean.banksy)$clust_M1_lam0.2_k50_res0.3 %in% c(1,2,4,5,6,7,8),]
plot(s$x, s$y, pch=16, cex=0.3, asp=1)
poly <- locator(type = "l")
inside <- sp::point.in.polygon(s$x, s$y, poly$x, poly$y) > 0
selected_cells <- s$cell[inside]
colData(young1.clean.banksy)$region[colData(young1.clean.banksy)$cell %in% selected_cells] = "Thalamus"

s = colData(young1.clean.banksy)[colData(young1.clean.banksy)$clust_M1_lam0.2_k50_res0.3 %in% c(1,2,3),]
group <- as.factor(s$region)
cols <- rainbow(length(levels(group)))
plot(s$x, s$y, pch=16, cex=0.3, asp=1,col=cols[group])
poly <- locator(type = "l")
inside <- sp::point.in.polygon(s$x, s$y, poly$x, poly$y) > 0
selected_cells <- s$cell[inside]
colData(young1.clean.banksy)$region[colData(young1.clean.banksy)$cell %in% selected_cells] = "SFO"

colData(young1.clean.banksy)$region[colData(young1.clean.banksy)$region %in% c(1,4)] = "Fimbra"

s = colData(young1.clean.banksy)[colData(young1.clean.banksy)$region %in% c(2),]
plot(s$x, s$y, pch=16, cex=0.3, asp=1)
poly <- locator(type = "l")
inside <- sp::point.in.polygon(s$x, s$y, poly$x, poly$y) > 0
selected_cells <- s$cell[inside]
colData(young1.clean.banksy)$region[colData(young1.clean.banksy)$cell %in% selected_cells] = "Vantricular"

colData(young1.clean.banksy)$region[colData(young1.clean.banksy)$region %in% c(2)] = "Fimbra"

s = colData(young1.clean.banksy)[colData(young1.clean.banksy)$region != "SFO",]
group <- as.factor(s$region)
cols <- rainbow(length(levels(group)))
plot(s$x, s$y, pch=16, cex=0.3, asp=1,col=cols[group])
poly <- locator(type = "l")
inside <- sp::point.in.polygon(s$x, s$y, poly$x, poly$y) > 0
selected_cells <- s$cell[inside]
colData(young1.clean.banksy)$region[colData(young1.clean.banksy)$cell %in% selected_cells] = "SFO"

colData(young1.clean.banksy)$region[colData(young1.clean.banksy)$region %in% c(5,6,7,8)] = "Fimbra"

s = colData(young1.clean.banksy)[colData(young1.clean.banksy)$region %in% c("Vantricular","Fimbra"),]
group <- as.factor(s$region)
cols <- rainbow(length(levels(group)))
plot(s$x, s$y, pch=16, cex=0.3, asp=1,col=cols[group])
poly <- locator(type = "l")
inside <- sp::point.in.polygon(s$x, s$y, poly$x, poly$y) > 0
selected_cells <- s$cell[inside]
colData(young1.clean.banksy)$region[colData(young1.clean.banksy)$cell %in% selected_cells] = "Vantricular"

p1 = ggplot(colData(young1.clean.banksy),aes(x = x, y = y,colour = region))+
  geom_point()+
  theme_classic()

s = colData(young1.clean.banksy)[colData(young1.clean.banksy)$cell %in% CellSelector(p1),]
group <- as.factor(s$region)
cols <- rainbow(length(levels(group)))
plot(s$x, s$y, pch=16, cex=0.3, asp=1,col=cols[group])
poly <- locator(type = "l")
inside <- sp::point.in.polygon(s$x, s$y, poly$x, poly$y) > 0
selected_cells <- s$cell[inside]
colData(young1.clean.banksy)$region[colData(young1.clean.banksy)$cell %in% selected_cells] = "SFO"

saveRDS(young1.clean.banksy,"/Volumes/SFO/The_SFO_project/data/processed/monkey_young1_banksy_anno.rds")

colData(old1.clean.banksy)$region = as.numeric(as.character(colData(old1.clean.banksy)$clust_M1_lam0.2_k50_res0.5))
colData(old1.clean.banksy)$cell = rownames(colData(old1.clean.banksy))

colData(old1.clean.banksy)$region[colData(old1.clean.banksy)$region %in% c(4)] = "Vantricular"

s = colData(old1.clean.banksy)[colData(old1.clean.banksy)$region %in% c(1,8),]
group <- as.factor(s$region)
cols <- rainbow(length(levels(group)))
plot(s$x, s$y, pch=16, cex=0.3, asp=1,col=cols[group])
poly <- locator(type = "l")
inside <- sp::point.in.polygon(s$x, s$y, poly$x, poly$y) > 0
selected_cells <- s$cell[inside]
colData(old1.clean.banksy)$region[colData(old1.clean.banksy)$cell %in% selected_cells] = "SFO"

s = colData(old1.clean.banksy)[colData(old1.clean.banksy)$region %in% c(1,2,5,7,8,9,11),]
group <- as.factor(s$region)
cols <- rainbow(length(levels(group)))
plot(s$x, s$y, pch=16, cex=0.3, asp=1,col=cols[group])
poly <- locator(type = "l")
inside <- sp::point.in.polygon(s$x, s$y, poly$x, poly$y) > 0
selected_cells <- s$cell[inside]
colData(old1.clean.banksy)$region[colData(old1.clean.banksy)$cell %in% selected_cells] = "Fimbra"

s = colData(old1.clean.banksy)[colData(old1.clean.banksy)$region %in% c(1,2,5,7,8,9,11,6),]
group <- as.factor(s$region)
cols <- rainbow(length(levels(group)))
plot(s$x, s$y, pch=16, cex=0.3, asp=1,col=cols[group])
poly <- locator(type = "l")
inside <- sp::point.in.polygon(s$x, s$y, poly$x, poly$y) > 0
selected_cells <- s$cell[inside]
colData(old1.clean.banksy)$region[colData(old1.clean.banksy)$cell %in% selected_cells] = "FWM"

s = colData(old1.clean.banksy)[colData(old1.clean.banksy)$region %in% c(1,2,5,7,8,9,11,6),]
group <- as.factor(s$region)
cols <- rainbow(length(levels(group)))
plot(s$x, s$y, pch=16, cex=0.3, asp=1,col=cols[group])
poly <- locator(type = "l")
inside <- sp::point.in.polygon(s$x, s$y, poly$x, poly$y) > 0
selected_cells <- s$cell[inside]
colData(old1.clean.banksy)$region[colData(old1.clean.banksy)$cell %in% selected_cells] = "Vantricular"

colData(old1.clean.banksy)$region[colData(old1.clean.banksy)$region %in% c(1,2,3,6,10,11,7,6,9,5,8)] = "Thalamus"

p2 = ggplot(colData(old1.clean.banksy)[colData(old1.clean.banksy)$region %in% c("Fimbra","Thalamus"),],aes(x = x, y = y,colour = region))+
  geom_point()+
  theme_classic()
old1.clean.banksy = readRDS("/Volumes/SFO/The_SFO_project/data/processed/monkey_old1_banksy_anno.rds")

saveRDS(old1.clean.banksy,"/Volumes/SFO/The_SFO_project/data/processed/monkey_old1_banksy_anno.rds")


monkey.all.bin50.seurat = readRDS("/Volumes/SFO/The_SFO_project/data/processed/monkey_all4_bin50_clean_clustered.rds")
clean.cells = c(rownames(colData(young1.clean.banksy)),rownames(colData(old1.clean.banksy)),
                rownames(colData(young2.clean.banksy)),rownames(colData(old2.clean.banksy)))
monkey.all.bin50.seurat@meta.data$cell = rownames(monkey.all.bin50.seurat@meta.data)
monkey.all.bin50.seurat.clean = subset(monkey.all.bin50.seurat,cell %in% clean.cells)


tmp = rbind(colData(young1.clean.banksy)[,c("sample_id","region")],colData(old1.clean.banksy)[,c("sample_id","region")])
tmp = rbind(tmp,colData(young2.clean.banksy)[,c("sample_id","region")])
tmp = rbind(tmp,colData(old2.clean.banksy)[,c("sample_id","region")])
tmp = tmp[rownames(monkey.all.bin50.seurat.clean@meta.data),]

monkey.all.bin50.seurat.clean@meta.data$region = tmp$region
monkey.all.bin50.seurat.clean@meta.data$region_big = monkey.all.bin50.seurat.clean@meta.data$region
monkey.all.bin50.seurat.clean@meta.data$region_big[monkey.all.bin50.seurat.clean@meta.data$region_big %in% c("ChP","Ependyma")] = "Vantricular"
Idents(monkey.all.bin50.seurat.clean) = "region_big"

DimPlot(monkey.all.bin50.seurat.clean,label = T,raster=FALSE)

#####adjustment
saveRDS(monkey.all.bin50.seurat.clean,"/Volumes/SFO/The_SFO_project/data/processed/monkey_all4_bin50_clean_clustered_clean.rds")
monkey.all.bin50.seurat.clean = readRDS("/Volumes/SFO/The_SFO_project/data/processed/monkey_all4_bin50_clean_clustered_clean.rds")
monkey.all.bin50.seurat.clean@meta.data$cell2 <- sub("^[^_]+_", "", monkey.all.bin50.seurat.clean@meta.data$cell)
write.csv(monkey.all.bin50.seurat.clean@meta.data,"/Volumes/SFO/The_SFO_project/data/external/monkey_meta.csv",quote = F)

monkey.all.bin50.seurat.clean = PrepSCTFindMarkers(monkey.all.bin50.seurat.clean)
region.markers = FindAllMarkers(monkey.all.bin50.seurat.clean,min.pct = 0.1,only.pos = T)


monkey.all.bin50.seurat.clean@meta.data$new = paste0(monkey.all.bin50.seurat.clean@meta.data$region_big,"_",monkey.all.bin50.seurat.clean@meta.data$orig.ident)
Idents(monkey.all.bin50.seurat.clean) = "new"

tmp = AggregateExpression(monkey.all.bin50.seurat.clean, assays = "RNA", return.seurat = T)
Idents(tmp)

library(DESeq2)
monkey.deg = data.frame()
for (i in c("SFO","Fimbra")){
    a.sample = paste0(i,c("-old1","-old2"))
    b.sample = paste0(i,c("-young1","-young2"))
    
    tmp.mtx = tmp@assays$RNA$counts[,c(a.sample,b.sample)]
    tmp.mtx = as.matrix(tmp.mtx)
    
    #keep.gene = rowSums(tmp.mtx > 0) >= 10
    tmp.mtx <- tmp.mtx[rowSums(tmp.mtx) >= 10, , drop = FALSE]
    
    sampleTable <- data.frame(
      sample = colnames(tmp.mtx),
      condition = c(rep("old",length(a.sample)),rep("young",length(b.sample)))  
    )
    
    dds <- DESeqDataSetFromMatrix(
      countData = tmp.mtx,
      colData = sampleTable,
      design = ~ condition
    )
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("condition", "old", "young"))
    tes.tab = as.data.frame(res)
    
    tes.tab$region = i
    tes.tab$gene = rownames(tes.tab)
    monkey.deg = rbind(monkey.deg,tes.tab)
}
monkey.deg$padj[is.na(monkey.deg$padj)] = 1

monkey.deg$sig = "none"
monkey.deg$sig[monkey.deg$log2FoldChange >= 0.5 & monkey.deg$pvalue <= 0.05] = "up"
monkey.deg$sig[monkey.deg$log2FoldChange <= (-0.5) & monkey.deg$pvalue <= 0.05] = "down"

saveRDS(monkey.deg,"/Volumes/SFO/The_SFO_project/data/external/monkey_deg_deseq2.rds")

monkey.deg = readRDS("/Volumes/SFO/The_SFO_project/data/external/monkey_deg_deseq2.rds")
p1 = ggplot(monkey.deg[monkey.deg$region == "SFO",],aes(x = log2FoldChange, y = -log10(pvalue),color = sig))+
  geom_point()+
  theme_classic()

monkey.all.bin50.seurat.clean@meta.data$new2 = monkey.all.bin50.seurat.clean@meta.data$orig.ident
monkey.all.bin50.seurat.clean@meta.data$new2[monkey.all.bin50.seurat.clean@meta.data$new2 %in% c("old1","old2")] = "old"
monkey.all.bin50.seurat.clean@meta.data$new2[monkey.all.bin50.seurat.clean@meta.data$new2 %in% c("young1","young2")] = "young"
monkey.all.bin50.seurat.clean@meta.data$new2 = paste0(monkey.all.bin50.seurat.clean@meta.data$region_big,"_",monkey.all.bin50.seurat.clean@meta.data$new2)
Idents(monkey.all.bin50.seurat.clean) = "new2"
tmp1 = FindMarkers(monkey.all.bin50.seurat.clean,ident.1 = "SFO_old",ident.2 = "SFO_young")
tmp1$gene = rownames(tmp1)
tmp1$region = "SFO"

tmp2 = FindMarkers(monkey.all.bin50.seurat.clean,ident.1 = "Fimbra_old",ident.2 = "Fimbra_young")
tmp2$gene = rownames(tmp2)
tmp2$region = "Fimbra"

monkey.deg.sc = rbind(tmp1,tmp2)
monkey.deg.sc$p_val_adj[is.na(monkey.deg.sc$p_val_adj)] = 1
monkey.deg.sc$sig = "none"
monkey.deg.sc$sig[monkey.deg.sc$avg_log2FC >= 0.5 & monkey.deg.sc$p_val_adj <= 0.05] = "up"
monkey.deg.sc$sig[monkey.deg.sc$avg_log2FC <= (-0.5) & monkey.deg.sc$p_val_adj <= 0.05] = "down"

saveRDS(monkey.deg.sc,"/Volumes/SFO/The_SFO_project/data/external/monkey_deg_sc.rds")
write.csv(monkey.deg.sc,"/Volumes/SFO/The_SFO_project/data/external/monkey_deg_sc.csv",quote = F)


monkey.deg.sc = readRDS("/Volumes/SFO/The_SFO_project/data/external/monkey_deg_sc.rds")

monkey.deg.sc$sig[monkey.deg.sc$region == "Fimbra"] = "non"
monkey.deg.sc$sig[monkey.deg.sc$avg_log2FC >= 0.8 & monkey.deg.sc$p_val_adj <= 0.05 & monkey.deg.sc$region == "Fimbra" ] = "up"
monkey.deg.sc$sig[monkey.deg.sc$avg_log2FC <= (-0.5) & monkey.deg.sc$p_val_adj <= 0.05 & monkey.deg.sc$region == "Fimbra"] = "down"

monkey.deg.sc$label = ""
sfo.label = c("LYZ","C1QC","CD163","SPP1","CD74","SERPINA3","B2M","PEG10","ULK4","CHGB","SLC1A2","RPL39","DLD")
fimbra.label = c("SPP1","C1QB","CD68","CD74","KLK6","MS4A7","PDE8B","SLC1A2","CST3","FTH1","SERPINA3")
monkey.deg.sc[monkey.deg.sc$region == "Fimbra" & monkey.deg.sc$gene %in% fimbra.label,]
for (i in sfo.label){
  monkey.deg.sc$label[monkey.deg.sc$gene == i & monkey.deg.sc$region == "SFO"] = i
}
for (i in fimbra.label){
  monkey.deg.sc$label[monkey.deg.sc$gene == i & monkey.deg.sc$region == "Fimbra"] = i
}

res_sfo_deg <- deg_resample_SFO(
  seu = monkey.all.bin50.seurat.clean,
  region_name = "SFO",
  assay = "RNA",
  slot = "counts",
  n_per_sample = 80,
  n_repeat = 100,
  min_cells_gene = 10,
  pseudocount = 1,
  seed = 123
)

res_sfo_deg.deg = monkey.deg[monkey.deg$region == "SFO",]
dim(res_sfo_deg.deg)
dim(res_sfo_deg)

rownames(res_sfo_deg) = res_sfo_deg$gene

tmp = merge(res_sfo_deg.deg,res_sfo_deg,by = "gene")
rownames(tmp) = tmp$gene
sfo.label = c("LYZ","C1QC","CD163","SPP1","CD74","SERPINA3","B2M","PEG10","ULK4","CHGB","SLC1A2","RPL39","DLD")
tmp[tmp$gene %in% sfo.label,]

tmp = tmp[!duplicated(tmp),]
rownames(tmp) = tmp$gene
tmp$sig = "non"
tmp$sig[tmp$log2FoldChange >= 0.6 & tmp$q_val <= 0.05 ] = "up"
tmp$sig[tmp$log2FoldChange <= (-0.6) & tmp$q_val <= 0.05] = "down"

tmp$label = ""
sfo.label = c("LYZ","C1QC","CD163","SPP1","CD74","SERPINA3","B2M","PEG10","ULK4","CHGB","SLC1A2","RPL39","DLD")
for (i in sfo.label){
  tmp$label[tmp$gene == i ] = i
}

res_Fimbra_deg <- deg_resample_SFO(
  seu = monkey.all.bin50.seurat.clean,
  region_name = "Fimbra",
  assay = "RNA",
  slot = "counts",
  n_per_sample = 200,
  n_repeat = 100,
  min_cells_gene = 10,
  pseudocount = 1,
  seed = 123
)

res_Fimbra_deg.deg = monkey.deg[monkey.deg$region == "Fimbra",]
tmp2 = merge(res_Fimbra_deg.deg,res_Fimbra_deg,by = "gene")

tmp2 = tmp2[!duplicated(tmp2),]
rownames(tmp2) = tmp2$gene
tmp2$sig = "non"
tmp2$sig[tmp2$log2FoldChange >= 0.6 & tmp2$q_val <= 0.05 ] = "up"
tmp2$sig[tmp2$log2FoldChange <= (-0.6) & tmp2$q_val <= 0.05] = "down"

tmp2$label = ""
fimbra.label = c("SPP1","C1QB","CD68","CD74","KLK6","MS4A7","PDE8B","SLC1A2","CST3","FTH1","SERPINA3")
for (i in fimbra.label){
  tmp2$label[tmp2$gene == i ] = i
}

all.deg = rbind(tmp,tmp2)

saveRDS(all.deg,"/Volumes/SFO/The_SFO_project/data/external/monkey_deg_deseq2.rds")
write.csv(all.deg,"/Volumes/SFO/The_SFO_project/data/external/monkey_deg_deseq2.csv",quote = F)

enrich.sfo = read.delim("/Volumes/SFO/The_SFO_project/data/external/sfo.csv",header = T,sep = ";")
enrich.sfo$Log.q.value. = c(7.953,5.826,4.246,0.730,2.997,3.452,4.454,2.748,6.029)
enrich.sfo = enrich.sfo[order(enrich.sfo$Log.q.value.,decreasing = T),]
enrich.sfo$Description = factor(enrich.sfo$Description,levels = enrich.sfo$Description,ordered = T)

p1 = ggplot(enrich.sfo,aes(x = Description, y = Log.q.value.,fill = sig))+
  geom_bar(stat = "identity")+
  geom_text(aes(y =0.2, label = Description), hjust =0, vjust =0.5, size =4.5,angle = 90,color = "grey8") +
  ylab("-log10(Q)")+xlab("")+
  theme_classic()+
  scale_fill_manual(values = c("#b30000","#5D90BA"),breaks = c("up","down"))+
  theme(axis.text.x = element_blank())

enrich.fimbra = read.delim("/Volumes/SFO/The_SFO_project/data/external/fimbria.csv",header = T,sep = ";")
enrich.fimbra$Log.q.value. = c(11.502,5.372,6.412,4.628,5.973,4.287,6.775,3.347)
enrich.fimbra = enrich.fimbra[order(enrich.fimbra$Log.q.value.,decreasing = T),]
enrich.fimbra$Description = factor(enrich.fimbra$Description,levels = enrich.fimbra$Description,ordered = T)
p2 = ggplot(enrich.fimbra,aes(x = Description, y = Log.q.value.,fill = sig))+
  geom_bar(stat = "identity")+
  geom_text(aes(y =0.2, label = Description), hjust =0, vjust =0.5, size =4.5,angle = 90,color = "grey8") +
  ylab("-log10(Q)")+xlab("")+
  theme_classic()+
  scale_fill_manual(values = c("#b30000","#5D90BA"),breaks = c("up","down"))+
  theme(axis.text.x = element_blank())
p = p1+ p2 
ggsave("/Users/ruoqing/Projects/sfo/fig_sn_merfish/monkey/enrichemnt.pdf",p,width = 10,height = 5)

region.color = c("#a14462","#356d67","#b0d45d","#ffe788","#9e6c69")
names(region.color) = c("SFO","Fimbra","Thalamus","Vantricular","FWM" )

for (i in unique(monkey.all.bin50.seurat.clean@meta.data$orig.ident)){
  print(i)
  p1 = ggplot(monkey.all.bin50.seurat.clean@meta.data[monkey.all.bin50.seurat.clean@meta.data$orig.ident == i,],aes(x = x, y = y,color = region_big))+
    geom_point(shape = 15,size = 0.5)+
    scale_color_manual(breaks = names(region.color),values = region.color)+
    theme_classic()+
    theme(legend.position = 'none',axis.title = element_blank())
  ggsave(paste0("/Users/ruoqing/Projects/sfo/fig_sn_merfish/monkey/",i,"_banksy.png"),p1,width = 6,height = 6)
}

library(tiff)
img <- readTIFF("/Volumes/SFO/monkey/round1/old/image/D03262C611_HE_regist.tif",native = FALSE)

img_h <- dim(img)[1]
img_w <- dim(img)[2]



p1 = plot_umap_legend(region.color,shape = 15)
ggsave("/Users/ruoqing/Projects/sfo/fig_sn_merfish/monkey/banksy_legend.pdf",p1)

tmp_clean = old2.clean
tmp_clean = young2.clean
tmp_clean = old1.clean
tmp_clean = young1.clean

tmp_clean = monkey.all.bin50.seurat.clean
target_gene <- "SERPINA3"
gene_expr <- tmp_clean@assays$SCT$data[target_gene, ]
plot_data = tmp_clean@meta.data
plot_data$Expression <- as.numeric(gene_expr)
p2 = ggplot(plot_data[plot_data$orig.ident %in% c("young2","old2"),], aes(x = x, y = y, color = Expression)) +
  geom_point(size = 0.5,shape = 15) +
  facet_grid(. ~ orig.ident)+
  scale_color_distiller(palette = "YlOrRd", direction = 1) +
  theme_classic()+
  coord_fixed() +xlab("")+ylab("")+
  ggtitle(target_gene)
ggsave("/Users/ruoqing/Projects/sfo/fig_sn_merfish/monkey/SERPINA3_show.png",p2,width = 10,height = 4,dpi = 300)
ggsave("/Users/ruoqing/Projects/sfo/fig_sn_merfish/monkey/SERPINA3_show.pdf",p2,width = 10,height = 4,dpi = 300)

target_gene <- "CD74"
gene_expr <- tmp_clean@assays$SCT$data[target_gene, ]
plot_data = tmp_clean@meta.data
plot_data$Expression <- as.numeric(gene_expr)
p3 = ggplot(plot_data[plot_data$orig.ident %in% c("young2","old2"),], aes(x = x, y = y, color = Expression)) +
  geom_point(size = 0.5,shape = 15) +
  facet_grid(. ~ orig.ident)+
  scale_color_distiller(palette = "YlOrRd", direction = 1) +
  theme_classic()+
  coord_fixed() +xlab("")+ylab("")+
  ggtitle(target_gene)
ggsave("/Users/ruoqing/Projects/sfo/fig_sn_merfish/monkey/cd74_show.png",p3,width = 10,height = 4,dpi = 300)
ggsave("/Users/ruoqing/Projects/sfo/fig_sn_merfish/monkey/cd74_show.pdf",p3,width = 10,height = 4,dpi = 300)

target_gene <- "C1QC"
gene_expr <- tmp_clean@assays$SCT$data[target_gene, ]
plot_data = tmp_clean@meta.data
plot_data$Expression <- as.numeric(gene_expr)
p4 = ggplot(plot_data[plot_data$orig.ident %in% c("young2","old2"),], aes(x = x, y = y, color = Expression)) +
  geom_point(size = 0.5,shape = 15) +
  facet_grid(. ~ orig.ident)+
  scale_color_distiller(palette = "YlOrRd", direction = 1) +
  theme_classic()+
  coord_fixed() +xlab("")+ylab("")+
  ggtitle(target_gene)
ggsave("/Users/ruoqing/Projects/sfo/fig_sn_merfish/monkey/c1qc_show.png",p4,width = 10,height = 4,dpi = 300)
ggsave("/Users/ruoqing/Projects/sfo/fig_sn_merfish/monkey/c1qc_show.pdf",p4,width = 10,height = 4,dpi = 300)

young2.clean@meta.data$clust_M1_lam0.2_k50_res0.5 = colData(young2.clean.banksy)$clust_M1_lam0.2_k50_res0.5
Idents(young2.clean) = "clust_M1_lam0.2_k50_res0.5"

markers = FindAllMarkers(young2.clean,min.pct = 0.1,only.pos = T)


monkey.cellbin.all = readRDS("/Volumes/SFO/The_SFO_project/data/processed/monkey_all4_cellbin_clean_clustered.rds")
Idents(monkey.cellbin.all) = "SCT_snn_res.0.2"
DimPlot(monkey.cellbin.all,label = T)
monkey.cellbin.all.markers = FindAllMarkers(monkey.cellbin.all,min.pct = 0.1,only.pos = T)

monkey.cellbin.all@meta.data$cluster = "Oligo"
monkey.cellbin.all@meta.data$cluster[monkey.cellbin.all@meta.data$SCT_snn_res.0.2 %in% c(2,14)] = "CPEC"
monkey.cellbin.all@meta.data$cluster[monkey.cellbin.all@meta.data$SCT_snn_res.0.2 %in% c(4,5,11,15)] = "Neuron"
monkey.cellbin.all@meta.data$cluster[monkey.cellbin.all@meta.data$SCT_snn_res.0.2 %in% c(6)] = "Ependymal"
monkey.cellbin.all@meta.data$cluster[monkey.cellbin.all@meta.data$SCT_snn_res.0.2 %in% c(1,8,10)] = "Astro"

Idents(monkey.cellbin.all) = "cluster"
library(harmony)
monkey.cellbin.oligo = subset(monkey.cellbin.all,cluster == "Oligo")

monkey.cellbin.oligo = SCTransform(monkey.cellbin.oligo,assay = "RNA",variable.features.n = 2000)
monkey.cellbin.oligo <- RunPCA(monkey.cellbin.oligo, features = VariableFeatures(monkey.cellbin.oligo), npcs = 50)
monkey.cellbin.oligo <- monkey.cellbin.oligo %>% 
  RunHarmony("orig.ident", plot_convergence = TRUE, dims.use = 1:30,reduction = "pca",reduction.save = "harmony", assay.use = "SCT")
monkey.cellbin.oligo <- FindNeighbors(monkey.cellbin.oligo, dims = 1:30,reduction = "harmony")
monkey.cellbin.oligo <- RunUMAP(monkey.cellbin.oligo, dims = 1:30,reduction = "harmony")
Idents(monkey.cellbin.oligo) = "orig.ident"
DimPlot(monkey.cellbin.oligo,label = T)

monkey.cellbin.oligo <- FindClusters(monkey.cellbin.oligo, resolution = 0.3)

monkey.cellbin.oligo.clean = subset(monkey.cellbin.oligo,SCT_snn_res.0.2 %in% c(0,6,7))
monkey.cellbin.oligo.clean <- RunUMAP(monkey.cellbin.oligo.clean, dims = 1:30,reduction = "harmony")

monkey.cellbin.oligo.clean <- FindClusters(monkey.cellbin.oligo.clean, resolution = 0.2)

DimPlot(monkey.cellbin.oligo.clean,label = T)
FeaturePlot(monkey.cellbin.oligo.clean,"SERPINA3")
monkey.cellbin.oligo.marker = FindAllMarkers(monkey.cellbin.oligo.clean,min.pct = 0.05,only.pos = T)

DimPlot(monkey.cellbin.oligo.clean,label = T,split.by = "orig.ident",ncol = 2)

monkey.cellbin.oligo.clean@meta.data$sub = "MOL"
monkey.cellbin.oligo.clean@meta.data$sub[monkey.cellbin.oligo.clean@meta.data$SCT_snn_res.0.2 %in% c(1,2)] = "ARO1"
monkey.cellbin.oligo.clean@meta.data$sub[monkey.cellbin.oligo.clean@meta.data$SCT_snn_res.0.2 %in% c(3)] = "ARO2"

Idents(monkey.cellbin.oligo.clean) = "sub"
levels(monkey.cellbin.oligo.clean) = c("MOL","ARO1","ARO2")

FeaturePlot(monkey.cellbin.oligo,"MOBP")
all <- FindClusters(all, resolution = 0.1, algorithm = 4)

library(dplyr)
monkey.cellbin.oligo.clean@meta.data$cell = rownames(monkey.cellbin.oligo.clean@meta.data)

the.table = data.frame()

for (i in unique(monkey.all.bin50.seurat.clean@meta.data$orig.ident)){
s = monkey.all.bin50.seurat.clean@meta.data[monkey.all.bin50.seurat.clean@meta.data$orig.ident == i,]
s = s[,c("x","y","region_big")]
s2 = monkey.cellbin.oligo.clean@meta.data[monkey.cellbin.oligo.clean@meta.data$orig.ident == i,]

x0 = min(s$x)
y0 = min(s$y)
s2$bin_x = floor((s2$x - x0) / 50) * 50 + x0
s2$bin_y = floor((s2$y - y0) / 50) * 50 + y0

s3 = left_join(
  s2,
  s,
  by = c("bin_x" = "x", "bin_y" = "y")
)
the.table = rbind(the.table,s3)
}
rownames(the.table) = the.table$cell
the.table = the.table[rownames(monkey.cellbin.oligo.clean@meta.data),]

monkey.cellbin.oligo.clean@meta.data$region_big = the.table$region_big


the.table = monkey.cellbin.oligo.clean@meta.data[monkey.cellbin.oligo.clean@meta.data$region_big %in% c("Fimbra","SFO"),]
the.table$group = "old"
the.table$group[the.table$orig.ident %in% c("young1","young2")] = "young"
the.table$new = paste0(the.table$group,"_",the.table$region_big)

the.table$new = factor(the.table$new,levels = c("young_Fimbra","young_SFO","old_Fimbra","old_SFO"),ordered = T)

monkey.cellbin.oligo.clean@meta.data$tmp = paste0(monkey.cellbin.oligo.clean@meta.data$region_big,"_",monkey.cellbin.oligo.clean@meta.data$sub)
Idents(monkey.cellbin.oligo.clean) = "region_big"
DimPlot(monkey.cellbin.oligo.clean,split.by = "region_big",ncol = 3)

saveRDS(monkey.cellbin.oligo.clean,"/Volumes/SFO/The_SFO_project/data/processed/monkey_all4_cellbin_clean_oligosub_clean.rds")
saveRDS(monkey.cellbin.all,"/Volumes/SFO/The_SFO_project/data/processed/monkey_all4_cellbin_clean_clustered.rds")

monkey.cellbin.all = readRDS("/Volumes/SFO/The_SFO_project/data/processed/monkey_all4_cellbin_clean_clustered.rds")
monkey.cellbin.oligo.clean = readRDS("/Volumes/SFO/The_SFO_project/data/processed/monkey_all4_cellbin_clean_oligosub_clean.rds")

monkey.cellbin.all@meta.data$cell = rownames(monkey.cellbin.all@meta.data)
drop.cells = setdiff(monkey.cellbin.all@meta.data$cell[monkey.cellbin.all@meta.data$cluster == "Oligo"],monkey.cellbin.oligo.clean@meta.data$cell)
monkey.cellbin.all.clean = subset(monkey.cellbin.all,cell %in% setdiff(monkey.cellbin.all@meta.data$cell,drop.cells))
monkey.cellbin.all.clean <- RunUMAP(monkey.cellbin.all.clean, dims = 1:30,reduction = "pca")