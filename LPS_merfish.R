### 2026-04-03
### Ruoqing Feng
### LPS treated merfish analysis with control at month 3 age

s5r3.segittal = load_vizgen_parquet(viz_dir = "/Volumes/Papers/202509291148_SFO-slide6_VMSC03901/region_R3/",fov_name = "s5r3",assay = "Vizgen",z = 3L,key = paste0("s5r3","_"))
s5r4.segittal = load_vizgen_parquet(viz_dir = "/Volumes/Papers/202509291148_SFO-slide6_VMSC03901/region_R4/",fov_name = "s5r4",assay = "Vizgen",z = 3L,key = paste0("s5r4","_"))
s5r5.segittal = load_vizgen_parquet(viz_dir = "/Volumes/Papers/202509291148_SFO-slide6_VMSC03901/region_R5/",fov_name = "s5r5",assay = "Vizgen",z = 3L,key = paste0("s5r5","_"))


lps.all = list(s5r3.segittal,s5r4.segittal,s5r5.segittal)
names(lps.all) = c("s5r3","s5r4","s5r5")

lps.all = annotate_samples(lps.all)
lps.all = merge_seurat_list(lps.all, keep_fov = TRUE)

blank_genes <- grep("^Blank-", rownames(lps.all), value = TRUE)
lps.all <- PercentageFeatureSet(
  lps.all, pattern = "^Blank-", assay = "Vizgen",
  col.name = "percent.blank"
)

lps.all$percent.blank[ is.nan(lps.all$percent.blank) ] <- 0
VlnPlot(lps.all,c("nCount_Vizgen","nFeature_Vizgen","volume","percent.blank"),pt.size = 0,ncol = 2)
lps.all = subset(lps.all,nFeature_Vizgen >= 10 & nCount_Vizgen >= 10 & nCount_Vizgen <= 5000 & volume < 6000 & percent.blank <= 10)

saveRDS(lps.all,"/Volumes/SFO/The_SFO_project/data/processed/lps_m3_all.rds")

lps.all.qc <- subset(lps.all, features = setdiff(rownames(lps.all), blank_genes))
saveRDS(lps.all.qc,"/Volumes/SFO/The_SFO_project/data/processed/lps_all_qced.rds")

lps.all.qc = SCTransform(lps.all.qc,clip.range = c(-10,10),assay = "Vizgen")
lps.all.qc = RunPCA(lps.all.qc,npcs = 50)
lps.all.qc = RunUMAP(lps.all.qc,dims = 1:40)


DimPlot(lps.all.qc,label = T)

for (i in unique(lps.all.qc$sample)){
  print(i)
  s = seurat_to_SpatialExperiment(obj = lps.all.qc,sample_id = i)
  saveRDS(s,paste0("/Volumes/SFO/The_SFO_project/data/processed/bamksyse/",i,"_se.rds"))
}


lps.all.qc = FindClusters(lps.all.qc,resolution = 0.8)


DimPlot(lps.all.qc,label = T)
FeaturePlot(lps.all.qc,"P2ry12")
 for (i in 1:3){
   slot(object = lps.all.qc@assays$SCT@SCTModel.list[[i]], name="umi.assay")<-"Vizgen"
 }
 
lps.all.qc@meta.data$cell = rownames(lps.all.qc@meta.data)
lps.all.qc = PrepSCTFindMarkers(lps.all.qc)
lps.marker = FindAllMarkers(lps.all.qc,only.pos = T,min.pct = 0.1)


lps.all.qc@meta.data$tmp = as.numeric(as.character(lps.all.qc@meta.data$SCT_snn_res.0.8))

lps.all.qc@meta.data$celltype = ""
lps.all.qc@meta.data$celltype[lps.all.qc@meta.data$tmp %in% c(0,3,7,18,23,29)] = "Oligo"
lps.all.qc@meta.data$celltype[lps.all.qc@meta.data$tmp %in% c(1,19,22,25)] = "Endo"
lps.all.qc@meta.data$celltype[lps.all.qc@meta.data$tmp %in% c(2,14,15,20)] = "Astro"
lps.all.qc@meta.data$celltype[lps.all.qc@meta.data$tmp %in% c(4,5,8,9,10,11,12,21,27,28)] = "Neuron"
lps.all.qc@meta.data$celltype[lps.all.qc@meta.data$tmp %in% c(6)] = "CPEC"
lps.all.qc@meta.data$celltype[lps.all.qc@meta.data$tmp %in% c(13,26)] = "Ependyma"
lps.all.qc@meta.data$celltype[lps.all.qc@meta.data$tmp %in% c(16)] = "Micro"
lps.all.qc@meta.data$celltype[lps.all.qc@meta.data$tmp %in% c(17)] = "OPC"
lps.all.qc@meta.data$celltype[lps.all.qc@meta.data$tmp %in% c(24)] = "SFO"


lps.endo.qc = subset(lps.all.qc,celltype == "Endo")
lps.endo.qc = SCTransform(lps.endo.qc,assay = "Vizgen")
lps.endo.qc = RunPCA(lps.endo.qc)
lps.endo.qc = RunUMAP(lps.endo.qc,dims = 1:30)
lps.endo.qc = FindClusters(lps.endo.qc,resolution = 0.5)

saveRDS(lps.endo.qc,"/Volumes/SFO/The_SFO_project/data/processed/lps_m3_endo.rds")

lps.all.qc@meta.data$celltype[lps.all.qc@meta.data$cell %in% lps.endo.qc@meta.data$cell[lps.endo.qc@meta.data$SCT_snn_res.0.5 == 6]] = "Mac"

sfo.major.color.cross = c("#D3B256","#4182CB","#D78203","#BF99E2","#5F6164","#373F89","#4CBDA8","#90D4A4","#EFA9AE","#5066a1")
names(sfo.major.color.cross) = c("Neuron","Oligo","Astro","Ependymal","CPEC","Micro","Endo","OPC","SFO","Mac")

p1 = DimPlot(lps.all.qc,cols = sfo.major.color.cross,pt.size = 0.3)+p.cleanumap+NoLegend()
ggsave("/Users/ruoqing/Projects/sfo/fig_sn_merfish/lps_merfish_cross_all.png",p1,width = 9,height = 6,dpi = 300)
the.legend = plot_umap_legend(sfo.major.color.cross)
ggsave("/Users/ruoqing/Projects/sfo/fig_sn_merfish/lps_cross_all_umap_legend.pdf",the.legend,width = 6,height = 4,dpi = 300)

Idents(lps.all.qc) = "celltype"
DimPlot(lps.all.qc,label = T)

lps.mac_micro.qc = subset(lps.all.qc,celltype %in% c("Mac","Micro"))
lps.mac_micro.qc = SCTransform(lps.mac_micro.qc,assay = "Vizgen")
lps.mac_micro.qc = RunPCA(lps.mac_micro.qc)
lps.mac_micro.qc = RunUMAP(lps.mac_micro.qc,dims = 1:30)

lps.mac_micro.qc = subset(lps.mac_micro.qc,celltype %in% c("Mac","Micro"))
lps.mac_micro.qc@meta.data$mac_sub = lps.mac_micro.qc@meta.data$celltype
saveRDS(lps.mac_micro.qc,"/Volumes/SFO/The_SFO_project/data/processed/lps_m3_mac_micro.rds")


mac_sub.ctrl = subset(sfo.merfish.qc.keep.update,celltype_clean == "Micro")
FeaturePlot(mac_sub.ctrl,"Cd68")
FeaturePlot(lps.mac_micro.qc,"Trem2")

DimPlot(lps.mac_micro.qc,label = T)
lps.mac_micro.qc = FindClusters(lps.mac_micro.qc,resolution = 0.3)

lps_ctrl.mac_micro.qc = merge(mac_sub.ctrl,lps.mac_micro.qc)
lps_ctrl.mac_micro.qc = subset(lps_ctrl.mac_micro.qc,sample %in% c("sfo2r3","sfo2r6","sfo4r3","sfo4r6","s5r3","s5r4","s5r5"))

lps_ctrl.mac_micro.qc@meta.data$group = "ctrl"
lps_ctrl.mac_micro.qc@meta.data$group[lps_ctrl.mac_micro.qc@meta.data$sample %in% c("s5r3","s5r4","s5r5")] = "lps"
lps_ctrl.mac_micro.qc = SCTransform(lps_ctrl.mac_micro.qc,assay = "Vizgen")
lps_ctrl.mac_micro.qc = RunPCA(lps_ctrl.mac_micro.qc)
lps_ctrl.mac_micro.qc = RunUMAP(lps_ctrl.mac_micro.qc,dims = 1:30)

lps_ctrl.mac_micro.qc <- FindNeighbors(lps_ctrl.mac_micro.qc, reduction = "pca", dims = 1:30)
lps_ctrl.mac_micro.qc <- FindClusters(lps_ctrl.mac_micro.qc, resolution = 0.3)
DimPlot(lps_ctrl.mac_micro.qc,label = T,split.by = "group")
tmp = FindAllMarkers(lps_ctrl.mac_micro.qc,min.pct = 0.1,only.pos = T)

FeaturePlot(lps_ctrl.mac_micro.qc,"Cd68",split.by = "sample")
ImageDimPlot(lps_ctrl.mac_micro.qc,fov = "sfo2r3",split.by = "SCT_snn_res.0.3")

Idents(lps_ctrl.mac_micro.qc) = "group"
for (i in 1:2){
  slot(object = lps_ctrl.mac_micro.qc@assays$SCT@SCTModel.list[[i]], name="umi.assay")<-"Vizgen"
}
lps_ctrl.mac_micro.qc = PrepSCTFindMarkers(lps_ctrl.mac_micro.qc)
micro.diff = FindMarkers(lps_ctrl.mac_micro.qc,ident.1 = "lps",ident.2 = "ctrl")

micro.diff$gene = rownames(micro.diff)
micro.diff$sig = "non"
micro.diff$sig[micro.diff$avg_log2FC >= 0.58 & micro.diff$p_val_adj <= 0.05] = "up"
micro.diff$sig[micro.diff$avg_log2FC <= (-0.58) & micro.diff$p_val_adj <= 0.05] = "down"
micro.diff$diff = micro.diff$pct.1 - micro.diff$pct.2
micro.diff$label = ""
label.gene = c("Socs3","Junb","Cxcl10","Icam1","Tnf","Cd83","Jak1","Irf7","Sall1","P2ry12","Cx3cr1")
for(i in label.gene){
  micro.diff$label[micro.diff$gene == i] = i
}


saveRDS(lps_ctrl.mac_micro.qc,"/Volumes/SFO/The_SFO_project/data/processed/lps_ctrl_mac_micro_qc.rds")
saveRDS(micro.diff,"/Volumes/SFO/The_SFO_project/data/interim/micro_diff_lps_ctrl.rds")


stress.micro = c("Socs3","Junb","Cxcl10","Icam1","Tnf","Cd83","Jak1","Irf7")

lps_ctrl.mac_micro.qc = AddModuleScore(
  object   = lps_ctrl.mac_micro.qc,
  features = list(stress.micro = stress.micro),
  name     = "Stress_Score",ctrl     = 10)

m3.merfish.ctrl_lps.seurat@meta.data$Stress_Score = -2
for (i in lps_ctrl.mac_micro.qc@meta.data$cell){
  m3.merfish.ctrl_lps.seurat@meta.data$Stress_Score[m3.merfish.ctrl_lps.seurat@meta.data$cell == i] = lps_ctrl.mac_micro.qc@meta.data$Stress_Score1[lps_ctrl.mac_micro.qc@meta.data$cell == i]
}
for (i in c("sfo2r3","sfo2r6","sfo4r3","sfo4r6","s5r3","s5r4","s5r5")){
  p1 = ImageFeaturePlot(m3.merfish.ctrl_lps.seurat,fov = c(i),features = "Stress_Score",boundaries = "centroids", axes = TRUE,size = 0.5,border.size = 0,dark.background = F)+scale_fill_gradientn(colours = c( "grey88","grey88","pink","red", "black"),values  = rescale(c(-2,  -0.8282041,0,0.5, 1.7)),limits = c(-2, 1.7))+ggtitle("")+
  theme_classic()+NoLegend()
  ggsave(paste0("/Users/ruoqing/Projects/sfo/fig_sn_merfish/str_micro/",i,".png"),p1,width = 6,height = 6,dpi = 300)
  p1 = ImageFeaturePlot(m3.merfish.ctrl_lps.seurat,fov = c(i),features = "Stress_Score",boundaries = "centroids", axes = TRUE,size = 1.5,border.size = 0,dark.background = F)+scale_fill_gradientn(colours = c( "grey90","grey90","blue", "red"),values  = rescale(c(-2,  -0.8282041,0, 1.7)),limits = c(-2, 1.7))+ggtitle("")+
    theme_classic()
  ggsave(paste0("/Users/ruoqing/Projects/sfo/fig_sn_merfish/str_micro/",i,"_legend.pdf"),p1,width = 6,height = 6)
}

tmp = data.frame()
for (i in unique(m3.merfish.ctrl_lps.seurat@meta.data$sample)){
  s = as.data.frame(m3.merfish.ctrl_lps.seurat@images[[i]]$centroids@coords)
  rownames(s) = m3.merfish.ctrl_lps.seurat@images[[i]]$centroids@cells
  tmp = rbind(tmp,s)
}
tmp = tmp[rownames(m3.merfish.ctrl_lps.seurat@meta.data),]
m3.merfish.ctrl_lps.seurat@meta.data$x = tmp$x
m3.merfish.ctrl_lps.seurat@meta.data$y = tmp$y

tmp = data.frame()
for (i in unique(lps_ctrl.mac_micro.qc@meta.data$sample)){
  s = as.data.frame(lps_ctrl.mac_micro.qc@images[[i]]$centroids@coords)
  rownames(s) = lps_ctrl.mac_micro.qc@images[[i]]$centroids@cells
  tmp = rbind(tmp,s)
}
tmp = tmp[rownames(lps_ctrl.mac_micro.qc@meta.data),]
lps_ctrl.mac_micro.qc@meta.data$x = tmp$x
lps_ctrl.mac_micro.qc@meta.data$y = tmp$y

}


m3.merfish.ctrl_lps.seurat@meta.data$cell = rownames(m3.merfish.ctrl_lps.seurat@meta.data)
saveRDS(m3.merfish.ctrl_lps.seurat,"/Volumes/SFO/The_SFO_project/data/processed/m3_merfish_ctrl_lps_seurat.rds")
m3.merfish.ctrl_lps.seurat = readRDS("/Volumes/SFO/The_SFO_project/data/processed/m3_merfish_ctrl_lps_seurat.rds")

DimPlot(m3.merfish.ctrl_lps.seurat)


lps.oligo.qc = subset(lps.all.qc,celltype %in% c("Oligo"))
lps.oligo.qc = SCTransform(lps.oligo.qc,assay = "Vizgen")
lps.oligo.qc = RunPCA(lps.oligo.qc)
lps.oligo.qc = RunUMAP(lps.oligo.qc,dims = 1:30)
lps.oligo.qc = FindClusters(lps.oligo.qc,resolution = 0.3)

saveRDS(lps.oligo.qc,"/Volumes/SFO/The_SFO_project/data/processed/lps_m3_oligo.rds")


oligo.ctrl = subset(sfo.merfish.qc.keep.update,celltype_clean == "Oligo")
oligo.ctrl = subset(oligo.ctrl,sample %in% c("sfo2r3","sfo2r6","sfo4r3","sfo4r6"))

lps_ctrl.oligo.qc = merge(oligo.ctrl,lps.oligo.qc)

lps_ctrl.oligo.qc@meta.data$group = "ctrl"
lps_ctrl.oligo.qc@meta.data$group[lps_ctrl.oligo.qc@meta.data$sample %in% c("s5r3","s5r4","s5r5")] = "lps"
lps_ctrl.oligo.qc = SCTransform(lps_ctrl.oligo.qc,assay = "Vizgen")
lps_ctrl.oligo.qc = RunPCA(lps_ctrl.oligo.qc)
lps_ctrl.oligo.qc = RunUMAP(lps_ctrl.oligo.qc,dims = 1:30)

lps_ctrl.oligo.qc <- FindNeighbors(lps_ctrl.oligo.qc, reduction = "pca", dims = 1:30)
lps_ctrl.oligo.qc <- FindClusters(lps_ctrl.oligo.qc, resolution = 1)
DimPlot(lps_ctrl.oligo.qc,label = T,split.by = "SCT_snn_res.1",ncol = 4)
DimPlot(lps_ctrl.oligo.qc,label = T)
for (i in 1:2){
  slot(object = lps_ctrl.oligo.qc@assays$SCT@SCTModel.list[[i]], name="umi.assay")<-"Vizgen"
}
lps_ctrl.oligo.qc = PrepSCTFindMarkers(lps_ctrl.oligo.qc)

FeaturePlot(lps_ctrl.oligo.qc,"Myrf",split.by = "group")

lps_ctrl.oligo.qc@meta.data$tmp = as.data.frame(as.character(lps_ctrl.oligo.qc@meta.data$SCT_snn_res.0.1))
Idents(lps_ctrl.oligo.qc) = "SCT_snn_res.0.1"
tmp.marker = FindAllMarkers(lps_ctrl.oligo.qc,min.pct = 0.1,only.pos = T)

lps_ctrl.oligo.qc@meta.data$oligo_sub = ""
lps_ctrl.oligo.qc@meta.data$oligo_sub[lps_ctrl.oligo.qc@meta.data$SCT_snn_res.0.1 == 1] = "COP"
lps_ctrl.oligo.qc@meta.data$oligo_sub[lps_ctrl.oligo.qc@meta.data$SCT_snn_res.0.1 == 0] = "MOL"
lps_ctrl.oligo.qc@meta.data$oligo_sub[lps_ctrl.oligo.qc@meta.data$SCT_snn_res.0.1 == 2] = "MOL"
lps_ctrl.oligo.qc@meta.data$oligo_sub[lps_ctrl.oligo.qc@meta.data$SCT_snn_res.0.1 == 3] = "STROL"
Idents(lps_ctrl.oligo.qc) = "oligo_sub"

oligo.color = c("darkorange", "purple", "cyan4","#BD6263")
names(oligo.color) = c("COP","MOL","OPC","STROL")

p1 = DimPlot(lps_ctrl.oligo.qc,cols = oligo.color)+p.cleanumap+NoLegend()
ggsave("/Users/ruoqing/Projects/sfo/fig_sn_merfish/ctrl_lps_oligo_umap.png",p1,width = 6,height = 4)
p1 = plot_umap_legend(oligo.color)
ggsave("/Users/ruoqing/Projects/sfo/fig_sn_merfish/ctrl_lps_oligo_umap_legend.pdf",p1,width = 6,height = 4)

saveRDS(lps_ctrl.oligo.qc,"/Volumes/SFO/The_SFO_project/data/processed/lps_ctrl_oligo_qc.rds")


use.marker = c("Mobp","Bcas1","Atp1a2",
               "Mog","Myrf","Mag","Opalin",
               "Serpina3n","C4b","Fos","Irf7")
tmp = AverageExpression(lps_ctrl.oligo.qc,assays = "SCT",features = use.marker,return.seurat = T)
DoHeatmap(tmp,features = use.marker)

tmp = tmp@assays$SCT$scale.data
bk = seq(-2,2,0.1)
p2 = pheatmap(tmp,breaks = bk,
              color = colorRampPalette(c("#5C88DA99","white", "#CC0C0099"))(length(bk)),
              cluster_cols = F,cluster_rows = F,border_color = NA,show_colnames = T)
ggsave("/Users/ruoqing/Projects/sfo/fig_sn_merfish/ctrl_lps_oligo_marker.pdf",p2,width = 2,height = 3)

oligo.color = c("darkorange", "purple", "cyan4","#BD6263","grey88")
names(oligo.color) = c("COP","MOL","OPC","STROL","others")

sfo.merfish.qc.keep.update@meta.data$oligo_sub = "others"
for (i in unique(lps_ctrl.oligo.qc@meta.data$oligo_sub)){
  sfo.merfish.qc.keep.update@meta.data$oligo_sub[sfo.merfish.qc.keep.update@meta.data$cell %in% lps_ctrl.oligo.qc@meta.data$cell[lps_ctrl.oligo.qc@meta.data$oligo_sub == i]] = i
}
Idents(sfo.merfish.qc.keep.update) = "oligo_sub"

for (i in c("sfo2r3","sfo2r6","sfo4r3","sfo4r6")){
p1 = ImageDimPlot(sfo.merfish.qc.keep.update, fov = i, boundaries = "centroids", axes = TRUE,size = 0.8,border.size = 0,
                  cols = oligo.color,dark.background = F)+
  theme_classic()+NoLegend()+
  theme(axis.title = element_blank())
ggsave(paste0("/Users/ruoqing/Projects/sfo/fig_sn_merfish/STROL_projection/",i,".png"),p1,width = 6,height = 6,dpi = 300)
}

sfo.merfish.qc.keep.update@meta.data$oligo_sub = "others"
sfo.merfish.qc.keep.update@meta.data$oligo_sub[sfo.merfish.qc.keep.update@meta.data$cell %in% lps_ctrl.oligo.qc@meta.data$cell[lps_ctrl.oligo.qc@meta.data$oligo_sub == "STROL"]] = "STROL"

Idents(sfo.merfish.qc.keep.update) = "oligo_sub"
for (i in c("sfo2r3","sfo2r6","sfo4r3","sfo4r6")){
  p1 = ImageDimPlot(sfo.merfish.qc.keep.update, fov = i, boundaries = "centroids", axes = TRUE,size = 0.8,border.size = 0,
                    cols = oligo.color,dark.background = F)+
    theme_classic()+NoLegend()+
    theme(axis.title = element_blank())
  ggsave(paste0("/Users/ruoqing/Projects/sfo/fig_sn_merfish/STROL_projection/",i,"_onlystr.png"),p1,width = 6,height = 6,dpi = 300)
}






for (i in c("s5r3","s5r4","s5r5")){
  p1 = ImageDimPlot(lps.all.qc, fov = i, boundaries = "centroids", axes = TRUE,size = 0.8,border.size = 0,
                    cols = oligo.color,dark.background = F)+
    theme_classic()+NoLegend()+
    theme(axis.title = element_blank())
  ggsave(paste0("/Users/ruoqing/Projects/sfo/fig_sn_merfish/STROL_projection/",i,".png"),p1,width = 6,height = 6,dpi = 300)
}


lps.all.qc@meta.data$oligo_sub = "others"
lps.all.qc@meta.data$oligo_sub[lps.all.qc@meta.data$cell %in% lps_ctrl.oligo.qc@meta.data$cell[lps_ctrl.oligo.qc@meta.data$oligo_sub == "STROL"]] = "STROL"
Idents(lps.all.qc) = "oligo_sub"

for (i in c("s5r3","s5r4","s5r5")){
  p1 = ImageDimPlot(lps.all.qc, fov = i, boundaries = "centroids", axes = TRUE,size = 0.8,border.size = 0,
                    cols = oligo.color,dark.background = F)+
    theme_classic()+NoLegend()+
    theme(axis.title = element_blank())
  ggsave(paste0("/Users/ruoqing/Projects/sfo/fig_sn_merfish/STROL_projection/",i,"_onlystr.png"),p1,width = 6,height = 6,dpi = 300)
}


DimPlot(lps.oligo.qc,label = T)
lps.oligo.marker = FindAllMarkers(lps.oligo.qc,min.pct = 0.1,only.pos = T)
FeaturePlot(lps.oligo.qc,"Bcas1")


ImageDimPlot(lps.oligo.qc,fov = "s5r3",split.by = "SCT_snn_res.0.3")

############################################################################################
i = "s5r3"
tmp = readRDS(paste0("/Volumes/SFO/The_SFO_project/data/processed/bamksyse/",i,"_se.rds"))
tmp = banksy_pipeline(tmp,resolution = 0.2)
colData(tmp)$region = as.character(colData(tmp)$clust_M1_lam0.8_k50_res0.2)
colData(tmp)$cells = rownames(colData(tmp))

colData(tmp)$region[colData(tmp)$clust_M1_lam0.8_k50_res0.2 == 8] = "SFO"
colData(tmp)$region[colData(tmp)$clust_M1_lam0.8_k50_res0.2 == 1] = "CTX"
colData(tmp)$region[colData(tmp)$clust_M1_lam0.8_k50_res0.2 %in% c(5,6)] = "ventricular"
colData(tmp)$region[colData(tmp)$clust_M1_lam0.8_k50_res0.2 %in% c(2)] = "fibro track"
colData(tmp)$region[colData(tmp)$clust_M1_lam0.8_k50_res0.2 %in% c(3)] = "TH"
colData(tmp)$region[colData(tmp)$clust_M1_lam0.8_k50_res0.2 %in% c(7)] = "HY"
colData(tmp)$region[colData(tmp)$clust_M1_lam0.8_k50_res0.2 %in% c(4)] = "CNU"

saveRDS(tmp,paste0("/Volumes/SFO/MERFISHRes/cleanV2/banksyobo/",i,"_se_banksyed_anno.rds"))
############################################################################################
i = "s5r4"
tmp = readRDS(paste0("/Volumes/SFO/The_SFO_project/data/processed/bamksyse/",i,"_se.rds"))
tmp = banksy_pipeline(tmp,resolution = 0.2)
colData(tmp)$region = as.character(colData(tmp)$clust_M1_lam0.8_k50_res0.2)
colData(tmp)$cells = rownames(colData(tmp))

colData(tmp)$region[colData(tmp)$clust_M1_lam0.8_k50_res0.2 == 7] = "SFO"
colData(tmp)$region[colData(tmp)$clust_M1_lam0.8_k50_res0.2 == 1] = "CTX"
colData(tmp)$region[colData(tmp)$clust_M1_lam0.8_k50_res0.2 %in% c(3)] = "ventricular"
colData(tmp)$region[colData(tmp)$clust_M1_lam0.8_k50_res0.2 %in% c(5)] = "CNU"
colData(tmp)$region[colData(tmp)$clust_M1_lam0.8_k50_res0.2 %in% c(2)] = "fibro track"
colData(tmp)$region[colData(tmp)$clust_M1_lam0.8_k50_res0.2 %in% c(4,6)] = "TH"

saveRDS(tmp,paste0("/Volumes/SFO/MERFISHRes/cleanV2/banksyobo/",i,"_se_banksyed_anno.rds"))
############################################################################################
i = "s5r5"
tmp = readRDS(paste0("/Volumes/SFO/The_SFO_project/data/processed/bamksyse/",i,"_se.rds"))
tmp = banksy_pipeline(tmp,resolution = 0.5)
colData(tmp)$region = as.character(colData(tmp)$clust_M1_lam0.8_k50_res0.5)
colData(tmp)$cells = rownames(colData(tmp))

colData(tmp)$region[colData(tmp)$clust_M1_lam0.8_k50_res0.5 %in% c(1,11,10)] = "fibro track"
colData(tmp)$region[colData(tmp)$clust_M1_lam0.8_k50_res0.5 %in% c(3,4,15,13)] = "CTX"
colData(tmp)$region[colData(tmp)$clust_M1_lam0.8_k50_res0.5 %in% c(9)] = "CNU"
colData(tmp)$region[colData(tmp)$clust_M1_lam0.8_k50_res0.5 %in% c(8,6,5,12)] = "ventricular"
colData(tmp)$region[colData(tmp)$clust_M1_lam0.8_k50_res0.5 == 14] = "SFO"
colData(tmp)$region[colData(tmp)$clust_M1_lam0.8_k50_res0.5 %in% c(2,7)] = "TH"

saveRDS(tmp,paste0("/Volumes/SFO/MERFISHRes/cleanV2/banksyobo/",i,"_se_banksyed_anno.rds"))



keep.section = c("s5r3","s5r4","s5r5")
coldata.merge = data.frame()
for (i in keep.section){
  tmp = readRDS(paste0("/Volumes/SFO/MERFISHRes/cleanV2/banksyobo/",i,"_se_banksyed_anno.rds"))
  tmp = as.data.frame(colData(tmp))
  tmp = tmp[,c("sample","region","cells")]
  coldata.merge = rbind(coldata.merge,tmp)
}

for (i in unique(coldata.merge$region)){
  print(i)
  m3.merfish.ctrl_lps.seurat@meta.data$region_new[m3.merfish.ctrl_lps.seurat@meta.data$cell %in% coldata.merge$cells[coldata.merge$region == i]] = i
}

for (i in unique(coldata.merge$region)){
  print(i)
  m3.merfish.ctrl_lps.seurat@meta.data$region_new_cc[m3.merfish.ctrl_lps.seurat@meta.data$cell %in% coldata.merge$cells[coldata.merge$region == i]] = i
}


######################################################################################################

sfo.merfish.qc.keep.update = readRDS("/Volumes/SFO/The_SFO_project/data/processed/sfo_merfish_cross_all_anno_banbksy_update_cc.rds")

for (i in unique(sfo.merfish.qc.keep.update@meta.data$region_new_cc)){
  print(i)
  m3.merfish.ctrl_lps.seurat@meta.data$region_new_cc[m3.merfish.ctrl_lps.seurat@meta.data$cell %in% sfo.merfish.qc.keep.update@meta.data$cell[sfo.merfish.qc.keep.update@meta.data$region_new_cc == i]] = i
}


lps_ctrl.oligo.qc = readRDS("/Volumes/SFO/The_SFO_project/data/processed/lps_ctrl_oligo_qc.rds")

m3.merfish.ctrl_lps.seurat@meta.data$oligo_sub = "other"
for (i in unique(lps_ctrl.oligo.qc@meta.data$oligo_sub)) {
  m3.merfish.ctrl_lps.seurat@meta.data$oligo_sub[m3.merfish.ctrl_lps.seurat@meta.data$cell %in% lps_ctrl.oligo.qc@meta.data$cell[lps_ctrl.oligo.qc@meta.data$oligo_sub == i]] = i
}


lps_ctrl.mac_micro.qc = readRDS("/Volumes/SFO/The_SFO_project/data/processed/lps_ctrl_mac_micro_qc.rds")

tmp = subset(m3.merfish.ctrl_lps.seurat,region_new_cc == "SFO")

Idents(tmp) = "oligo_sub"

for (i in unique(tmp$sample)){
  print(i)
  p1 = ImageDimPlot(tmp,fov = i,cols = merfish.oligo.sub.cols,dark.background = F,axes = T,border.size = 0.1)+
  theme_classic()+NoLegend()+xlab("")+ylab("")
  ggsave(paste0("/Users/ruoqing/Projects/sfo/fig_sn_merfish/lps/",i,".pdf"),p1,width = 5,height = 5)
}

saveRDS(m3.merfish.ctrl_lps.seurat,"/Volumes/SFO/The_SFO_project/data/processed/m3_merfish_ctrl_lps_seurat.rds")
saveRDS(lps_ctrl.mac_micro.qc,"/Volumes/SFO/The_SFO_project/data/processed/lps_ctrl_mac_micro_qc.rds")
saveRDS(lps_ctrl.oligo.qc,"/Volumes/SFO/The_SFO_project/data/processed/lps_ctrl_oligo_qc.rds")

m3.merfish.ctrl_lps.seurat = readRDS("/Volumes/SFO/The_SFO_project/data/processed/m3_merfish_ctrl_lps_seurat.rds")
lps_ctrl.oligo.qc = readRDS("/Volumes/SFO/The_SFO_project/data/processed/lps_ctrl_oligo_qc.rds")