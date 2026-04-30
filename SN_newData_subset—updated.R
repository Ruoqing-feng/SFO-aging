### Ruoqing-Feng
library(Seurat)

l1 = Read10X_h5("/Volumes/SFO/SN/mtx/l1/filtered_feature_bc_matrix.h5")
l2 = Read10X_h5("/Volumes/SFO/SN/mtx/l2/filtered_feature_bc_matrix.h5")
l3 = Read10X_h5("/Volumes/SFO/SN/mtx/l3/filtered_feature_bc_matrix.h5")
l4 = Read10X_h5("/Volumes/SFO/SN/mtx/l4/filtered_feature_bc_matrix.h5")
l5 = Read10X_h5("/Volumes/SFO/SN/mtx/l5/filtered_feature_bc_matrix.h5")
l6 = Read10X_h5("/Volumes/SFO/SN/mtx/l6/filtered_feature_bc_matrix.h5")

colnames(l1) = paste0(colnames(l1),"-l1")
colnames(l2) = paste0(colnames(l2),"-l2")
colnames(l3) = paste0(colnames(l3),"-l3")
colnames(l4) = paste0(colnames(l4),"-l4")
colnames(l5) = paste0(colnames(l5),"-l5")
colnames(l6) = paste0(colnames(l6),"-l6")

tmp = cbind(l1,l2)
tmp = cbind(tmp,l3)
tmp = cbind(tmp,l4)
tmp = cbind(tmp,l5)
tmp = cbind(tmp,l6)

tmp = CreateSeuratObject(tmp)


tail_part <- sub(".*-(.*)$", "\\1", colnames(tmp))
tmp@meta.data$sample = tail_part
Idents(tmp) = "sample"

VlnPlot(tmp,c("nCount_RNA","nFeature_RNA"),pt.size = 0)


tmp.clean = subset(tmp,nCount_RNA <= 15000 & nFeature_RNA <= 5000 & nFeature_RNA >= 400)

VlnPlot(tmp.clean,c("nCount_RNA","nFeature_RNA"),pt.size = 0)


tmp.clean = NormalizeData(tmp.clean)
tmp.clean = FindVariableFeatures(tmp.clean,nfeatures = 3000)
tmp.clean = ScaleData(tmp.clean)

tmp.clean = RunPCA(tmp.clean)
ElbowPlot(tmp.clean,ndims = 50)

tmp.clean = RunUMAP(tmp.clean,dims = 1:30,reduction = 'pca')
Idents(tmp.clean) = "sample"

DimPlot(tmp.clean,label = T,split.by = "sample",ncol = 3)

FeaturePlot(tmp.clean,"C3")
saveRDS(tmp.clean,"/Volumes/SFO/snAnalysis/new_seqdata_clean.rds")

tmp.clean  = readRDS("/Volumes/SFO/snAnalysis/new_seqdata_clean.rds")


write.csv(tmp.clean@meta.data,"/Volumes/SFO/snAnalysis/new_sn_clean.csv",quote = F)

FeaturePlot(tmp.clean,"Rbfox1")
DimPlot(tmp.clean,label = T)

tmp.clean = FindNeighbors(tmp.clean,assay = "RNA",reduction = "pca",dims = 1:30)
tmp.clean <- FindClusters(tmp.clean, resolution = 0.3)

tmp.markers = FindAllMarkers(tmp.clean,min.pct = 0.1,only.pos = T)

tmp.clean@meta.data$tmp = as.numeric(as.character(tmp.clean@meta.data$RNA_snn_res.0.3))

tmp.clean@meta.data$RNA_cluster = "Neuron"
tmp.clean@meta.data$RNA_cluster[tmp.clean@meta.data$tmp %in% c(0)] = "Oligo"
tmp.clean@meta.data$RNA_cluster[tmp.clean@meta.data$tmp %in% c(1,21)] = "ChP"
tmp.clean@meta.data$RNA_cluster[tmp.clean@meta.data$tmp %in% c(6,19)] = "Astro"
tmp.clean@meta.data$RNA_cluster[tmp.clean@meta.data$tmp %in% c(11)] = "Tancyte-Astro"
tmp.clean@meta.data$RNA_cluster[tmp.clean@meta.data$tmp %in% c(8)] = "PVF"
tmp.clean@meta.data$RNA_cluster[tmp.clean@meta.data$tmp %in% c(12,24)] = "Micro/Mac"
tmp.clean@meta.data$RNA_cluster[tmp.clean@meta.data$tmp %in% c(13)] = "OPC"
tmp.clean@meta.data$RNA_cluster[tmp.clean@meta.data$tmp %in% c(14)] = "Ependymal"
tmp.clean@meta.data$RNA_cluster[tmp.clean@meta.data$tmp %in% c(20)] = "Endo"
tmp.clean@meta.data$RNA_cluster[tmp.clean@meta.data$tmp %in% c(27)] = "MSC"

Idents(tmp.clean) = "RNA_cluster"
p1 = DimPlot(tmp.clean,label = T)+NoLegend()
ggsave("/Users/ruoqing/Projects/sfo/sn_new/2d.png",p1,width = 8,height = 6,dpi = 300)

saveRDS(tmp.clean,"/Volumes/SFO/snAnalysis/new/all_clean_seurat.rds")

sfo.new.all.seurat = readRDS("/Volumes/SFO/snAnalysis/new/all_clean_seurat.rds")
DimPlot(sfo.new.all.seurat,label = T)

sfo.new.neuron.seurat = subset(sfo.new.all.seurat,RNA_cluster == "Neuron")
sfo.new.neuron.seurat = NormalizeData(sfo.new.neuron.seurat)
sfo.new.neuron.seurat = FindVariableFeatures(sfo.new.neuron.seurat,nfeatures = 2000)
sfo.new.neuron.seurat = ScaleData(sfo.new.neuron.seurat)

sfo.new.neuron.seurat = RunPCA(sfo.new.neuron.seurat)
sfo.new.neuron.seurat = RunUMAP(sfo.new.neuron.seurat,dims = 1:30,reduction = 'pca')

sfo.new.neuron.seurat = FindNeighbors(sfo.new.neuron.seurat,assay = "RNA",reduction = "pca",dims = 1:30)
sfo.new.neuron.seurat <- FindClusters(sfo.new.neuron.seurat, resolution = 0.5)

DimPlot(sfo.new.neuron.seurat,label = T)+NoLegend()
saveRDS(sfo.new.neuron.seurat,"/Volumes/SFO/snAnalysis/new/neuron_clean_seurat.rds")

sfo.new.neuron.seurat = readRDS("/Volumes/SFO/snAnalysis/new/neuron_clean_seurat.rds")
###########################################################################################

sfo.new.oligo.seurat = subset(sfo.new.all.seurat,RNA_cluster == "Oligo")
sfo.new.oligo.seurat = NormalizeData(sfo.new.oligo.seurat)
sfo.new.oligo.seurat = FindVariableFeatures(sfo.new.oligo.seurat,nfeatures = 2000)
sfo.new.oligo.seurat = ScaleData(sfo.new.oligo.seurat)

sfo.new.oligo.seurat = RunPCA(sfo.new.oligo.seurat)
sfo.new.oligo.seurat = RunUMAP(sfo.new.oligo.seurat,dims = 1:15,reduction = 'pca')

sfo.new.oligo.seurat = FindNeighbors(sfo.new.oligo.seurat,assay = "RNA",reduction = "pca",dims = 1:15)
sfo.new.oligo.seurat <- FindClusters(sfo.new.oligo.seurat, resolution = 1)

DimPlot(sfo.new.oligo.seurat,label = T,split.by = "sample",ncol = 3)+NoLegend()
FeaturePlot(sfo.new.oligo.seurat,"Serpina3n")
tmp.marker = FindAllMarkers(sfo.new.oligo.seurat,min.pct = 0.1,only.pos = T)

saveRDS(sfo.new.oligo.seurat,"/Volumes/SFO/snAnalysis/new/oligo_clean_seurat.rds")

###########################################################################################

sfo.new.micro.seurat = subset(sfo.new.all.seurat,RNA_cluster == "Micro/Mac")
sfo.new.micro.seurat = NormalizeData(sfo.new.micro.seurat)
sfo.new.micro.seurat = FindVariableFeatures(sfo.new.micro.seurat,nfeatures = 2000)
sfo.new.micro.seurat = ScaleData(sfo.new.micro.seurat)

sfo.new.micro.seurat = RunPCA(sfo.new.micro.seurat)
sfo.new.micro.seurat = RunUMAP(sfo.new.micro.seurat,dims = 1:15,reduction = 'pca')

sfo.new.micro.seurat = FindNeighbors(sfo.new.micro.seurat,assay = "RNA",reduction = "pca",dims = 1:15)
sfo.new.micro.seurat <- FindClusters(sfo.new.micro.seurat, resolution = 0.3)


sfo.new.micro.seurat = readRDS("/Volumes/SFO/snAnalysis/new/imm_clean_seurat.rds")
sfo.new.micro.seurat.clean = subset(sfo.new.micro.seurat,RNA_snn_res.0.3 %in% c(1,2,3,4))
sfo.new.micro.seurat.clean = NormalizeData(sfo.new.micro.seurat.clean)
sfo.new.micro.seurat.clean = FindVariableFeatures(sfo.new.micro.seurat.clean,nfeatures = 1000)
sfo.new.micro.seurat.clean = ScaleData(sfo.new.micro.seurat.clean)
sfo.new.micro.seurat.clean = RunPCA(sfo.new.micro.seurat.clean)
sfo.new.micro.seurat.clean = FindNeighbors(sfo.new.micro.seurat.clean,assay = "RNA",reduction = "pca",dims = 1:15)
sfo.new.micro.seurat.clean  = RunUMAP(sfo.new.micro.seurat.clean,dims = 1:15,reduction = 'pca')
sfo.new.micro.seurat.clean <- FindClusters(sfo.new.micro.seurat.clean, resolution = 0.5)

sfo.new.micro.seurat.clean@meta.data$tmp = sfo.new.micro.seurat.clean@meta.data$RNA_snn_res.0.5
sfo.new.micro.seurat.clean@meta.data$tmp[sfo.new.micro.seurat.clean@meta.data$tmp == 3] = 1
Idents(sfo.new.micro.seurat.clean) = "tmp"
DimPlot(sfo.new.micro.seurat.clean,label = T)+NoLegend()
DimPlot(sfo.new.micro.seurat.clean,label = T,split.by = "sample",ncol = 3)+NoLegend()

FeaturePlot(sfo.new.micro.seurat.clean,"F13a1")



sfo.new.micro.seurat.clean@meta.data$sub = ""
sfo.new.micro.seurat.clean@meta.data$sub[sfo.new.micro.seurat.clean@meta.data$RNA_snn_res.0.5 %in% c(0,3)] = "Micro.Hemeo"
sfo.new.micro.seurat.clean@meta.data$sub[sfo.new.micro.seurat.clean@meta.data$RNA_snn_res.0.5 %in% c(1)] = "pre-DAM"
sfo.new.micro.seurat.clean@meta.data$sub[sfo.new.micro.seurat.clean@meta.data$RNA_snn_res.0.5 %in% c(2)] = "DAM"
sfo.new.micro.seurat.clean@meta.data$sub[sfo.new.micro.seurat.clean@meta.data$RNA_snn_res.0.5 %in% c(4)] = "Mac.MHChigh"
sfo.new.micro.seurat.clean@meta.data$sub[sfo.new.micro.seurat.clean@meta.data$RNA_snn_res.0.5 %in% c(5)] = "Mac.MHClow"

sfo.new.micro.seurat.clean@meta.data$cell = rownames(sfo.new.micro.seurat.clean@meta.data)

tmp.marker.clean = FindAllMarkers(sfo.new.micro.seurat.clean,min.pct = 0.1,only.pos = T)

FeaturePlot(sfo.new.micro.seurat.clean,"Gpnmb")

saveRDS(sfo.new.micro.seurat.clean,"/Volumes/SFO/snAnalysis/new/microglia_clean_seurat_final.rds")

sfo.new.micro.seurat.clean = readRDS("/Volumes/SFO/snAnalysis/new/microglia_clean_seurat_final.rds")
sfo.new.micro.seurat.clean.clean = subset(sfo.new.micro.seurat.clean,sample %in% c("l1","l2","l3","l4","l5"))
sfo.new.micro.seurat.clean.clean  = RunUMAP(sfo.new.micro.seurat.clean.clean,dims = 1:15,reduction = 'pca')

sfo.new.micro.seurat.clean.clean = readRDS("/Volumes/SFO/The_SFO_project/data/processed/microglia_clean_seurat_final_clean5.rds")
sfo.new.micro.seurat.clean.clean@meta.data$sub[sfo.new.micro.seurat.clean.clean@meta.data$sub %in% c("DAM","pre-DAM")] = "Micro.Active"
Idents(sfo.new.micro.seurat.clean.clean) = "sub"

p1 = DimPlot(sfo.new.micro.seurat.clean.clean,pt.size = 1,cols = sn.micro.color,label = F)+NoLegend()+p.cleanumap
ggsave("/Users/ruoqing/Projects/sfo/fig_sn_merfish/sn_micro_clean_umap.png",p1,width = 6,height = 4,dpi = 300)

levels(sfo.new.micro.seurat.clean.clean) = c("Micro.Hemeo","Micro.Active","Mac.MHClow","Mac.MHChigh")

micro.marker = FindAllMarkers(sfo.new.micro.seurat.clean.clean,min.pct = 0.1,only.pos = T)
write.csv(micro.marker,"/Volumes/SFO/The_SFO_project/data/interim/DEG/SN_micro_clean_marker.csv",quote = F)
use.gene = c("Tmem119","P2ry12","Csmd3","Hexb","Itgam","Selplg","Csf1r",
             "Trem2","Ccl3","Lpl","Cxcl10","Spp1",
             "C1qa","Lyve1",
             "Cd163","F13a1","Mrc1","Ms4a7","Clec12a","Cd74","H2-Ab1","Fxyd5")
tmp = AverageExpression(sfo.new.micro.seurat.clean.clean,assays = "RNA",features = use.gene,return.seurat = T)
DoHeatmap(tmp,features = use.gene)

saveRDS(sfo.new.micro.seurat.clean.clean,"/Volumes/SFO/The_SFO_project/data/processed/microglia_clean_seurat_final_clean5.rds")

sfo.new.micro.seurat.clean.clean = readRDS("/Volumes/SFO/The_SFO_project/data/processed/microglia_clean_seurat_final_clean5.rds")

DimPlot(sfo.new.micro.seurat.clean.clean)

tmp = sfo.new.micro.seurat.clean.clean@meta.data

tmp$age = "m3"
tmp$age[tmp$sample %in% c("l2","l5")] = "m15"
tmp$age[tmp$sample %in% c("l3")] = "m24"

tmp$age = factor(tmp$age,levels = c("m3","m15","m24"),ordered = T)


library(dplyr)
library(tidyr)


sub_micro <- c("Micro.Hemeo","Micro.Active")

df_micro_mean <- tmp %>%
  filter(sub %in% sub_micro) %>%
  dplyr::count(age, sample, sub, name = "n_cell") %>%           # replicate-level count
  group_by(age, sub) %>%
  summarise(n_cell = mean(n_cell), .groups = "drop") %>% # mean per replicate
  mutate(sub = factor(sub, levels = sub_micro))

p1 <- ggplot(df_micro_mean, aes(x = age, y = n_cell, fill = sub)) +
  geom_col(position = "stack") +
  theme_classic() +
  scale_fill_manual(breaks = names(sn.micro.color), values = sn.micro.color) +
  ylab("Mean # cells per biological replicate")

sub_mac = c("Mac.MHClow","Mac.MHChigh")

df_mac_mean <- tmp %>%
  filter(sub %in% sub_mac) %>%
  dplyr::count(age, sample, sub, name = "n_cell") %>%           # replicate-level count
  group_by(age, sub) %>%
  summarise(n_cell = mean(n_cell), .groups = "drop") %>% # mean per replicate
  mutate(sub = factor(sub, levels = sub_mac))

p2 <- ggplot(df_mac_mean, aes(x = age, y = n_cell, fill = sub)) +
  geom_col(position = "stack") +
  theme_classic() +
  scale_fill_manual(breaks = names(sn.micro.color), values = sn.micro.color) +
  ylab("Mean # cells per biological replicate")

p = p1+p2
ggsave("/Users/ruoqing/Projects/sfo/fig_sn_merfish/sn_mac_micro_clean_bar.pdf",p,width = 7,height = 4)


tmp = sfo.new.micro.seurat.clean.clean@meta.data
tmp = tmp[,c("sample","sub")]
tmp = reshape2::dcast(tmp,sample ~ sub)
tmp = as.data.frame(tmp)

tmp$age = "m3"
tmp$age[tmp$sample %in% c("l2","l5")] = "m15"
tmp$age[tmp$sample %in% c("l3")] = "m24"

tmp = tmp[,c("sample","age","Mac.MHChigh","Mac.MHClow","Micro.Active","Micro.Hemeo")]
tmp = tmp[tmp$age %in% c("m3","m15"),]

write.csv(tmp,"/Volumes/SFO/The_SFO_project/data/interim/sccoda/mac_cell_num.csv",row.names = F,quote = F)

tmp[,3:6] = tmp[,3:6]/rowSums(tmp[,3:6])
tmp = reshape2::melt(tmp)
age.color = c("#ff796c","#75bbfd","#677a04")
names(age.color) = c("m3","m15","m24")
tmp$age = factor(tmp$age,levels = c("m3","m15"),ordered = T)
tmp$variable = factor(tmp$variable,levels = c("Micro.Hemeo","Micro.Active","Mac.MHClow","Mac.MHChigh"),ordered = T)

p1 = ggplot(tmp,aes(x = variable, y = value,fill = age))+
  geom_boxplot(size = 0.1)+
  scale_fill_manual(breaks = names(age.color),values = age.color)+
  theme_classic()
ggsave("/Users/ruoqing/Projects/sfo/fig_sn_merfish/micro_mac_sccoda.pdf",p1,width = 4,height = 2)

sfo.new.micro.seurat.clean.clean@meta.data$age = "m3"
sfo.new.micro.seurat.clean.clean@meta.data$age[sfo.new.micro.seurat.clean.clean@meta.data$sample %in% c("l2","l5")] = "m15"
sfo.new.micro.seurat.clean.clean@meta.data$age[sfo.new.micro.seurat.clean.clean@meta.data$sample %in% c("l3")] = "m24"

sfo.new.micro.seurat.clean.clean@meta.data$new_anno = "Microglia"
sfo.new.micro.seurat.clean.clean@meta.data$new_anno[sfo.new.micro.seurat.clean.clean@meta.data$sub %in% c("Mac.MHClow","Mac.MHChigh")] = "BAM"
Idents(sfo.new.micro.seurat.clean.clean) = "new_anno"
micro_mac_deg_aging = data.frame()
for (i in unique(sfo.new.micro.seurat.clean.clean@meta.data$new_anno)) {
  print(i)

tmp = subset(sfo.new.micro.seurat.clean.clean,new_anno == i)
tmp = subset(tmp,age %in% c("m3","m24"))
DefaultAssay(tmp) <- "RNA"
expr <- as.matrix(GetAssayData(tmp, slot = "data")) 

detected_freq <- rowSums(expr > 0) / ncol(expr)
keep_genes <- detected_freq >= 0.02
expr <- expr[keep_genes, ]

cd <- tmp@meta.data[colnames(expr), , drop = FALSE]
cd$age_month <- recode(cd$age,
                       "m3"  = 3,
#                       "m15" = 15
                       "m24" = 24
                       ) |> as.numeric()
cd$gc_z  <- scale(log1p(colSums(expr > 0))) |> as.numeric()
cd$umi_z <- scale(log1p(cd$nCount_RNA))     |> as.numeric()

sca <- FromMatrix(
  exprs = expr,
  cData = cd,
  fData = data.frame(gene = rownames(expr), row.names = rownames(expr))
)

zfit <- zlm(~ age_month + gc_z + umi_z, sca)

s  <- summary(zfit, doLRT = "age_month")
dt <- as.data.table(s$datatable)

lrt <- dt[component == "H" & contrast == "age_month",
          .(gene = primerid, pval = `Pr(>Chisq)`)]
lfc <- dt[component == "logFC" & contrast == "age_month",
          .(gene = primerid, log2FC_per_month = coef)]

deg_age <- merge(lrt, lfc, by = "gene")
deg_age$padj_bonf <- p.adjust(deg_age$pval, method = "bonferroni")
deg_age$sub = i
micro_mac_deg_aging = rbind(micro_mac_deg_aging,deg_age)

}

library(Seurat)
library(MAST)
library(dplyr)

deg_list <- list()

for (i in unique(sfo.new.micro.seurat.clean.clean@meta.data$new_anno)) {
  
  cat("\n=== ", i, " ===\n")
  
  tmp <- subset(sfo.new.micro.seurat.clean.clean, subset = new_anno == i & age %in% c("m3","m24"))
  DefaultAssay(tmp) <- "RNA"
  
  tab_age <- table(tmp$age)
  if (!all(c("m3","m24") %in% names(tab_age))) {
    cat("Skip: missing m3 or m24\n")
    next
  }
  if (tab_age["m3"] < 50 || tab_age["m24"] < 50) {
    cat("Skip: not enough cells. m3=", tab_age["m3"], " m24=", tab_age["m24"], "\n")
    next
  }
  
  set.seed(1)
  cells_m3  <- WhichCells(tmp, expression = age == "m3")
  cells_m24 <- WhichCells(tmp, expression = age == "m24")
  if (length(cells_m3)  > 2000) cells_m3  <- sample(cells_m3,  2000)
  if (length(cells_m24) > 2000) cells_m24 <- sample(cells_m24, 2000)
  tmp <- subset(tmp, cells = c(cells_m3, cells_m24))
  
  expr <- as.matrix(GetAssayData(tmp, assay = "RNA", slot = "data"))
  
  detected_freq <- rowSums(expr > 0) / ncol(expr)
  keep_genes <- detected_freq > 0.10
  expr <- expr[keep_genes, , drop = FALSE]
  if (nrow(expr) < 10) {
    cat("Skip: too few genes after filtering\n")
    next
  }
  
  cd <- tmp@meta.data[colnames(expr), , drop = FALSE]
  cd$age <- factor(cd$age, levels = c("m3","m24"))
  cd$gc_z <- as.numeric(scale(log1p(colSums(expr > 0))))
  cd$wellKey <- rownames(cd)
  
  fdat <- data.frame(
    primerid = rownames(expr),
    gene = rownames(expr),
    row.names = rownames(expr)
  )
  
  sca <- FromMatrix(exprsArray = expr, cData = cd, fData = fdat)
  zfit <- zlm(~ age + gc_z, sca)
  
  age_term <- grep("^age", colnames(zfit@coefD), value = TRUE)  # e.g. "agem24"
  if (length(age_term) != 1) {
    stop("Unexpected age term(s): ", paste(age_term, collapse = ", "))
  }
  
  s  <- summary(zfit, doLRT = age_term)
  dt <- s$datatable
  
  ptab <- dt %>%
    filter(component == "H", contrast == age_term) %>%
    transmute(gene = primerid, p_value = `Pr(>Chisq)`)
  
  ctab <- dt %>%
    filter(component == "logFC", contrast == age_term) %>%
    transmute(gene = primerid, coef = coef)
  
  res <- ptab %>%
    left_join(ctab, by = "gene") %>%
    mutate(
      p_fdr = p.adjust(p_value, method = "BH"),  # <- FDR (Benjamini-Hochberg)
      log2FC = coef / log(2),
      new_anno = i,
      n_cells_m3  = as.integer(tab_age["m3"]),
      n_cells_m24 = as.integer(tab_age["m24"])
    ) %>%
    arrange(p_fdr, p_value)
  
  deg_list[[i]] <- res
}

deg_table <- bind_rows(deg_list)

write.csv(deg_table,"/Volumes/SFO/The_SFO_project/data/interim/DEG/microglia_bam_aging_deg.csv",quote = F)


deg_table = readRDS("/Volumes/SFO/The_SFO_project/data/interim/DEG/microglia_bam_aging_deg.rds")
deg_table$sig = "none"
deg_table$sig[deg_table$p_fdr < 0.05 & deg_table$log2FC >= 0.2] = "up"
deg_table$sig[deg_table$p_fdr < 0.05 & deg_table$log2FC <= (-0.2)] = "down"

write.csv(deg_table,"/Users/ruoqing/Projects/sfo/microglia_bam_aging_deg.csv",quote = F)

#######################################################

sfo.new.oligo.seurat = readRDS("/Volumes/SFO/snAnalysis/new/oligo_clean_seurat.rds")

sfo.new.oligo.seurat@meta.data$age = ""
sfo.new.oligo.seurat@meta.data$age[sfo.new.oligo.seurat@meta.data$sample %in% c("l1","l4")] = "m3"
sfo.new.oligo.seurat@meta.data$age[sfo.new.oligo.seurat@meta.data$sample %in% c("l2","l5")] = "m15"
sfo.new.oligo.seurat@meta.data$age[sfo.new.oligo.seurat@meta.data$sample %in% c("l3","l6")] = "m24"

Idents(sfo.new.oligo.seurat) = "age"

sfo.sn.15vs3 = FindMarkers(sfo.new.oligo.seurat,ident.1 = "m15",ident.2 = "m3")
sfo.sn.24vs3 = FindMarkers(sfo.new.oligo.seurat,ident.1 = "m24",ident.2 = "m3")

sfo.sn.15vs3$gene = rownames(sfo.sn.15vs3)
sfo.sn.24vs3$gene = rownames(sfo.sn.24vs3)

sfo.sn.15vs3$diff = sfo.sn.15vs3$pct.1 - sfo.sn.15vs3$pct.2
sfo.sn.24vs3$diff = sfo.sn.24vs3$pct.1 - sfo.sn.24vs3$pct.2

sfo.sn.15vs3$sig = "non"
sfo.sn.15vs3$sig[sfo.sn.15vs3$avg_log2FC >= 0.58 & sfo.sn.15vs3$p_val_adj <= 0.05] = "up"
sfo.sn.15vs3$sig[sfo.sn.15vs3$avg_log2FC <= (-0.58) & sfo.sn.15vs3$p_val_adj <= 0.05] = "down"

sfo.sn.15vs3$label = ""
for (i in c("C4b","S100b","Ogdhl","Apoe","Serpina3n","Ighm","Klk6","Ifi27","Prkg1","Rbfox1","Iglc2")){
  sfo.sn.15vs3$label[sfo.sn.15vs3$gene == i] = i
}

sfo.sn.24vs3$sig = "non"
sfo.sn.24vs3$sig[sfo.sn.24vs3$avg_log2FC >= 0.58 & sfo.sn.24vs3$p_val_adj <= 0.05] = "up"
sfo.sn.24vs3$sig[sfo.sn.24vs3$avg_log2FC <= (-0.58) & sfo.sn.24vs3$p_val_adj <= 0.05] = "down"

sfo.sn.24vs3$label = ""
for (i in c("C4b","S100b","Ogdhl","Apoe","Serpina3n","Ighm","Klk6","Ifi27","Prkg1","Rbfox1","Igkc")){
  sfo.sn.24vs3$label[sfo.sn.24vs3$gene == i] = i
}

#oligo sub 
DimPlot(sfo.new.oligo.seurat,split.by = "sample",ncol = 3)
Idents(sfo.new.oligo.seurat) = "RNA_snn_res.0.3"

DimPlot(sfo.new.oligo.seurat,label = T,split.by = "sample",ncol = 3)
FeaturePlot(sfo.new.oligo.seurat,"C4b",split.by = "age")
DimPlot(sfo.new.oligo.seurat,label = T)

sfo.new.oligo.seurat <- FindClusters(sfo.new.oligo.seurat, resolution = 0.8)

sfo.new.oligo.seurat@meta.data$oligo_sub = ""
sfo.new.oligo.seurat@meta.data$oligo_sub[sfo.new.oligo.seurat@meta.data$RNA_snn_res.0.8 %in% c(2,4,6)] = "oligo1"
sfo.new.oligo.seurat@meta.data$oligo_sub[sfo.new.oligo.seurat@meta.data$RNA_snn_res.0.8 %in% c(0,1,8)] = "oligo2"
sfo.new.oligo.seurat@meta.data$oligo_sub[sfo.new.oligo.seurat@meta.data$RNA_snn_res.0.8 %in% c(3,5)] = "oligo3"
sfo.new.oligo.seurat@meta.data$oligo_sub[sfo.new.oligo.seurat@meta.data$RNA_snn_res.0.8 %in% c(10)] = "oligo4"
sfo.new.oligo.seurat@meta.data$oligo_sub[sfo.new.oligo.seurat@meta.data$RNA_snn_res.0.8 %in% c(9)] = "oligo5"
sfo.new.oligo.seurat@meta.data$oligo_sub[sfo.new.oligo.seurat@meta.data$RNA_snn_res.0.8 %in% c(7)] = "oligo6"

Idents(sfo.new.oligo.seurat) = "oligo_sub"
DimPlot(sfo.new.oligo.seurat,label = T,split.by = "age",ncol = 3)

oligo.markers = FindAllMarkers(sfo.new.oligo.seurat,min.pct = 0.1,only.pos = T)

FeaturePlot(sfo.new.oligo.seurat,"Bcas1")


sfo.new.oligo.seurat = readRDS("/Volumes/SFO/The_SFO_project/data/processed/sfo_sn_oligo_update.rds")
library(Seurat)
DimPlot(sfo.new.oligo.seurat)

markers = FindAllMarkers(sfo.new.oligo.seurat,min.pct = 0.1,only.pos = T)

Idents(sfo.new.oligo.seurat) = "RNA_snn_res.0.3"

DimPlot(sfo.new.oligo.seurat,split.by = "sample",label = T,ncol = 3)

sfo.new.oligo.seurat.clean = subset(sfo.new.oligo.seurat, sample %in% c("l1","l2","l3","l4","l5"))

sfo.new.oligo.seurat.clean = RunUMAP(sfo.new.oligo.seurat.clean,dims = 1:15,reduction = 'pca')

DimPlot(sfo.new.oligo.seurat.clean,split.by = "sample",label = T,ncol = 3)

VlnPlot(sfo.new.oligo.seurat.clean,"Olig2")

sfo.new.oligo.seurat.clean@meta.data$oligo_sub2 = ""
sfo.new.oligo.seurat.clean@meta.data$oligo_sub2[sfo.new.oligo.seurat.clean@meta.data$RNA_snn_res.0.3 == 0] = "MOL"
sfo.new.oligo.seurat.clean@meta.data$oligo_sub2[sfo.new.oligo.seurat.clean@meta.data$RNA_snn_res.0.3 == 3] = "MOL"
sfo.new.oligo.seurat.clean@meta.data$oligo_sub2[sfo.new.oligo.seurat.clean@meta.data$RNA_snn_res.0.3 == 7] = "MOL"
sfo.new.oligo.seurat.clean@meta.data$oligo_sub2[sfo.new.oligo.seurat.clean@meta.data$RNA_snn_res.0.3 == 6] = "COP"
sfo.new.oligo.seurat.clean@meta.data$oligo_sub2[sfo.new.oligo.seurat.clean@meta.data$RNA_snn_res.0.3 == 4] = "ARO1"
sfo.new.oligo.seurat.clean@meta.data$oligo_sub2[sfo.new.oligo.seurat.clean@meta.data$RNA_snn_res.0.3 == 1] = "ARO2"
sfo.new.oligo.seurat.clean@meta.data$oligo_sub2[sfo.new.oligo.seurat.clean@meta.data$RNA_snn_res.0.3 == 5] = "ARO2"
sfo.new.oligo.seurat.clean@meta.data$oligo_sub2[sfo.new.oligo.seurat.clean@meta.data$RNA_snn_res.0.3 == 2] = "ARO3"

Idents(sfo.new.oligo.seurat.clean) = "oligo_sub2"
DimPlot(sfo.new.oligo.seurat.clean,label = T)
FeaturePlot(sfo.new.oligo.seurat.clean,c("C4b","Serpina3n"))

saveRDS(sfo.new.oligo.seurat.clean,"/Volumes/SFO/The_SFO_project/data/processed/sfo_sn_oligo_clean_update.rds")

sfo.new.oligo.seurat.clean = readRDS("/Volumes/SFO/The_SFO_project/data/processed/sfo_sn_oligo_clean_update.rds")

sfo.new.oligo.seurat.clean@meta.data$oligo_sub2[sfo.new.oligo.seurat.clean@meta.data$oligo_sub2 == "ARO1"] = "IRO"
sfo.new.oligo.seurat.clean@meta.data$oligo_sub2[sfo.new.oligo.seurat.clean@meta.data$oligo_sub2 == "ARO2"] = "ARO1"
sfo.new.oligo.seurat.clean@meta.data$oligo_sub2[sfo.new.oligo.seurat.clean@meta.data$oligo_sub2 == "ARO3"] = "ARO2"

Idents(sfo.new.oligo.seurat.clean) = "oligo_sub2"
p1 = DimPlot(sfo.new.oligo.seurat.clean,label = T,cols = merfish.oligo.sub.cols)+p.cleanumap+NoLegend()
ggsave("/Users/ruoqing/Projects/sfo/fig_sn_merfish/sn_oligo_newdata_sub.png",width = 6,height = 4,dpi = 300)

sfo.new.oligo.seurat.clean@meta.data$cell = rownames(sfo.new.oligo.seurat.clean@meta.data)

sfo.new.oligo_opc.seurat.clean = subset(sfo.new.all.seurat,RNA_cluster %in% c("Oligo","OPC"))
sfo.new.oligo_opc.seurat.clean = subset(sfo.new.oligo_opc.seurat.clean,sample %in% c("l1","l2","l3","l4","l5"))
sfo.new.oligo_opc.seurat.clean = NormalizeData(sfo.new.oligo_opc.seurat.clean)
sfo.new.oligo_opc.seurat.clean = FindVariableFeatures(sfo.new.oligo_opc.seurat.clean,nfeatures = 2000)
sfo.new.oligo_opc.seurat.clean = ScaleData(sfo.new.oligo_opc.seurat.clean)

sfo.new.oligo_opc.seurat.clean = RunPCA(sfo.new.oligo_opc.seurat.clean)
sfo.new.oligo_opc.seurat.clean = RunUMAP(sfo.new.oligo_opc.seurat.clean,dims = 1:20,reduction = 'pca')

sfo.new.oligo_opc.seurat.clean@meta.data$cell = rownames(sfo.new.oligo_opc.seurat.clean@meta.data)
dirty.cell = setdiff(sfo.new.oligo_opc.seurat.clean@meta.data$cell[sfo.new.oligo_opc.seurat.clean@meta.data$RNA_cluster == "Oligo"],rownames(sfo.new.oligo.seurat.clean@meta.data))

DimPlot(sfo.new.oligo_opc.seurat.clean)
FeaturePlot(sfo.new.oligo_opc.seurat.clean,"Bcas1")

sfo.new.oligo_opc.seurat.clean = FindNeighbors(sfo.new.oligo_opc.seurat.clean,assay = "RNA",reduction = "pca",dims = 1:20)
sfo.new.oligo_opc.seurat.clean <- FindClusters(sfo.new.oligo_opc.seurat.clean, resolution = 0.3)
 
Idents(sfo.new.oligo_opc.seurat.clean) = "RNA_snn_res.0.3"
DimPlot(sfo.new.oligo_opc.seurat.clean,label = T)

clean.cells = c(sfo.new.oligo_opc.seurat.clean@meta.data$cell[sfo.new.oligo_opc.seurat.clean@meta.data$RNA_snn_res.0.3 == 4],
                rownames(sfo.new.oligo.seurat.clean@meta.data))
sfo.new.oligo_opc.seurat.clean = subset(sfo.new.oligo_opc.seurat.clean,cell %in% clean.cells)
sfo.new.oligo_opc.seurat.clean = RunUMAP(sfo.new.oligo_opc.seurat.clean,dims = 1:20,reduction = 'pca')

sfo.new.oligo_opc.seurat.clean@meta.data$sub = "OPC"
for (i in unique(sfo.new.oligo.seurat.clean@meta.data$oligo_sub2)){
  print(i)
  sfo.new.oligo_opc.seurat.clean@meta.data$sub[sfo.new.oligo_opc.seurat.clean@meta.data$cell %in% sfo.new.oligo.seurat.clean@meta.data$cell[sfo.new.oligo.seurat.clean@meta.data$oligo_sub2 == i]] = i
}
Idents(sfo.new.oligo_opc.seurat.clean) = "sub"
p1 = DimPlot(sfo.new.oligo_opc.seurat.clean,label = F,cols = merfish.oligo.sub.cols)+NoLegend()+p.cleanumap
ggsave("/Users/ruoqing/Projects/sfo/fig_sn_merfish/sn_oligo_opc_newdata_sub.png",width = 6,height = 4,dpi = 300)

tmp = FindAllMarkers(sfo.new.oligo_opc.seurat.clean,min.pct = 0.1,only.pos = T)
saveRDS(tmp,"/Volumes/SFO/The_SFO_project/data/interim/DEG/sn_all_oligo_opc_sub_findmarkers.rds")

library(MAST)

DefaultAssay(sfo.new.oligo_opc.seurat.clean) <- "RNA"
expr <- as.matrix(GetAssayData(sfo.new.oligo_opc.seurat.clean, slot = "data")) 

detected_freq <- rowSums(expr > 0) / ncol(expr)
keep_genes <- detected_freq >= 0.02
expr <- expr[keep_genes, ]

cd <- sfo.new.oligo_opc.seurat.clean@meta.data[colnames(expr), , drop = FALSE]
cd$age_month <- recode(cd$age,
                       "m3"  = 3,
                       "m15" = 15,
                       "m24" = 24) |> as.numeric()
cd$gc_z  <- scale(log1p(colSums(expr > 0))) |> as.numeric()
cd$umi_z <- scale(log1p(cd$nCount_RNA))     |> as.numeric()

sca <- FromMatrix(
  exprs = expr,
  cData = cd,
  fData = data.frame(gene = rownames(expr), row.names = rownames(expr))
)

zfit <- zlm(~ age_month + gc_z + umi_z, sca)

s  <- summary(zfit, doLRT = "age_month")
dt <- as.data.table(s$datatable)

lrt <- dt[component == "H" & contrast == "age_month",
          .(gene = primerid, pval = `Pr(>Chisq)`)]
lfc <- dt[component == "logFC" & contrast == "age_month",
          .(gene = primerid, log2FC_per_month = coef)]

deg_age <- merge(lrt, lfc, by = "gene")
deg_age$padj_bonf <- p.adjust(deg_age$pval, method = "bonferroni")


#####
sfo.new.all.seurat = readRDS("/Volumes/SFO/snAnalysis/new/all_clean_seurat.rds")

sfo.new.micro.seurat.clean.clean = readRDS("/Volumes/SFO/The_SFO_project/data/processed/microglia_clean_seurat_final_clean5.rds")

sfo.new.neuron.seurat = readRDS("/Volumes/SFO/snAnalysis/new/neuron_clean_seurat.rds")
sfo.new.neuron.seurat.clean = readRDS("/Volumes/SFO/snAnalysis/sfo_neuron_clean.rds")


sfo.new.all.seurat.clean = subset(sfo.new.all.seurat,sample %in% c("l1","l2","l3","l4","l5"))
DimPlot(sfo.new.all.seurat.clean)

sfo.new.all.seurat.clean@meta.data$cell = rownames(sfo.new.all.seurat.clean@meta.data)

setdiff(sfo.new.all.seurat.clean@meta.data$cell[sfo.new.all.seurat.clean@meta.data$RNA_cluster == "Neuron"],sfo.new.neuron.seurat.clean@meta.data$cell_id[sfo.new.neuron.seurat.clean@meta.data$sample %in% c("l1","l2","l3","l4","l5")])

setdiff(sfo.new.all.seurat.clean@meta.data$cell[sfo.new.all.seurat.clean@meta.data$RNA_cluster == "Micro/Mac"],sfo.new.micro.seurat.clean.clean@meta.data$cell)

sfo.new.all.seurat.clean = subset(sfo.new.all.seurat.clean,cell %in% setdiff(sfo.new.all.seurat.clean@meta.data$cell,bad.cell))
saveRDS(sfo.new.all.seurat.clean,"/Volumes/SFO/The_SFO_project/data/processed/sn_all_clean_seurat_clean5.rds")


sfo.new.all.seurat.clean@meta.data$age = ""
sfo.new.all.seurat.clean@meta.data$age[sfo.new.all.seurat.clean@meta.data$sample %in% c("l1","l4")] = "m3"
sfo.new.all.seurat.clean@meta.data$age[sfo.new.all.seurat.clean@meta.data$sample %in% c("l2","l5")] = "m15"
sfo.new.all.seurat.clean@meta.data$age[sfo.new.all.seurat.clean@meta.data$sample %in% c("l3")] = "m24"