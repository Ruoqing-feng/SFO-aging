
library(MAST)

sfo.all.seurat = readRDS("/Volumes/SFO/visume/visium10_clean.rds")

fiber <- subset(sfo.all.seurat, subset = region_name == "fiber tracts")
DefaultAssay(fiber) <- "SCT"
expr <- as.matrix(GetAssayData(fiber, slot = "data")) 

detected_freq <- rowSums(expr > 0) / ncol(expr)
keep_genes <- detected_freq >= 0.02
expr <- expr[keep_genes, ]

cd <- fiber@meta.data[colnames(expr), , drop = FALSE]
cd$age_month <- recode(cd$age,
                       "m3"  = 3,
                       "m15" = 15,
                       "m24" = 24) |> as.numeric()
cd$gc_z  <- scale(log1p(colSums(expr > 0))) |> as.numeric()
cd$umi_z <- scale(log1p(cd$nCount_SCT))     |> as.numeric()

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

saveRDS(deg_age,"/Volumes/SFO/The_SFO_project/data/interim/Visium_WM_DEG_relatetoAGE.rds")

big.region.pseudo = AggregateExpression(sfo.all.seurat, assays = "SCT", return.seurat = T, group.by = c("age", "sample", "region_name"))
big.region.pseudo$new = paste0(big.region.pseudo$region_name,"_",big.region.pseudo$age)
Idents(big.region.pseudo) <- "new"


a = "m15"
b = "m3"

a.sample = unique(sfo.all.seurat@meta.data$sample[sfo.all.seurat@meta.data$region_name == "fiber tracts" & sfo.all.seurat@meta.data$age == a])
b.sample = unique(sfo.all.seurat@meta.data$sample[sfo.all.seurat@meta.data$region_name == "fiber tracts" & sfo.all.seurat@meta.data$age == b])

a.sample.list = paste0(a,"_",a.sample,"_","fiber tracts")
b.sample.list = paste0(b,"_",b.sample,"_","fiber tracts")

tmp.mtx = big.region.pseudo@assays$SCT$counts[,c(a.sample.list,b.sample.list)]
tmp.mtx = as.matrix(tmp.mtx)

tmp.mtx <- tmp.mtx[rowSums(tmp.mtx) >= 10, , drop = FALSE]

sampleTable <- data.frame(
  sample = colnames(tmp.mtx),
  condition = c(rep(a,length(a.sample.list)),rep(b,length(b.sample.list)))  
)

dds <- DESeqDataSetFromMatrix(
  countData = tmp.mtx,
  colData = sampleTable,
  design = ~ condition
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", a, b))
tes.tab = as.data.frame(res)

fibro.track.pseudo.m15 = tes.tab
fibro.track.pseudo.m15$gene = rownames(fibro.track.pseudo.m15)

a = "m24"
a.sample = unique(sfo.all.seurat@meta.data$sample[sfo.all.seurat@meta.data$region_name == "fiber tracts" & sfo.all.seurat@meta.data$age == a])
b.sample = unique(sfo.all.seurat@meta.data$sample[sfo.all.seurat@meta.data$region_name == "fiber tracts" & sfo.all.seurat@meta.data$age == b])

a.sample.list = paste0(a,"_",a.sample,"_","fiber tracts")
b.sample.list = paste0(b,"_",b.sample,"_","fiber tracts")

tmp.mtx = big.region.pseudo@assays$SCT$counts[,c(a.sample.list,b.sample.list)]
tmp.mtx = as.matrix(tmp.mtx)

tmp.mtx <- tmp.mtx[rowSums(tmp.mtx) >= 10, , drop = FALSE]

sampleTable <- data.frame(
  sample = colnames(tmp.mtx),
  condition = c(rep(a,length(a.sample.list)),rep(b,length(b.sample.list)))  
)

dds <- DESeqDataSetFromMatrix(
  countData = tmp.mtx,
  colData = sampleTable,
  design = ~ condition
)
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", a, b))
tes.tab = as.data.frame(res)

fibro.track.pseudo.m24 = tes.tab
fibro.track.pseudo.m24$gene = rownames(fibro.track.pseudo.m24)
fibro.track.pseudo.m24[is.na(fibro.track.pseudo.m24)] = 0

fibro.track.pseudo.m15[is.na(fibro.track.pseudo.m15)] = 0


fibro.track.pseudo.m24$sig = "non"
fibro.track.pseudo.m24$sig[fibro.track.pseudo.m24$log2FoldChange >= 0.9 & fibro.track.pseudo.m24$padj <= 0.05] = "up"
fibro.track.pseudo.m24$sig[fibro.track.pseudo.m24$log2FoldChange <= (-0.9) & fibro.track.pseudo.m24$padj <= 0.05] = "down"


fibro.track.pseudo.m15$sig = "non"
fibro.track.pseudo.m15$sig[fibro.track.pseudo.m15$log2FoldChange >= 0.9 & fibro.track.pseudo.m15$padj <= 0.05] = "up"
fibro.track.pseudo.m15$sig[fibro.track.pseudo.m15$log2FoldChange <= (-0.9) & fibro.track.pseudo.m15$padj <= 0.05] = "down"

saveRDS(fibro.track.pseudo.m15,"/Volumes/SFO/The_SFO_project/data/interim/Visium_WM_DEG_m15.rds")
saveRDS(fibro.track.pseudo.m24,"/Volumes/SFO/The_SFO_project/data/interim/Visium_WM_DEG_m24.rds")


pseudo.res = unique(c(fibro.track.pseudo.m15$gene[fibro.track.pseudo.m15$sig == "up"],fibro.track.pseudo.m24$gene[fibro.track.pseudo.m24$sig == "up"]))

intersect(pseudo.res,"C4b")

deg_age$pseudo = "no"
deg_age$pseudo[deg_age$gene %in% pseudo.res] = "yes"

s = deg_age[deg_age$padj_bonf <= 0.05 & deg_age$pseudo == "yes" & deg_age$log2FC_per_month > 0,]
saveRDS(s,"/Volumes/SFO/The_SFO_project/data/interim/aging_deg_3time_visium2_274.rds")
write.table(s,"/Volumes/SFO/The_SFO_project/data/interim/aging_deg_3time_visium2.txt",quote = F,row.names = F)

s = readRDS("/Volumes/SFO/The_SFO_project/data/interim/aging_deg_3time_visium2_274.rds")
pseudo.res.down = unique(c(fibro.track.pseudo.m15$gene[fibro.track.pseudo.m15$sig == "down"],fibro.track.pseudo.m24$gene[fibro.track.pseudo.m24$sig == "down"]))
deg_age$gene[deg_age$gene %in% pseudo.res.down] = "yes_down"

s2 = deg_age[deg_age$padj_bonf <= 0.05 & deg_age$pseudo == "yes_down" & deg_age$log2FC_per_month < 0,]

sss = union(fibro.track.pseudo.m15$gene[fibro.track.pseudo.m15$sig != "down"],
            fibro.track.pseudo.m24$gene[fibro.track.pseudo.m24$sig != "down"])

tmp = subset(sfo.all.seurat,region_name == "fiber tracts")
tmp = SCTransform(tmp,assay = "Spatial")
Idents(tmp) = "age"

tmp = AverageExpression(tmp,assays = "SCT",features = s$gene,return.seurat = T)
DoHeatmap(tmp,features = sss)

tmp = tmp@assays$SCT$scale.data
bk = seq(-2,2,0.1)
p1 = pheatmap::pheatmap(tmp,breaks = bk,
                        color = colorRampPalette(c("#5C88DA99","white", "#CC0C0099"))(length(bk)),
                        cluster_cols = F,cluster_rows = T,border_color = NA)

mat <- tmp 
thr <- 0.1  

white_idx   <- apply(mat, 1, function(x) max(abs(x)) < thr)
white_genes <- rownames(mat)[white_idx]

length(white_genes)
head(white_genes)
mat_filt <- mat
p1 = pheatmap::pheatmap(mat_filt,breaks = bk,
                        color = colorRampPalette(c("#5C88DA99","white", "#CC0C0099"))(length(bk)),
                        cluster_cols = F,cluster_rows = T,border_color = NA)

intersect(white_genes,s$gene)

table(s$tmp)

s$tmp = -log10(s$padj_bonf+1e-300)
s = s[order(s$tmp,decreasing = T),]


use.term = c("inflammatory response",
             "immune effector process",
             "leukocyte activation",
             "antigen processing and presentation",
             "Tyrobp causal network in microglia",
             "Microglia pathogen phagocytosis pathway",
             "Phagosome - Mus musculus (house mouse)",
             "positive regulation of cytokine production")

all.terms = read.delim("/Volumes/SFO/The_SFO_project/data/external/metascape_result.txt",header = T)
#setdiff(use.term,all.terms$Description)
all.terms = all.terms[all.terms$Description %in% use.term,]
all.terms = all.terms[!duplicated(all.terms$Description),]

all.terms$Log.q.value. = c(30.811,28.048,22.092,21.959,19.753,19.029,18.761,15.645)
all.terms = all.terms[order(all.terms$Log.q.value.,decreasing = T),]
all.terms$Description = factor(all.terms$Description,levels = all.terms$Description,ordered = T)


laura.selected.region = readRDS("/Volumes/SFO/The_SFO_project/data/processed/Visium_laura_selected_regions6.rds")
scanpy.meta = read.csv("/Volumes/SFO/The_SFO_project/data/external/Visium_scanpy_meta.csv",header = T,row.names = 1)

tmp = scanpy.meta[rownames(laura.selected.region@meta.data),]

laura.selected.region@meta.data$Aging_score2 = tmp$Aging_score2
tmp = laura.selected.region@meta.data[,c("age","sample","selectedregion","Aging_score2")]
tmp = tmp[tmp$selectedregion %in% c("cc","fimbria","vhc"),]

s = matrix(100,ncol = length(unique(tmp$selectedregion)),nrow = length(unique(tmp$sample)))
colnames(s) = c("cc","fimbria","vhc")
rownames(s) = unique(tmp$sample)

for (i in colnames(s)){
  for (j in rownames(s)){
    s[j,i] = length(tmp$Aging_score2[tmp$sample == j & tmp$selectedregion == i & tmp$Aging_score2 >= 1])/length(tmp$Aging_score2[tmp$sample == j & tmp$selectedregion == i])
  }
}

s = as.data.frame(s)
s$sample = rownames(s)
s$age = "m3"
s$age[s$sample %in% c("r151","r152","n15m1")] = "m15"
s$age[s$sample %in% c("r241","r243","n22m1")] = "m24"

s = reshape2::melt(s)
s$age = factor(s$age,levels = c("m3","m15","m24"),ordered = T)
s$variable = factor(s$variable,levels = c("cc","vhc","fimbria"),ordered = T)

colnames(s)[4] = "ratio"

s2 = s[,c("age","variable")]
s2 = s2[!duplicated(s2),]
s2$ratio = 100
s2$means = 100

for (i in unique(s2$age)){
  for (j in unique(s2$variable)){
    s2$ratio[s2$age == i & s2$variable == j] = length(tmp$Aging_score2[tmp$age == i & tmp$selectedregion == j & tmp$Aging_score2 >= 1])/length(tmp$Aging_score2[tmp$age == i & tmp$selectedregion == j])
    s2$means[s2$age == i & s2$variable == j] = mean(tmp$Aging_score2[tmp$age == i & tmp$selectedregion == j])
  }
}
