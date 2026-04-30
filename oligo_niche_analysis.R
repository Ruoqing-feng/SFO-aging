
sfo.merfish.qc.keep.update = readRDS("/Volumes/SFO/The_SFO_project/data/processed/sfo_merfish_cross_all_anno_banbksy_update_cc.rds")
merfish.all.micro_mac.clean = readRDS("/Volumes/SFO/The_SFO_project/data/processed/merfish_all_mac_micro_updated_clean.rds")


big.nerigbor = data.frame()
for (i in unique(sfo.merfish.qc.keep.update@meta.data$sample)){
  print(i)
  tmp = sfo.merfish.qc.keep.update@meta.data[sfo.merfish.qc.keep.update@meta.data$sample == i,]
  tmp = tmp[,c("x","y")]
  knn_result <- nn2(tmp, k = 51)
  knn_cell_ids <- apply(knn_result$nn.idx, 2, function(idx) rownames(tmp)[idx])
  colnames(knn_cell_ids) = 1:51
  rownames(knn_cell_ids) = rownames(tmp)
  
  big.nerigbor = rbind(big.nerigbor,knn_cell_ids)
}

big.nerigbor = big.nerigbor[,c(2:51)]
big.nerigbor$centercell = rownames(big.nerigbor)
saveRDS(big.nerigbor,"/Volumes/SFO/The_SFO_project/data/interim/all_neighbor_analysis_sections_updated.rds")

big.nerigbor.long = reshape2::melt(big.nerigbor,"centercell")
colnames(big.nerigbor.long) = c("centercell","position","cell")
big.nerigbor.long$range[big.nerigbor.long$position %in% c(2:11)] = "neighbor10_1"
big.nerigbor.long$range[big.nerigbor.long$position %in% c(12:21)] = "neighbor10_2"
big.nerigbor.long$range[big.nerigbor.long$position %in% c(22:31)] = "neighbor10_3"
big.nerigbor.long$range[big.nerigbor.long$position %in% c(32:41)] = "neighbor10_4"
big.nerigbor.long$range[big.nerigbor.long$position %in% c(42:51)] = "neighbor10_5"

big.nerigbor.long$sample = ""
for (i in unique(sfo.merfish.qc.keep.update@meta.data$sample)){
  print(i)
  big.nerigbor.long$sample[big.nerigbor.long$centercell %in% sfo.merfish.qc.keep.update@meta.data$cell[sfo.merfish.qc.keep.update@meta.data$sample == i]] = i
}

big.nerigbor.long$sample = ""
for (i in unique(sfo.merfish.qc.keep.update@meta.data$sample)){
  print(i)
  big.nerigbor.long$sample[big.nerigbor.long$centercell %in% sfo.merfish.qc.keep.update@meta.data$cell[sfo.merfish.qc.keep.update@meta.data$sample == i]] = i
}

big.nerigbor.long$age = ""
big.nerigbor.long$age[big.nerigbor.long$sample %in% c("sfo2r3","sfo2r6","sfo4r3","sfo4r6")] = "m3"
big.nerigbor.long$age[big.nerigbor.long$sample %in% c("sfo1r3","sfo2r2","sfo2r5","sfo4r2","sfo4r5")] = "m15"
big.nerigbor.long$age[big.nerigbor.long$sample %in% c("sfo1r2","sfo4r1","sfo4r4")] = "m24"

big.nerigbor.long$region_new_cc_new = ""
for (i in unique(sfo.merfish.qc.keep.update@meta.data$region_new_cc_new)){
  print(i)
  big.nerigbor.long$region_new_cc_new[big.nerigbor.long$centercell %in% sfo.merfish.qc.keep.update@meta.data$cell[sfo.merfish.qc.keep.update@meta.data$region_new_cc_new == i]] = i
}

merfish.all.oligo.clean = readRDS("/Volumes/SFO/The_SFO_project/data/processed/merfish_all_oligo_sub_clean.rds")
Idents(merfish.all.oligo.clean) = "oligo_sub"

Oligo.cells = merfish.all.oligo.clean@meta.data$cell
big.nerigbor.long.oligo = big.nerigbor.long[big.nerigbor.long$centercell %in% Oligo.cells,]
big.nerigbor.long.oligo = big.nerigbor.long.oligo[big.nerigbor.long.oligo$region_new_cc_new %in% c("cc","fimbria","SFO"),]

big.nerigbor.long.oligo$micro_sub = "others"
for (i in unique(merfish.all.micro_mac.clean@meta.data$sub_anno)){
  big.nerigbor.long.oligo$micro_sub[big.nerigbor.long.oligo$cell %in% merfish.all.micro_mac.clean@meta.data$cell[merfish.all.micro_mac.clean@meta.data$sub_anno == i]] = i
}

big.nerigbor.long.oligo$oligo_sub = "others"
for (i in unique(merfish.all.oligo.clean@meta.data$oligo_sub)){
  big.nerigbor.long.oligo$oligo_sub[big.nerigbor.long.oligo$centercell %in% merfish.all.oligo.clean@meta.data$cell[merfish.all.oligo.clean@meta.data$oligo_sub == i]] = i
}

big.nerigbor.long.oligo$oligo_sub2[big.nerigbor.long.oligo$oligo_sub %in% c("MOL","COP")] = "Oligo"
big.nerigbor.long.oligo$oligo_sub2[big.nerigbor.long.oligo$oligo_sub %in% c("ARO1","ARO2","ARO3")] = "ARO"

big.nerigbor.long.oligo$new = paste0(big.nerigbor.long.oligo$oligo_sub2,"_",big.nerigbor.long.oligo$age,"_",big.nerigbor.long.oligo$region_new_cc_new)

ggplot(big.nerigbor.long.oligo[big.nerigbor.long.oligo$micro_sub != "others",],aes(x = new, fill = micro_sub))+
  geom_bar(stat = "count",position = "fill")+
  theme_classic()+theme(axis.text = element_text(angle = 90))

big.nerigbor.long.oligo$major = "others"
for (i in unique(sfo.merfish.qc.keep.update@meta.data$celltype_clean_big)){
  big.nerigbor.long.oligo$major[big.nerigbor.long.oligo$cell %in% sfo.merfish.qc.keep.update@meta.data$cell[sfo.merfish.qc.keep.update@meta.data$celltype_clean_big == i]] = i
}

ggplot(big.nerigbor.long.oligo,aes(x = new2, fill = major))+
  geom_bar(stat = "count",position = "fill")+
  theme_classic()+theme(axis.text = element_text(angle = 90))

head(big.nerigbor.long.oligo)


tmp = matrix(100,ncol = length(unique(big.nerigbor.long.oligo$major)),nrow = length(unique(big.nerigbor.long.oligo$new2)))
colnames(tmp) = unique(big.nerigbor.long.oligo$major)
rownames(tmp) = unique(big.nerigbor.long.oligo$new2)

for (i in rownames(tmp)){
  print(i)
  n = nrow(big.nerigbor.long.oligo[big.nerigbor.long.oligo$new2 == i,])
  for (j in colnames(tmp)){
    m = nrow(big.nerigbor.long.oligo[big.nerigbor.long.oligo$new2 == i & big.nerigbor.long.oligo$major == j,])
    tmp[i,j] = m/n
  }
}

tmp = t(tmp)
tmp2 = tmp/tmp[,5]



tmp3 = matrix(100,ncol = length(unique(big.nerigbor.long.oligo$major)),nrow = length(unique(big.nerigbor.long.oligo$new)))
colnames(tmp3) = unique(big.nerigbor.long.oligo$major)
rownames(tmp3) = unique(big.nerigbor.long.oligo$new)

for (i in rownames(tmp3)){
  print(i)
  n = nrow(big.nerigbor.long.oligo[big.nerigbor.long.oligo$new == i,])
  for (j in colnames(tmp3)){
    m = nrow(big.nerigbor.long.oligo[big.nerigbor.long.oligo$new == i & big.nerigbor.long.oligo$major == j,])
    tmp3[i,j] = m/n
  }
}

tmp3 = t(tmp3)
tmp4 = tmp3/tmp3[,13]


big.nerigbor = readRDS("/Volumes/SFO/The_SFO_project/data/interim/all_neighbor_analysis_sections_updated.rds")

big.nerigbor.long = reshape2::melt(big.nerigbor,"centercell")
colnames(big.nerigbor.long) = c("centercell","position","cell")

merfish.all.oligo.clean = readRDS("/Volumes/SFO/The_SFO_project/data/processed/merfish_all_oligo_sub_clean.rds")
Idents(merfish.all.oligo.clean) = "oligo_sub"

big.nerigbor.long$sample = ""
for (i in unique(sfo.merfish.qc.keep.update@meta.data$sample)){
  print(i)
  big.nerigbor.long$sample[big.nerigbor.long$centercell %in% sfo.merfish.qc.keep.update@meta.data$cell[sfo.merfish.qc.keep.update@meta.data$sample == i]] = i
}

big.nerigbor.long$sample = ""
for (i in unique(sfo.merfish.qc.keep.update@meta.data$sample)){
  print(i)
  big.nerigbor.long$sample[big.nerigbor.long$centercell %in% sfo.merfish.qc.keep.update@meta.data$cell[sfo.merfish.qc.keep.update@meta.data$sample == i]] = i
}

big.nerigbor.long$age = ""
big.nerigbor.long$age[big.nerigbor.long$sample %in% c("sfo2r3","sfo2r6","sfo4r3","sfo4r6")] = "m3"
big.nerigbor.long$age[big.nerigbor.long$sample %in% c("sfo1r3","sfo2r2","sfo2r5","sfo4r2","sfo4r5")] = "m15"
big.nerigbor.long$age[big.nerigbor.long$sample %in% c("sfo1r2","sfo4r1","sfo4r4")] = "m24"

big.nerigbor.long$region_new_cc_new = ""
for (i in unique(sfo.merfish.qc.keep.update@meta.data$region_new_cc_new)){
  print(i)
  big.nerigbor.long$region_new_cc_new[big.nerigbor.long$centercell %in% sfo.merfish.qc.keep.update@meta.data$cell[sfo.merfish.qc.keep.update@meta.data$region_new_cc_new == i]] = i
}

Oligo.cells = merfish.all.oligo.clean@meta.data$cell
big.nerigbor.long.oligo = big.nerigbor.long[big.nerigbor.long$centercell %in% Oligo.cells,]
big.nerigbor.long.oligo = big.nerigbor.long.oligo[big.nerigbor.long.oligo$region_new_cc_new %in% c("cc","fimbria","SFO"),]

big.nerigbor.long.oligo$micro_sub = "others"
for (i in unique(merfish.all.micro_mac.clean@meta.data$sub_anno)){
  big.nerigbor.long.oligo$micro_sub[big.nerigbor.long.oligo$cell %in% merfish.all.micro_mac.clean@meta.data$cell[merfish.all.micro_mac.clean@meta.data$sub_anno == i]] = i
}

big.nerigbor.long.oligo$oligo_sub = "others"
for (i in unique(merfish.all.oligo.clean@meta.data$oligo_sub)){
  big.nerigbor.long.oligo$oligo_sub[big.nerigbor.long.oligo$centercell %in% merfish.all.oligo.clean@meta.data$cell[merfish.all.oligo.clean@meta.data$oligo_sub == i]] = i
}

big.nerigbor.long.oligo$oligo_sub2[big.nerigbor.long.oligo$oligo_sub %in% c("MOL","COP")] = "Oligo"
big.nerigbor.long.oligo$oligo_sub2[big.nerigbor.long.oligo$oligo_sub %in% c("ARO1","ARO2","ARO3")] = "ARO"


big.nerigbor.long.oligo$new = paste0(big.nerigbor.long.oligo$oligo_sub2,"_",big.nerigbor.long.oligo$age,"_",big.nerigbor.long.oligo$region_new_cc_new)

tmp = reshape2::dcast(big.nerigbor.long.oligo,centercell ~ micro_sub)

tmp$oligo_sub = "others"
for (i in unique(merfish.all.oligo.clean@meta.data$oligo_sub)){
  tmp$oligo_sub[tmp$centercell %in% merfish.all.oligo.clean@meta.data$cell[merfish.all.oligo.clean@meta.data$oligo_sub == i]] = i
}

tmp$oligo_sub2[tmp$oligo_sub %in% c("MOL","COP")] = "Oligo"
tmp$oligo_sub2[tmp$oligo_sub %in% c("ARO1","ARO2","ARO3")] = "ARO"

tmp = tmp[,-7]
tmp$sample = ""
for (i in unique(sfo.merfish.qc.keep.update@meta.data$sample)){
  print(i)
  tmp$sample[tmp$centercell %in% sfo.merfish.qc.keep.update@meta.data$cell[sfo.merfish.qc.keep.update@meta.data$sample == i]] = i
}

tmp$age = ""
tmp$age[tmp$sample %in% c("sfo2r3","sfo2r6","sfo4r3","sfo4r6")] = "m3"
tmp$age[tmp$sample %in% c("sfo1r3","sfo2r2","sfo2r5","sfo4r2","sfo4r5")] = "m15"
tmp$age[tmp$sample %in% c("sfo1r2","sfo4r1","sfo4r4")] = "m24"

tmp$region_new_cc_new = ""
for (i in unique(sfo.merfish.qc.keep.update@meta.data$region_new_cc_new)){
  print(i)
  tmp$region_new_cc_new[tmp$centercell %in% sfo.merfish.qc.keep.update@meta.data$cell[sfo.merfish.qc.keep.update@meta.data$region_new_cc_new == i]] = i
}

sfo.merfish.qc.keep.update@meta.data$micro_sub = "others"
for (i in unique(merfish.all.micro_mac.clean@meta.data$sub_anno)){
  sfo.merfish.qc.keep.update@meta.data$micro_sub[sfo.merfish.qc.keep.update@meta.data$cell %in% merfish.all.micro_mac.clean@meta.data$cell[merfish.all.micro_mac.clean@meta.data$sub_anno == i]] = i
}

for (i in unique(merfish.all.oligo.clean@meta.data$oligo_sub)){
  sfo.merfish.qc.keep.update@meta.data$micro_sub[sfo.merfish.qc.keep.update@meta.data$cell %in% merfish.all.oligo.clean@meta.data$cell[merfish.all.oligo.clean@meta.data$oligo_sub == i]] = i
}

sfo.merfish.qc.keep.update@meta.data$micro_sub[sfo.merfish.qc.keep.update@meta.data$micro_sub %in% c("MOL","COP")] = "Oligo"
sfo.merfish.qc.keep.update@meta.data$micro_sub[sfo.merfish.qc.keep.update@meta.data$micro_sub %in% c("ARO1","ARO2","ARO3")] = "ARO"

Idents(sfo.merfish.qc.keep.update) = "micro_sub"

the_niche = subset(sfo.merfish.qc.keep.update,region_new_cc_new %in% c("cc","fimbria","SFO"))
Idents(the_niche) = "micro_sub"
ImageDimPlot(the_niche,fov = "sfo4r5",cols = merfish.micro.color,dark.background = F,border.size = 0,axes = T,boundaries = "centroids",size = 1)+theme_classic()

#ARO_m15_fimbria
i = "3602320900024100685"
use.cell = c(i,big.nerigbor.long.oligo$cell[big.nerigbor.long.oligo$centercell == i])

the_niche = subset(sfo.merfish.qc.keep.update,cell %in% use.cell)
the_niche@meta.data$micro_sub[the_niche@meta.data$cell == i] = "centercell"
Idents(the_niche) = "micro_sub"
p1 = ImageDimPlot(the_niche,cols = merfish.micro.color,dark.background = F,border.size = 0,axes = T)+
  theme_classic()+xlab("")+ylab("")+NoLegend()
ggsave("/Users/ruoqing/Projects/sfo/fig_sn_merfish/fig4_show_case.pdf",p1,width = 6,height = 6)
p1 = ImageDimPlot(the_niche,cols = merfish.micro.color,dark.background = F,border.size = 0,axes = T)+
  theme_classic()+xlab("")+ylab("")
ggsave("/Users/ruoqing/Projects/sfo/fig_sn_merfish/fig4_show_case_legend.pdf",p1,width = 6,height = 6)