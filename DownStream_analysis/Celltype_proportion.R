library(ggsci)
library(ggplot2)
library(MAESTRO)
library(Seurat)
SeuratObj=readRDS("./Bladder_cancer/ImmStromal_analysis/reBulid_Seurat/Immune/Lymphocyte/BLCA_reBuild_Lymphocyte_seurat.rds")
SeuratObj@meta.data$curated_anno=droplevels(SeuratObj@meta.data$curated_anno)
cluster_sample=table(SeuratObj@meta.data$curated_anno,SeuratObj@meta.data$Sample)
source_pro=t(apply(cluster_sample,2,function(x) x/sum(x)))

data=data.frame(Celltype=rep(colnames(source_pro),each=10),Percentage=as.vector(source_pro),Sample=rep(rownames(source_pro),13))
data$Celltype=factor(data$Celltype,levels = c("B_MS4A1","CD8Teff_GZMA","CD8Teff_GZMB","CD8Teff_PRKCH","NK_KLRD1","Plasma_IGHG2","Plasma_IGLC2","T_CXCL13","T_IL7R","Tnaive_TCF7","Tprolif_MKI67","Treg_FOXP3","Treg_TNFRSF4"))
color_man=c("#276D9F","#2F8AC4", "#A0C9D9","#52BCA3","#ACD48A",  "#24796C","#CC61B0", "#F5BC6E", "#EA945A","#F3D4F4", "#BCBD22", "#C7AED5",  "#726BAE")

p=ggplot(data=data, aes(x=Sample, y=Percentage, fill=Celltype)) +
  geom_bar(stat="identity",alpha=0.8)+theme_classic()+theme(axis.text.x = element_text(angle = 90),legend.position ="bottom" )+scale_fill_manual(values=color_man)

ggsave(file.path("BLCA_Lymphocyte_celltype_sample_propor.pdf"),p,height = 8,width = 6)

SeuratObj@meta.data$Source="Primary"

library(ggpubr)

p=ggplot(data, aes(x=Celltype, y=Percentage, color=Source,fill=Source)) + geom_boxplot(position=position_dodge(0.85),outlier.size =0,outlier.color ="white",alpha=0.9)+ theme_bw()+theme(axis.text.x = element_text(angle = 90,vjust = 0.5),panel.background = element_rect( colour = "white"),legend.position = "top")+ scale_color_manual(values = rev(c("#518496","#facf5a")))+scale_fill_manual(values = rev(c("#518496","#facf5a")))+stat_compare_means( aes(label = ..p.signif..))
#+scale_x_discrete(limits=c("T_naive","CD8T_naive","CD8T_mem","CD8T_em","CD8T_exh","CD8T_ISG","CD8T_NK","NK","Tprolif","CD4T_naive","Th1","Th17","Tfh","CD4T_mem","Treg","T_B","T_Mono"))
ggsave(file="Lymphocyte_distribution_across_source.pdf",p,width = 8,height = 4)
