Inflammation=unique(c("VIM","FAP","COL3A1","DES","IL6","CXCL12","CXCL1","IGF1","FIGF","PDGFD","CXCL12","CXCL13","PDPN","CXCL12","CXCL14","PDGFRA","CXCL12","CXCL14","CXCL1","CXCL2","PDGFRA"))#Inflammation AND chemotaxis
Immu_reg=unique(c("C5","C5AR1","TGFB1","TGFBR1","TGFB2","TGFBR1","CXCL12","CXCR4","CXCL13","CXCR5","C7","CFB","CFH","CFI","BCAM","F11R","IRF5","CCL2","IL6","CXCL2","CXCL1","CXCL3"))
AP=unique(c("CD74","HLA-DQA1","CD83","HLA-DRA","HLA-DPA1","HLA-DQA","CD74","HLA-DRA","HLA-DRB1","H2-Q4"))
ECM=unique(c("ACTA2","TAGLN","BGN","COL8A1","COL15A1","IGFBP7","TPM1","TPM2","COL10A1","POSTN","MYL9,COL13A1","COL14A1","ACTA2","TAGLN","MYH11","MYLK","ACTG2","POSTN","FN1","LUM","DCN","VCAN","COL5A1","COL5A2","COL63A","ACTA2","TAGLN and PDGFA","MMP2","DCN","COL1A2","RGS5","MYL9","MYH11"))
Signature_list=read.table("./Pancancer_MacroMono_signature_score.txt",header=F,sep="\t")
colnames(Signature_list)=c("Symbol","Signature")
Signature_list$Symbol=gsub(" ","",Signature_list$Symbol)
Angiogenesis=Signature_list$Symbol[which(Signature_list$Signature == "Angiogenesis" )]
Tip_sig=unique(c(c("ANGPT2","APLN","FSCN1","PGF","PLXND1","ADM","PDGFB","CXCR4"),c("ADM","ANGPT2","APLN","CXCR4","ESM1","PXDN","PGF","LXN","ANGPTL2","INSR")))


Gene_list=c(list(Inflammation),list(Immu_reg),list(AP),list(Angiogenesis),list(ECM),list(Tip_sig))
names(Gene_list)=c("Inflammation","Immu_reg","AP","Angiogenesis","ECM","Tip")


DefaultAssay(SeuratObj)<-"RNA"
SeuratObj <- AddModuleScore(
  object = SeuratObj,
  features = Gene_list,
  assay	="RNA",
  ctrl = 5,
  name = names(Gene_list)
)
colnames(SeuratObj@meta.data)[13:18]=names(Gene_list)

library(ggpubr)
library(ggplot2)
library(ggsci)


Sign_score=SeuratObj@meta.data
Sign_score_source_celltype<<-NULL
get_celltype_AverExp<-function(x){
  temp_data=aggregate(Sign_score[,x]~curated_anno,Sign_score,mean)
  temp_data$Cell_type=x
  colnames(temp_data)=c("Celltype","Signature_score","Signature")
  Sign_score_source_celltype<<-rbind(Sign_score_source_celltype,temp_data)
}
result=apply(matrix(names(Gene_list)),1,get_celltype_AverExp)


cell_type_score<<-NULL
batch_generate_celltype<-function(x){
  temp_score=Sign_score_source_celltype$Signature_score[which(Sign_score_source_celltype$Celltype ==x)]
  cell_type_score<<-rbind(cell_type_score,temp_score)
  
}
result=apply(matrix(unique(Sign_score_source_celltype$Celltype)),1,batch_generate_celltype)

colnames(cell_type_score)=names(Gene_list)
rownames(cell_type_score)=unique(Sign_score_source_celltype$Celltype)


library(ComplexHeatmap)
library(pheatmap)
library(circlize)
range(cell_type_score)
pdf("cell_type_signateru.pdf",height = 5.5,width = 5)
col_fun = colorRamp2(c( 0,2), c("white", "#276D9F"))
pheatmap(cell_type_score,cluster_rows = FALSE,cluster_cols=FALSE,scale="none",border_color=NA,color =c( rep("#461D64",20),colorRampPalette(colors = c("#461D64","#299888","#E7E331"))(50),rep("#E7E331",20)))
dev.off()
