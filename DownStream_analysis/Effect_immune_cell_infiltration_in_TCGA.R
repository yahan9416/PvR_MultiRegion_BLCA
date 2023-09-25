#First load the estimated immune cell infiltration by TIMER website
estimation_tcga=read.table("./Survival_cohorts/immuneEstimation.txt",header = TRUE,sep="\t",row.names = 1)
patient_order=substring(rownames(estimation_tcga),1,12)
index=match(Survival_info[which(Survival_info$cancer == "BLCA"),1],patient_order)

estimation_blca=estimation_tcga[na.omit(index),]
#load the differential expression genes
marker_gene=read.table("./Bladder_cancer/ImmStromal_analysis/Stromal/CCA/BLCA_reBuild_Stromal_DEG_Celltype.txt",header = TRUE,sep="\t")

CellSig_CTL<<-NULL
batch_test_influence_Immune_infiltration<-function(clus){
  temp_marker_gene=marker_gene[which(marker_gene$cluster == clus),]
  DEgenes_cluster=temp_marker_gene$gene[order(temp_marker_gene$avg_logFC,decreasing = T)]
  
  if(length(DEgenes_cluster)>200){
    DEgenes_cluster=DEgenes_cluster[1:200]
  }else{
    DEgenes_cluster=DEgenes_cluster
  }
  cancer_expMat=t(expr[["BLCA"]])
  inter_gene=na.omit(intersect(DEgenes_cluster,rownames(cancer_expMat)))
  gsva_score=gsva(cancer_expMat,list(inter_gene))
  
  index=match(rownames(estimation_blca),colnames(gsva_score))
  temp_celltype=gsva_score[index]
  result=apply(as.matrix(estimation_blca),2,function(x){
    unlist(cor.test(x,temp_celltype))[c(3,4)]
  })
  result=as.vector(result)
  result=c(clus,result)
  names(result)=c("Celltype",paste(rep(colnames(estimation_blca),each=2),c("Pvalue","Cor_eff"),sep="|"))
  CellSig_CTL<<-rbind(CellSig_CTL,result)
}
result=apply(matrix(unique(marker_gene$cluster)),1,batch_test_influence_Immune_infiltration)


write.table(CellSig_CTL,file="Corr_TimerEsti_Imm_Stromal_Signature.txt",col.names = TRUE,sep="\t",quote = FALSE)
saveRDS(CellSig_CTL,file="Corr_TimerEsti_Imm_Stromal_Signature.rds")



temp_marker_gene=marker_gene[which(marker_gene$cluster == "Endo_ACKR1"),]
DEgenes_cluster=temp_marker_gene$gene[order(temp_marker_gene$avg_logFC,decreasing = T)]

if(length(DEgenes_cluster)>200){
  DEgenes_cluster=DEgenes_cluster[1:200]
}else{
  DEgenes_cluster=DEgenes_cluster
}
cancer_expMat=t(expr[["BLCA"]])
inter_gene=na.omit(intersect(DEgenes_cluster,rownames(cancer_expMat)))
gsva_score=gsva(cancer_expMat,list(inter_gene))

index=match(rownames(estimation_blca),colnames(gsva_score))
estimation_blca$Endo_ACKR1=gsva_score[1,index]

p=ggplot(estimation_blca, aes(x=Fibro_CCL5, y=Dendritic)) +
  geom_point() + 
  geom_smooth(method=lm,fullrange=TRUE)+theme_bw()
ggsave(file="Endo_ACKR1_signature_Timer_Dendritic_correlation.pdf",p,width = 5,height = 5)

