library("survival")
library("survminer")
library(GSVA)

load('./Survival_cohorts/TCGAexpr.RData')
load("./Survival_cohorts/TCGAclin.RData")
selecte_cancertype="BLCA"
selected_index=match(selecte_cancertype,names(expr))
expr=expr[selected_index]
Survival_info=clin
marker_gene=read.table("./Bladder_cancer/ImmStromal_analysis/Myeloid/CCA/BLCA_Myeloid_DEG_Celltype.txt",header = TRUE,sep="\t")
#load the differential expressed genes

Survial_cancer_celltype<<-NULL
cluster_DEgene_effect_one_dataset<-function(clus){
  print(clus)
  temp_marker_gene=marker_gene[which(marker_gene$cluster == clus),]
  DEgenes_cluster=temp_marker_gene$gene[order(temp_marker_gene$avg_logFC,decreasing = T)]
  
  if(length(DEgenes_cluster)>50){
    DEgenes_cluster=DEgenes_cluster[1:50]
  }else{
    DEgenes_cluster=DEgenes_cluster
  }
  
  
  cancer_expMat=t(expr[["BLCA"]])
  inter_gene=na.omit(intersect(DEgenes_cluster,rownames(cancer_expMat)))
  gsva_score=gsva(cancer_expMat,list(inter_gene))
  
  
  colnames(gsva_score)=substring(colnames(gsva_score),1,12)
  survival_infor=Survival_info
  length(intersect(survival_infor$patient,colnames(gsva_score)))
  sample=intersect(survival_infor$patient,colnames(gsva_score))
  index1=match(sample,survival_infor$patient)
  index2=match(sample,colnames(gsva_score))
  #print(head(survival_infor))
  survival_gsva=cbind(survival_infor[index1,],gsva_score[index2])
  #survival_infor=cbind(survival_infor,gsva_es)
  colnames(survival_gsva)[10]="GSVA_score"
  survival_gsva$GSVA_score=as.numeric(as.vector(survival_gsva$GSVA_score))
  #survival_gsva$OS[which(survival_gsva$OS > 1825)]=1825
  res.cox <- coxph(Surv(OS, EVENT) ~ GSVA_score, data = survival_gsva)
  temp_result=c(clus,"BLCA",summary(res.cox)[[7]],summary(res.cox)[[9]][3])
  Survial_cancer_celltype<<-rbind(Survial_cancer_celltype,temp_result)
  
  survival_gsva$GSVA_score_class="High"
  survival_gsva$GSVA_score_class[which(survival_gsva$GSVA_score < median(survival_gsva$GSVA_score))]="Low"
  fit<-survfit(Surv(OS,EVENT)~GSVA_score_class, data=survival_gsva)
  surv_diff<-NULL
  surv_diff <- survdiff(Surv(OS,EVENT)~GSVA_score_class, data=survival_gsva)
  p.value <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)
  
  ggsurv <- ggsurvplot(fit,palette = "jco", pval = FALSE,risk.table = T)+ggtitle(paste0("TCGA_BLCA: ",clus,"P_value=",format(p.value,digits = 3)))
  setwd("./Bladder_cancer/ImmStromal_analysis/Lymphocyte/CCA/Survival")
  p=ggarrange(ggsurv$plot, ggsurv$table, heights = c(1.5, 0.8),ncol = 1, nrow = 2)
  ggsave(paste0("TCGA_BLCA",clus,"_top50_Correaltion.pdf"),p,width = 7,height = 10.5)
}
result=apply(as.matrix(unique(as.vector(unlist(marker_gene$cluster)))),1,cluster_DEgene_effect_one_dataset)

colnames(Survial_cancer_celltype)=c("Celltype","Cancertype","coef", "HR","se(coef)", "Z_score", "Pr(>|z|)","P_value")
Survial_cancer_celltype=as.data.frame(Survial_cancer_celltype)
length(unique(Survial_cancer_celltype$Celltype))
Survial_cancer_celltype$P_value=as.numeric(Survial_cancer_celltype$P_value)
Survial_cancer_celltype$Z_score=as.numeric(Survial_cancer_celltype$Z_score)
Survial_cancer_celltype$HR=as.numeric(Survial_cancer_celltype$HR)
saveRDS(Survial_cancer_celltype,file="BLCA_Lymphocyte_Celltype_TOP50_GSVA_Coxregression_survival.rds")