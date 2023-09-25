library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
library(Gmisc)
plan("multiprocess", workers = 12)
options(future.globals.maxSize = 10*1024^4)
library(fgsea)
library(GSVA)
pathway_file="./KEGG_activity_score/KEGG_Metabolism.gmt"
pathways <- gmtPathways(pathway_file)
pathway_names <- names(pathways)
load("./Metabolic_result/Enzyme_relationship.RData")
Seurat_obj=readRDS("./BLCA_process_Seurat_res.rds")
#Seurat_obj=Seurat_obj$RNA
DefaultAssay(Seurat_obj)<-"RNA"
Seurat_obj <- FindVariableFeatures(object = Seurat_obj, selection.method = "vst",   nfeatures = 5000)
metabolic_gene=unique(unlist(pathways))
length(metabolic_gene)
length(intersect(metabolic_gene,VariableFeatures(Seurat_obj)))
select_gene=unique(c(VariableFeatures(Seurat_obj),metabolic_gene))
Imputation_tpm=GetAssayData(Seurat_obj)
Imputation_tpm=as.matrix(Imputation_tpm)
Imputation_tpm=Imputation_tpm[match(select_gene,rownames(Imputation_tpm)),]
#Imputation_tpm=tpm
#################################################
#Correct or enzyme
Enzyme_or_matrix=Imputation_tpm[na.omit(match(unlist(Enzyme_or_list),rownames(Imputation_tpm))),]
gene_name_or<<-NULL
Enzyme_or_matrix_correct<<-NULL
batch_process_Enzyme_or_gene_exp<-function(x){
  index=na.omit(match(x,rownames(Enzyme_or_matrix)))
  if(length(index) == 1){
    gene_name_or<<-c(gene_name_or,rownames(Enzyme_or_matrix)[index])
    Enzyme_or_matrix_correct<<-rbind(Enzyme_or_matrix_correct,Enzyme_or_matrix[index,])
  }
  if(length(index) > 1){
    temp_expr=Enzyme_or_matrix[na.omit(match(x,rownames(Enzyme_or_matrix))),]
    gene_name_or<<-c(gene_name_or,rownames(Enzyme_or_matrix)[index[1]])
    result=apply(temp_expr,2,function(y) max(y))
    Enzyme_or_matrix_correct<<-rbind(Enzyme_or_matrix_correct,as.vector(unlist(result)))
  }
}
result=lapply(Enzyme_or_list,batch_process_Enzyme_or_gene_exp)
rownames(Enzyme_or_matrix_correct)<-gene_name_or

################################################
#Correct or enzyme
Enzyme_and_matrix=Imputation_tpm[na.omit(match(unlist(Enzyme_and_list),rownames(Imputation_tpm))),]
gene_name_and<<-NULL
Enzyme_and_matrix_correct<<-NULL
batch_process_Enzyme_or_gene_exp<-function(x){
  index=na.omit(match(x,rownames(Enzyme_and_matrix)))
  if(length(index) == 1){
    gene_name_and<<-c(gene_name_and,rownames(Enzyme_and_matrix)[index])
    Enzyme_and_matrix_correct<<-rbind(Enzyme_and_matrix_correct,Enzyme_and_matrix[index,])
    #return(Enzyme_or_matrix[index,])
  }
  if(length(index) > 1){
    temp_expr=Enzyme_and_matrix[na.omit(match(x,rownames(Enzyme_and_matrix))),]
    gene_name_and<<-c(gene_name_and,rownames(Enzyme_and_matrix)[index[1]])
    result=apply(temp_expr,2,function(y) max(y))
    Enzyme_and_matrix_correct<<-rbind(Enzyme_and_matrix_correct,as.vector(unlist(result)))
    #return(result)
  }
}
result=lapply(Enzyme_or_list,batch_process_Enzyme_or_gene_exp)
rownames(Enzyme_and_matrix_correct)<-gene_name_and

replace_gene_index=unique(c(na.omit(match(unlist(Enzyme_and_list),rownames(Imputation_tpm))),na.omit(match(unlist(Enzyme_or_list),rownames(Imputation_tpm)))))

Imputation_tpm=Imputation_tpm[-1*replace_gene_index,]
Imputation_tpm=rbind(Imputation_tpm,Enzyme_and_matrix_correct)
Imputation_tpm=rbind(Imputation_tpm,Enzyme_or_matrix_correct)
saveRDS(Imputation_tpm,file="./BLCA_Selectgeneexp_Imputation_tpm.rds")
