library(fgsea)
library(GSVA)
setwd("/Storage/hanya/BLCA")
Imputation_tpm=readRDS("BLCA_Selectgeneexp_Imputation_tpm.rds")
library(GSVA)
pathway_file="/Storage/hanya/BLCA/KEGG_Metabolism.gmt"
pathways <- gmtPathways(pathway_file)
pathway_names <- names(pathways)
pathways_gene=unique(unlist(pathways))
length(pathways_gene) #totally 1667 genes
length(intersect(rownames(Imputation_tpm),pathways_gene)) # intersect with top 10000 or all (1135genes)
pathway_activity<-gsva(Imputation_tpm,pathways)
saveRDS(pathway_activity,file="BLCA_Selectgeneexp_seurat_metabolic_gsva_score.rds")