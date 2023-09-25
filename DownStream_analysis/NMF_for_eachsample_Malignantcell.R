library(MAESTRO)
library(Seurat)
library(Biobase)
library(NMF)

Seurat_obj=readRDS("/fs/home/hanya/Project/Bladder_cancer/MalEpi_analysis/BLCA_reBuild_MaliEpi_seurat.rds")
head(Seurat_obj@meta.data)
table(Seurat_obj@meta.data$Sample,Seurat_obj@meta.data$curated_anno)
#是对Sample 做NMF呢，还是Patient 呢？
#check 一下每个sample malignant cell 的数目, 最少的Sample 有1832个细胞
MaliEpi_cells=rownames(Seurat_obj@meta.data)[which(Seurat_obj@meta.data$curated_anno  == "Malignant")]
MaliEpi_seurat=subset(Seurat_obj,cells=MaliEpi_cells)
NMF_list<<-list()
batch_for_cacnertype<-function(temp_sample){
  temp_cells=rownames(MaliEpi_seurat@meta.data)[which(MaliEpi_seurat@meta.data$Sample == temp_sample)]
  temp_seurat=subset(MaliEpi_seurat,cells=temp_cells)
  temp_seurat=CreateSeuratObject(temp_seurat@assays$RNA@counts,meta.data=temp_seurat@meta.data)
  temp_seurat=NormalizeData(temp_seurat)
  temp_seurat=FindVariableFeatures(temp_seurat)
  temp_seurat=ScaleData(temp_seurat,do.center = FALSE)
  expmat=temp_seurat@assays$RNA@scale.data
  res <- nmf(expmat, 12,  seed = 10,method="snmf/r")
  #NMF_list<<-c(NMF_list,res)
  saveRDS(res,file=paste0("/fs/home/hanya/Project/Bladder_cancer/MalEpi_analysis/NMF/",temp_sample,"MaliEpi_Sample_Prog12_NMF.rds"))
}
result=apply(matrix(unique(MaliEpi_seurat@meta.data$Sample)),1,batch_for_cacnertype)

#saveRDS(NMF_list,"/fs/home/hanya/Project/Bladder_cancer/NMF/BLCA_MaliEpi_Sample_Prog10_NMF.rds")