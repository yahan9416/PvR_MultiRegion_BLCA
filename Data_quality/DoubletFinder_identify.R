library(DoubletFinder)
#library(MAESTRO)
library(Seurat)
library(ggplot2)
library(future)
files=list.files("/fs/home/hanya/Project/Bladder_cancer/After_process",recursive=T)
files=files[grep(".h5",files)]


Batch_identify_Doublet_for_each_sample<-function(x){
  expr = Read10X_h5(paste0("/fs/home/hanya/Project/Bladder_cancer/After_process/",x))
  seu_kidney <- CreateSeuratObject(expr)
  seu_kidney <- NormalizeData(seu_kidney)
  seu_kidney <- FindVariableFeatures(seu_kidney, selection.method = "vst", nfeatures = 2000)
  seu_kidney <- ScaleData(seu_kidney)
  seu_kidney <- RunPCA(seu_kidney)
  seu_kidney <- RunUMAP(seu_kidney, dims = 1:10)
  ## pK Identification (no ground-truth)
  sweep.res.list_kidney <- paramSweep_v3(seu_kidney, PCs = 1:10, sct = FALSE)
  sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
  bcmvn_kidney <- find.pK(sweep.stats_kidney)
  
  ## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
  seu_kidney <- FindNeighbors(seu_kidney, dims = 1:10)
  seu_kidney <- FindClusters(seu_kidney, resolution = 0.5)
  annotations <- seu_kidney@meta.data$seurat_clusters 
  homotypic.prop <- modelHomotypic(annotations)           ## ex: 
  nExp_poi <- round(0.075*dim(seu_kidney@meta.data)[1])  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  
  ## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
  seu_kidney <- doubletFinder_v3(seu_kidney, PCs = 1:10, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  return(cbind(rownames(seu_kidney@meta.data),seu_kidney@meta.data[,7]))
}
result=apply(matrix(files),1,Batch_identify_Doublet_for_each_sample)
sample=lapply(strsplit(files,"/"),function(x) x[1])
names(result)=unlist(sample)
saveRDS(result,file="/fs/home/hanya/Project/Bladder_cancer/Doublets/DoubletFinder/DoubletFinder_result.rds")
