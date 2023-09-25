library(MAESTRO)
library(Seurat)
library(ggpubr)
library(future)
library(Gmisc)
plan("multiprocess", workers = 12)
options(future.globals.maxSize = 10*1024^4)
SeuratObj=readRDS(",/BLCA_reBuild_Stromal_seurat.rds")
SeuratObj@meta.data$curated_anno=droplevels(SeuratObj@meta.data$curated_anno)

RNA=SeuratObj
batch=SeuratObj@meta.data$Patient
nfeatures = 3000
dims.use = 1:40
cluster.res =1
only.pos = FALSE
genes.test.use = "presto"
genes.cutoff = 1e-05
genes.pct = 0.1
genes.logfc = 0.25
runpca.agrs = list()
findneighbors.args = list()
findclusters.args = list()
runpca.agrs = list(npcs = 50)
only.pos = TRUE
RNA@meta.data$batch <- batch
data.list <- SplitObject(RNA, split.by = "batch")


data.list <- SplitObject(RNA, split.by = "batch")
for(i in 1:length(data.list)){
  data.list[[i]] <- NormalizeData(data.list[[i]], verbose = FALSE)
  data.list[[i]] <- FindVariableFeatures(data.list[[i]], selection.method = "vst", nfeatures = 3000, verbose = FALSE)
  
}
anchors <- FindIntegrationAnchors(object.list = data.list, 
                                  dims = dims.use, anchor.features = nfeatures)


RNA.integrated <- IntegrateData(anchorset = anchors, dims = dims.use)
RNA.integrated@project.name <- "BLCA_Stromal_CCA"
DefaultAssay(RNA.integrated) <- "integrated"
RNA.integrated <- ScaleData(RNA.integrated, verbose = FALSE)
RNA.integrated <- fastDoCall("RunPCA", c(object = RNA.integrated, 
                                         runpca.agrs))
p = ElbowPlot(object = RNA.integrated, ndims = RNA.integrated@commands$RunPCA.integrated@params$npcs)
ggsave(file.path(paste0(RNA.integrated@project.name, "_PCElbowPlot.png")), 
       p, width = 10, height = 4)
RNA.integrated <- RunUMAP(object = RNA.integrated, reduction = "pca", 
                          dims = dims.use)
RNA.integrated <- fastDoCall("FindNeighbors", c(object = RNA.integrated, 
                                                reduction = "pca", dims = dims.use, findneighbors.args))
RNA.integrated <- fastDoCall("FindClusters", c(object = RNA.integrated, 
                                               resolution = 1.2, findclusters.args))


p = DimPlot(object = RNA.integrated, label = TRUE, pt.size = 0.2)
ggsave(file.path(paste0(RNA.integrated@project.name, "_cluster_new.png")), 
       p, width = 5, height = 4)
p = DimPlot(object = RNA.integrated, group = "batch", label = TRUE, 
            pt.size = 0.2)
ggsave(file.path(paste0(RNA.integrated@project.name, "_batch.png")), 
       p, width = 8, height = 4)
p = DimPlot(object = RNA.integrated, group = "Sample", label = TRUE, 
            pt.size = 0.2,label.size = 2)
ggsave(file.path(paste0(RNA.integrated@project.name, "_Sample.png")), 
       p, width = 5, height = 4)

p = DimPlot(object = RNA.integrated, group = "curated_anno", label = TRUE, 
            pt.size = 0.2,label.size = 2)
ggsave(file.path(paste0(RNA.integrated@project.name, "_curated_anno.png")), 
       p, width = 6.5, height = 4)

SeuratObj=FindClusters(SeuratObj,resolution=1.1)
p = DimPlot(object = SeuratObj, label = TRUE, pt.size = 0.2)
ggsave(file.path(paste0(SeuratObj@project.name, "_cluster_res.png")), 
       p, width = 5, height = 4)

saveRDS(SeuratObj,file="BLCA_Stromal_CCA_seurat.rds")
