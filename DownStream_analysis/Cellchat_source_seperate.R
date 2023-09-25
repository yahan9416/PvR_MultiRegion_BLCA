library(Seurat) # please update to Seurat V4
library(tidyverse)
library(MAESTRO)
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 10*1024^4)
Seurat_obj=readRDS("BLCA_process_Seurat_res.rds")
setwd("/.Cell_Cell_interaction")
meta.data=readRDS("./Cell_Cell_interaction/BLCA_meta_data.rds")
DefaultAssay(Seurat_obj)<-"RNA"
expmat=GetAssayData(Seurat_obj)




cells_primary=rownames(meta.data)[which(meta.data$Source == "Primary")]
Seurat_primary=subset(Seurat_obj,cells=cells_primary)



primary_exp=GetAssayData(Seurat_primary,slot="data")
meta.data$curated_anno=droplevels(meta.data$curated_anno)
meta_infp=meta.data[match(colnames(primary_exp),rownames(meta_info)),]
cellchat <- createCellChat(object = primary_exp, meta = meta_infp, group.by = "curated_anno")

groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database


future::plan("multiprocess", workers = 10)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 1)
cellchat <- computeCommunProbPathway(cellchat)


cellchat <- aggregateNet(cellchat)
saveRDS(cellchat,file="BLCA_primary_cellchat.rds")
