library(copykat)
library(MAESTRO)
library(Seurat)
library(ggplot2)


Merge_h5=Read10X_h5("/fs/home/hanya/Project/Bladder_cancer/Raw_data/After_process/HHZ_B/HHZ_B_filtered_feature_bc_matrix/HHZ_B_gene_count.h5")
Merge_h5=as.matrix(Merge_h5)
setwd("/fs/home/hanya/Project/Bladder_cancer/Basic_analysis/CopyKat")
copykat.test <- copykat(rawmat=Merge_h5, id.type="S", ngene.chr=1, win.size=25, KS.cut=0.1, sam.name="HHZ_B", distance="euclidean", n.cores=20)
saveRDS(copykat.test,file="HHZ_B_copykat_clustering_results")