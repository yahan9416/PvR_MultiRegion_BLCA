meta_info=readRDS("./Metabolic/BLCA_all_lineage_subcelltype_meta_information.rds")
Metabolic_result=readRDS("./Metabolic/BLCA_Selectgeneexp_seurat_metabolic_gsva_score.rds")

pathway=c("Glycolysis / Gluconeogenesis","Citrate cycle (TCA cycle)","Pyruvate metabolism","Oxidative phosphorylation","Fatty acid degradation")
index=match(pathway,rownames(Metabolic_result))
temp_Metabolic_result=Metabolic_result[index,]
epi_mal_cells=rownames(meta_info)[which(meta_info[,1] %in% c("Malignant","Epithelial"))]
temp_Metabolic_result=temp_Metabolic_result[,match(epi_mal_cells,colnames(temp_Metabolic_result))]

patient_level=unique(meta_info[,2])
batch_calculate_Energy_metabolic_pathway_corr<-function(temp_pat){
  temp_pat_index=grep(temp_pat,colnames(temp_Metabolic_result))
  temp_patient_Metabolic_result=temp_Metabolic_result[,temp_pat_index]
  
  Correlation_result<<-NULL
  Calculate_correlation_between_malicell<-function(temp_cell_index){
    temp_exp=temp_patient_Metabolic_result[,temp_cell_index]
    result=apply(matrix(temp_cell_index:length(temp_pat_index)),1,function(pair_index){
      unlist(cor.test(temp_exp,temp_patient_Metabolic_result[,pair_index]))[4]
    })
    result=cbind(colnames(temp_patient_Metabolic_result)[temp_cell_index],cbind(colnames(expmat)[temp_cell_index:length(temp_pat_index)],result))
    Correlation_result<<-rbind(Correlation_result,result)
  }
  result=apply(matrix(1:ncol(temp_patient_Metabolic_result)),1,Calculate_correlation_between_malicell)
   saveRDS(Correlation_result,file=paste0(temp_pat,"_Malignant_cell_MetabolicEnergy_correlation.rds"))#Glutamine metabolism provides anaplerotic fuel, and restrains glucose-dependent differentia- tion and the function of macrophages and T cells
  
}
result=apply(matrix(patient_level),1,batch_calculate_Energy_metabolic_pathway_corr)