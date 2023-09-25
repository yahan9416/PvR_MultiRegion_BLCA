#Load in the Cell-cell interaction result within each sample
file_List=list.files("./Cell_Cell_interaction/CellphoneDB/out",recursive = TRUE,full.names = TRUE)
file_List=file_List[grep("significant_means",file_List)]


Interaction_gene_pair<<-NULL
#First pay attention to the interaction between immune and tumor.
batch_process_patient_communitation<-function(temp_file){
  print(temp_file)
  comm_result=read.table(temp_file,header = TRUE,sep="\t")
  sample_name=unlist(strsplit(dirname(temp_file),"\\/"))[length(unlist(strsplit(dirname(temp_file),"\\/")))]
  sample_name=gsub("_BLCA_SampleMalig_allcell","",sample_name)
  #only fource on the interaction betwen Malignant/Epi and other cells
  comm_result=comm_result[,c(2,grep(sample_name,colnames(comm_result)))]
  
  sample_list=lapply(strsplit(colnames(comm_result),"\\."),function(x){ return(x[grep(sample_name,x)]) })
  sample_list_length=unlist(lapply(sample_list, function(x) length(x)))
  Epithe_interaction=c(which(sample_list_length == 2))
  comm_result=comm_result[,-1*Epithe_interaction]
  if(length(which(duplicated(comm_result[,1]))) > 0){
    comm_result=comm_result[-1*which(duplicated(comm_result[,1])),]
  }
  rownames(comm_result)=comm_result[,1]
  comm_result=comm_result[,-1]
  comm_result=as.matrix(comm_result)
  comm_result[which(is.na(comm_result))]=0
  gene_pair=rownames(comm_result)[which(rowSums(comm_result)>0)]
  gene_pair=list(gene_pair)
  names(gene_pair)=sample_name
  Interaction_gene_pair<<-c(Interaction_gene_pair,gene_pair)
}
result=apply(matrix(file_List),1,batch_process_patient_communitation)

table(table(unlist(Interaction_gene_pair)))
sort(table(unlist(Interaction_gene_pair)),decreasing = TRUE)

##########calculate the similarity of cell communication between different sample ############

All_unique_gene_pair=unique(unlist(Interaction_gene_pair))
#generate a gene pair 
file_List=list.files("/fs/home/hanya/Project/Bladder_cancer/Cell_Cell_interaction/CellphoneDB/out",recursive = TRUE,full.names = TRUE)
file_List=file_List[grep("significant_means",file_List)]

Sampe_genepair<<-NULL
Sample_list<<-NULL
batch_generate_patient_interaction_network<-function(temp_file){
  comm_result=read.table(temp_file,header = TRUE,sep="\t")
  sample_name=unlist(strsplit(dirname(temp_file),"\\/"))[length(unlist(strsplit(dirname(temp_file),"\\/")))]
  sample_name=gsub("_BLCA_SampleMalig_allcell","",sample_name)
  comm_result=comm_result[,c(2,grep(sample_name,colnames(comm_result)))]
  
  sample_list=lapply(strsplit(colnames(comm_result),"\\."),function(x){ return(x[grep(sample_name,x)]) })
  sample_list_length=unlist(lapply(sample_list, function(x) length(x)))
  Epithe_interaction=c(which(sample_list_length == 2))
  comm_result=comm_result[,-1*Epithe_interaction]
  if(length(which(duplicated(comm_result[,1]))) > 0){
    comm_result=comm_result[-1*which(duplicated(comm_result[,1])),]
  }
  rownames(comm_result)=comm_result[,1]
  comm_result=comm_result[,-1]
  

  sample_list=na.omit(unique(unlist(sample_list)))
  comm_result=as.matrix(comm_result)
  comm_result[which(is.na(comm_result))]=0
  
  batch_generate_state_result<-function(temp_state){
    index=grep(temp_state,colnames(comm_result))
    gene_pair=rownames(comm_result)[which(rowSums(comm_result[,index])>0)]
    temp_gene_pair=rep(0,length(All_unique_gene_pair))
    temp_gene_pair[which(All_unique_gene_pair %in% gene_pair)]=1
    Sample_list<<-c(Sample_list,temp_state)
    Sampe_genepair<<-cbind(Sampe_genepair,temp_gene_pair)
  } 
  result=apply(matrix(sample_list),1,batch_generate_state_result)
}
result=apply(matrix(file_List),1,batch_generate_patient_interaction_network)

rownames(Sampe_genepair)=All_unique_gene_pair
colnames(Sampe_genepair)=Sample_list[1:10]
colSums(Sampe_genepair)

state_pair=merge(Sample_list,Sample_list)
state_pair=as.matrix(state_pair)
Sample_CCI_correlation<<-NULL
calculate_state_correlation<-function(temp_sate){
  index=match(temp_sate,colnames(Sampe_genepair))
  cor_result=unlist(cor.test(Sampe_genepair[,index[1]],Sampe_genepair[,index[2]]))[c(4,3)]
  Sample_CCI_correlation<<-rbind(Sample_CCI_correlation,c(temp_sate,cor_result))
}
result=apply(state_pair,1,calculate_state_correlation)


correlation_matrix=matrix(as.numeric(Sample_CCI_correlation[,3]),byrow=TRUE,ncol=length(Sample_list))
rownames(correlation_matrix)=Sample_list
colnames(correlation_matrix)=Sample_list

correlation_matrix=correlation_matrix[c(1,2,5,6,3,4,7:10),c(1,2,5,6,3,4,7:10)]
library(ComplexHeatmap)
library(pheatmap)
library(circlize)
setwd("./CellphoneDB/Sample_level_malignant")
pdf("MalEpi_sample_interaction_similarity.pdf")
pheatmap(correlation_matrix,cluster_rows = FALSE,cluster_cols = FALSE, color = c(rep("#414393",3),colorRampPalette(colors = c("#414393","white","#E64536"))(10),rep("#E64536",3)))
dev.off()
