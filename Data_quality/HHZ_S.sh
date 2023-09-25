celescope rna sample --outdir ./00.sample --sample HHZ_S --thread 15 --chemistry auto  --fq1 ./HHZ-230301-S_hbwoz2vovo_R1.fastq.gz 
celescope rna barcode --outdir ./01.barcode --sample HHZ_S --thread 15 --chemistry auto --lowNum 2  --fq1 ./HHZ-230301-S_hbwoz2vovo_R1.fastq.gz --fq2 ./HHZ-230301-S_hbwoz2vovo_R2.fastq.gz
celescope rna cutadapt --outdir ./02.cutadapt --sample HHZ_S --thread 15 --minimum_length 20 --nextseq_trim 20 --overlap 10 --insert 150  --fq ./01.barcode/HHZ_S_2.fq
celescope rna star --outdir ./03.star --sample HHZ_S --thread 15 --genomeDir /fs/home/hanya/CeleScope/hs_ensembl_99 --outFilterMultimapNmax 1 --starMem 30  --fq ./02.cutadapt/HHZ_S_clean_2.fq 
celescope rna featureCounts --outdir ./04.featureCounts --sample HHZ_S --thread 15 --gtf_type gene --genomeDir /fs/home/hanya/CeleScope/hs_ensembl_99  --input ./03.star/HHZ_S_Aligned.sortedByCoord.out.bam 
celescope rna count --outdir ./05.count --sample HHZ_S --thread 15 --genomeDir /fs/home/hanya/CeleScope/hs_ensembl_99 --expected_cell_num 3000 --cell_calling_method EmptyDrops_CR  --bam ./04.featureCounts/HHZ_S_name_sorted.bam --force_cell_num None 
celescope rna analysis --outdir ./06.analysis --sample ZXW_HB --thread 15 --genomeDir /fs/home/hanya/CeleScope/hs_ensembl_99  --matrix_file ./05.count/HHZ_S_filtered_feature_bc_matrix 
