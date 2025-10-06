library(TcellExTRECT);library(foreach);library(data.table);library(tidyverse)
setwd("/proj/sens2022005/nobackup/TcellExTRECT/result_purity_CN_adjusted")
bed.file <-"/castor/project/proj_nobackup/nf-core/nf-core-sarek-2.7.1/Sarek-data/Twist_Comprehensive_Exome_Covered_Targets_GRCh38.pad.sort.merge.bed"
exons_hg38_custom <- createExonDFBed(bed.file,'hg38')
colnames(exons_hg38_custom)=c("chr","start","stop")
exons_hg38_custom=as.tibble(exons_hg38_custom)
meta=as.data.frame(fread("/proj/sens2022005/nobackup/TcellExTRECT/data/samplesheet_WES_baseline_tumor.csv"))
sampleID=meta$sampleID
data(tcra_seg_hg38)
Tcell= foreach (i=sampleID,.combine=rbind) %do% {
  cat(i, "...\n")
  bam=meta[meta$sampleID==i,"bam"]
  Purity=meta[meta$sampleID==i,"Purity"]
  CN=meta[meta$sampleID==i,"TCRA_CN"]
  cov.file <- getCovFromBam(bamPath = bam,
                          outPath ='',
                          vdj.seg = tcra_seg_hg38)
  cov_df <- loadCov(cov.file)                       
  TCRA.out <- as.data.frame(runTcellExTRECT(vdj.region.df=cov_df,exons.selected=exons_hg38_custom,vdj.seg=tcra_seg_hg38,'hg38',GC_correct = TRUE,sample_name=i)) 
  TCRA.out <- as.data.frame(adjustTcellExTRECT(TCRA.out, purity =Purity,TCRA.cn =CN))
}
write.table(Tcell,file="/proj/sens2022005/nobackup/TcellExTRECT/result_purity_CN_adjusted/TcellExTRECT_tumor_PREDIX_HER2_baseline_purity_CN_adjusted.txt",quote = F,row.names =T,sep = "\t")
