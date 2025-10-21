###############################################
# CUTseq data
library(readxl)
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
#Cutseq=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/all_data_by_genes.txt")%>%as.data.frame()
Cutseq=readRDS("E:/Projects/PREDIX_HER2/CUTseq/data/cna_information/PREDIX_HER2_baseline_copyratio.rds")
clinical_anno=read_excel("E:/Projects/PREDIX_HER2/Clin/data_tidy/PREDIX_HER2_23Jan2023.xlsx")%>%filter(TREAT%in%c("Experimentell","Standard"))%>%mutate(patientID=as.character(patientID))
seg_count=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/sample_seg_counts.txt")
seg_count=seg_count[order(seg_count$segment_count),]
seg_count=seg_count[!duplicated(substr(seg_count$sample,1,4)),]
seg_count=seg_count[substr(seg_count$sample,1,4)%in%clinical_anno$patientID,]
sampleid=seg_count$sample
Cutseq=Cutseq[,sampleid]
colnames(Cutseq)=substr(colnames(Cutseq),1,4)
#saveRDS(Cutseq,file="E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_CUTseq_gistic2_gene_level.rds")
#saveRDS(Cutseq,file="E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline_copyratio_gene_level.rds")
###############################################
library(data.table);library(tidyverse)
# CUTseq data
swgs=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_CUTseq_gistic2_gene_level.rds")
swgs=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline_copyratio_gene_level.rds")
# rna-seq
tpm=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/TMM-normalized-TPM.rds")
colnames(tpm)=substr(colnames(tpm),9,12)
# MS
meta=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Meta_adjust_V2.tsv')%>%as.data.frame()
MS=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Quant_adjust_V2.tsv')%>%as.data.frame()
#MS=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Quant_adjust_HarmonizR.tsv')%>%as.data.frame()
#MS=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Quant_dup_V2.tsv')%>%as.data.frame()
row.names(MS)=MS$Gene_name;MS$Gene_name=NULL
all.equal(meta$Sample_ID,colnames(MS))
colnames(MS)=meta$patientID
range(MS[,3])
####annotate######
library(biomaRt)
####intersect#####
sampleID_dna_rna=intersect(colnames(swgs),colnames(tpm))
sampleID_dna_ms=intersect(colnames(swgs),colnames(MS))
sampleID_rna_ms=intersect(colnames(tpm),colnames(MS))
gene_rna=intersect(row.names(tpm),row.names(MS))
gene=intersect(intersect(row.names(swgs),row.names(tpm)),row.names(MS))
sample=intersect(intersect(colnames(swgs),colnames(tpm)),colnames(MS))
### load function ##
library(multiOmicsViz)
compute_correlation <- function(cna_mat_filter, expr_mat_match_filter) {
  # ref the code from multiOmicsViz
  library(dplyr)
  # 计算相关矩阵
  corrArray <- cor(t(cna_mat_filter), t(expr_mat_match_filter), method = "spearman")
  corrArray[is.na(corrArray)] <- 0
  
  # 计算有效样本数
  n <- t(!is.na(t(cna_mat_filter))) %*% (!is.na(t(expr_mat_match_filter)))
  
  # 计算 t 值和 p 值
  t <- (corrArray * sqrt(n - 2)) / sqrt(1 - corrArray^2)
  corrP <- 2 * (1 - pt(abs(t), (n - 2)))
  
  # 计算 cis 相关性
  cis_cor <- data.frame(
    gene = rownames(corrP),
    r = diag(corrArray),
    p = diag(corrP)
  ) %>%
    mutate(FDR = p.adjust(p, "BH"))
  
  rownames(cis_cor) <- NULL
  cis_cor$cis_group <- "NS"
  cis_cor$cis_group[cis_cor$FDR < 0.05 & cis_cor$r > 0] <- "cis_cor"
  
  # 计算 trans 相关性
  count_fdr_below_0_05 <- apply(corrP, 1, function(row) {
    p_adjusted <- p.adjust(row, method = "BH")
    sum(p_adjusted < 0.05)
  })
  
  trans_cor <- data.frame(
    gene = rownames(corrP),
    trans_cor_freq = count_fdr_below_0_05
  )
  
  trans_cor$trans_cor_freq[trans_cor$gene %in% cis_cor$gene[cis_cor$cis_group == "cis_cor"]] <-
    trans_cor$trans_cor_freq[trans_cor$gene %in% cis_cor$gene[cis_cor$cis_group == "cis_cor"]] - 1
  trans_cor$trans_cor_freq[trans_cor$trans_cor_freq==-1]=0
  
  trans_cor$trans_group <- "NS"
  trans_cor$trans_group[trans_cor$trans_cor_freq > 50] <- "trans_cor"
  
  # 合并结果
  results <- left_join(cis_cor, trans_cor, by = "gene")
  
  return(results)
}
######CNA-mRNA######
cna=swgs[gene,sampleID_dna_rna]
rna=tpm[gene,sampleID_dna_rna]
CNA_mRNA=compute_correlation(cna,rna)
table(CNA_mRNA$cis_group)
table(CNA_mRNA$trans_group)

######CNA-prot######
cna=swgs[gene,sampleID_dna_ms]
prot=MS[gene,sampleID_dna_ms]
CNA_prot=compute_correlation(cna,prot)
table(CNA_prot$cis_group)
table(CNA_prot$trans_group)


######mRNA-prot######
rna=tpm[gene,sampleID_rna_ms]
prot=MS[gene,sampleID_rna_ms]
mRNA_prot=compute_correlation(rna,prot)
table(mRNA_prot$cis_group)
median(mRNA_prot$r,na.rm =T)
mean(mRNA_prot$r,na.rm =T)

sheet=list("CNA_mRNA_cor"=CNA_mRNA,
           "CNA_prot_cor"=CNA_prot,
           "mRNA_prot_cor"=mRNA_prot)
openxlsx::write.xlsx(sheet,file="E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/CNA_mRNA_protein_cor.xlsx")

######################
######mRNA-prot######
rna=tpm[gene_rna,sampleID_rna_ms]
prot=MS[gene_rna,sampleID_rna_ms]
mRNA_prot=compute_correlation(rna,prot)
mRNA_prot$trans_cor_freq=NULL
mRNA_prot$trans_group=NULL
mRNA_prot$cis_group=NULL
mRNA_prot$group="NS"
mRNA_prot$group[mRNA_prot$FDR<0.05]="FDR<0.05"
sheet=list("mRNA_prot_cor"=mRNA_prot)
openxlsx::write.xlsx(sheet,file="E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/mRNA_protein_cor.xlsx")




####################################
#             Run OmicsEV         #
####################################
library(data.table);library(tidyverse)
tpm=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/TMM-normalized-TPM.rds")
colnames(tpm)=paste0("P",substr(colnames(tpm),9,12)) 
tpm=round(2^tpm-1,digits = 1) 
tpm=tpm+0.1
# MS
meta=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Meta_adjust_V2.tsv')%>%as.data.frame()
meta$patientID=paste0("P",meta$patientID) 
MS=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Quant_adjust_V2.tsv')%>%as.data.frame()
MS[,2:ncol(MS)]=round(2^MS[,2:ncol(MS)],digits = 1) 
MS[,2:ncol(MS)]=MS[,2:ncol(MS)]+0.1
colnames(MS)[1]="ID"
all.equal(meta$Sample_ID,colnames(MS)[2:ncol(MS)])
colnames(MS)[2:ncol(MS)]=meta$patientID
sampleID_rna_ms=intersect(colnames(tpm),colnames(MS))
MS=MS[,c("ID",sampleID_rna_ms)]
rna=tpm[MS$ID,sampleID_rna_ms]
rna$ID=row.names(rna)
rna=rna[,c("ID",sampleID_rna_ms)]
row.names(rna)=NULL
table(rna==0)
# missing
MS_HarmonizR=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Quant_adjust_HarmonizR.tsv')%>%as.data.frame()
MS_HarmonizR[,2:ncol(MS_HarmonizR)]=round(2^MS_HarmonizR[,2:ncol(MS_HarmonizR)],digits = 1) 
MS_HarmonizR[,2:ncol(MS_HarmonizR)]=MS_HarmonizR[,2:ncol(MS_HarmonizR)]+0.1
colnames(MS_HarmonizR)[1]="ID"
all.equal(meta$Sample_ID,colnames(MS_HarmonizR)[2:ncol(MS_HarmonizR)])
colnames(MS_HarmonizR)[2:ncol(MS_HarmonizR)]=meta$patientID
MS_HarmonizR=MS_HarmonizR[MS_HarmonizR$ID%in%intersect(row.names(tpm),MS_HarmonizR$ID),c("ID",sampleID_rna_ms)]
sampleID_rna_ms=intersect(colnames(tpm),colnames(MS_HarmonizR))
MS_HarmonizR=MS_HarmonizR[,c("ID",sampleID_rna_ms)]
rna=tpm[MS_HarmonizR$ID,sampleID_rna_ms]
rna$ID=row.names(rna)
rna=rna[,c("ID",sampleID_rna_ms)]
row.names(rna)=NULL
table(rna==0)

# meta data
meta=meta[meta$patientID%in%intersect(meta$patientID,colnames(MS)),]
predic_meta=meta[,c("patientID","ER","Biol_Status")]
colnames(predic_meta)=c("sample","class","batch")
predic_meta$sample=as.character(predic_meta$sample)
predic_meta$batch[predic_meta$batch=="Original"]=1
predic_meta$batch[predic_meta$batch=="Duplicate"]=2
predic_meta$order=seq(1:nrow(predic_meta))

meta=meta[,c("patientID","PAM50_sspbc","Biol_Status")]
colnames(meta)=c("sample","class","batch")
meta$sample=as.character(meta$sample)
meta$batch[meta$batch=="Original"]=1
meta$batch[meta$batch=="Duplicate"]=2
meta$order=seq(1:nrow(meta))
# export tsv
write.table(MS,file="E:/Projects/PREDIX_HER2/Multimodal/Analyses/OmicsEV/data/protein/proteomics.tsv",quote = F,row.names =F,sep="\t")
write.table(MS_HarmonizR,file="E:/Projects/PREDIX_HER2/Multimodal/Analyses/OmicsEV/data/protein_missing/proteomics_HarmonizR.tsv",quote = F,row.names =F,sep="\t")
write.table(rna,file="E:/Projects/PREDIX_HER2/Multimodal/Analyses/OmicsEV/data/rna.tsv",quote = F,row.names =F,sep="\t")
write.table(meta,file="E:/Projects/PREDIX_HER2/Multimodal/Analyses/OmicsEV/data/sample_list_v2.tsv",quote = F,row.names =F,sep="\t")
write.table(predic_meta,file="E:/Projects/PREDIX_HER2/Multimodal/Analyses/OmicsEV/data/sample_ML.tsv",quote = F,row.names =F,sep="\t")

str(predic_meta)
