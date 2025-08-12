###############################################
#fig4a correlation between cna-rna cna-protein#
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
saveRDS(Cutseq,file="E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline_copyratio_gene_level.rds")

######CNA-mRNA######
library(multiOmicsViz)
# Define the source and target omics data
sourceOmics=swgs[gene,sample]
# Create a list of target omics data
targetOmicsList <- list()
tpm=tpm[gene,sample]
protein=MS[gene,sample]
targetOmicsList[[1]] <- tpm
targetOmicsList[[2]] <- protein

# Define the output file path
outputfile <- "E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure4/CNA-mRNA_heatmap"

# Run the multiOmicsViz function
sig <- calculateCorForTwoMatrices(matrix1 = sourceOmics,matrix2 = protein,fdr=0.05)
multiOmicsViz(sourceOmics, "CNA", "All", targetOmicsList, "mRNA", "All", 0.05, outputfile,nThreads=8)

## define function
setwd('E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure4/')
custom.function = function(x) { multiOmicsViz(sourceOmics,"CNA","All",targetOmicsList,
                                              c("RNA","protein"),"All", x, 
                                              paste("CNA correlation with protein", x, sep = "_"), 
                                              nThreads = 2, legend=TRUE)
}

ilitern_num = c(0.05)
# --------------------- parallel compute ------------------------
library(parallel)
custom.function(ilitern_num)
multiOmicsViz(sourceOmics,"CNA","All",targetOmicsList,
              "RNA","All", 0.05, "CCA_CNA correlation with RNA_0.05_scale", nThreads = NULL, legend=TRUE)


targetOmicsList <- list()
tpm=tpm[gene,]
protein=MS[gene,]
cna=swgs[gene,]
targetOmicsList[[1]] <- protein
multiOmicsViz(cna,"CNA","All",targetOmicsList,
              "protein","All", 0.05, 
              paste("CNA correlation with proteins", 0.05, sep = "_"), 
              legend=TRUE)

###############################################
#write.table(Cutseq,file="E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_CUTseq_gistic2_gene_level.txt",quote = F,row.names =F,sep="\t",)
library(data.table);library(tidyverse)
# CUTseq data
swgs=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_CUTseq_gistic2_gene_level.rds")
#swgs=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline_copyratio_gene_level.rds")
# rna-seq
tpm=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/TMM-normalized-TPM.rds")
colnames(tpm)=substr(colnames(tpm),9,12)
# MS
meta=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Meta_adjust_V2.tsv')%>%as.data.frame()
MS=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Quant_adjust_V2.tsv')%>%as.data.frame()
#MS=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Quant_dup_V2.tsv')%>%as.data.frame()
row.names(MS)=MS$Gene_name;MS$Gene_name=NULL
all.equal(meta$Sample_ID,colnames(MS))
colnames(MS)=meta$patientID
range(MS[,3])
####intersect#####
sampleID_dna_rna=intersect(colnames(swgs),colnames(tpm))
sampleID_dna_ms=intersect(colnames(swgs),colnames(MS))
sampleID_rna_ms=intersect(colnames(tpm),colnames(MS))
gene_dna_rna=intersect(row.names(swgs),row.names(tpm))
gene_dna_ms=intersect(row.names(swgs),row.names(MS))
gene_rna_ms=intersect(row.names(tpm),row.names(MS))
gene=intersect(intersect(row.names(swgs),row.names(tpm)),row.names(MS))
sample=intersect(intersect(colnames(swgs),colnames(tpm)),colnames(MS))
######CNA-mRNA######
library(multiOmicsViz)
# Define the source and target omics data
sourceOmics=swgs[gene_dna_rna,sampleID_dna_rna]
# Create a list of target omics data
targetOmicsList <- list()
targetOmicsList[[1]] <- tpm[gene_dna_rna,sampleID_dna_rna]
# Define the output file path
outputfile <- "E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS7/"
# Run the multiOmicsViz function
multiOmicsViz(sourceOmics, "CNA", "All", targetOmicsList, "mRNA", "All", 0.05, outputfile)


###############################################
#fig4a correlation between cna-rna cna-protein#
###############################################
library(tidyverse);library(data.table)
#-----------------load data---------------------#
#swgs=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline_copyratio_gene_level.rds")
swgs=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_CUTseq_gistic2_gene_level.rds")
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
####intersect#####
sampleID_dna_rna=intersect(colnames(swgs),colnames(tpm))
sampleID_dna_ms=intersect(colnames(swgs),colnames(MS))
sampleID_rna_ms=intersect(colnames(tpm),colnames(MS))
gene=intersect(intersect(row.names(swgs),row.names(tpm)),row.names(MS))
sample=intersect(intersect(colnames(swgs),colnames(tpm)),colnames(MS))
###################################################
a=calculateCorForTwoMatrices(swgs[gene,sampleID_dna_ms],MS[gene,sampleID_dna_ms],0.05)
dim(a)
table(diag(a))

#-------------------------------------------------#
# remotes::install_github("WangKang-Leo/multiOmicsViz_mod")
library(multiOmicsViz)
#------------ cna-rna ---------------#
sourceOmics=swgs[gene,sampleID_dna_rna]
# Create a list of target omics data
targetOmicsList <- list()
targetOmicsList[[1]] <- tpm[gene,sampleID_dna_rna]

# Define the output file path
outputfile <- "E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS7/CNA-mRNA"
# Run the multiOmicsViz function
multiOmicsViz(sourceOmics, "CNA", "All", targetOmicsList, "mRNA", "All", 0.05, outputfile,legend=F)
#------------ cna-protein ---------------#
sourceOmics=swgs[gene,sampleID_dna_ms]
# Create a list of target omics data
targetOmicsList <- list()
targetOmicsList[[1]] <- MS[gene,sampleID_dna_ms]
# Define the output file path
outputfile <- "E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS7/CNA-protein"
# Run the multiOmicsViz function
multiOmicsViz(sourceOmics, "CNA", "All", targetOmicsList, "protein", "All", 0.05, outputfile,legend=F)


##################################
############Figure4B##############
##################################
# function prediction
library(ggplot2)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Analyses/OmicsEV/final_res.rds")
df=data$fun_pred$data
x <- df %>% select(term,AUC,dataSet) %>% spread(key = dataSet,value = AUC)
colnames(x)[2]='Protein'
plot(x$rna,x$proteome,xlim=c(0.5,1),ylim=c(0.5,1),xlab="RNA",ylab="Protein");abline(a = 0,b=1)
d <- x %>% filter(!is.na(RNA),!is.na(Protein))
d=d[d$RNA<0.90,]
d$col <- "black"
d$col[d$Protein >1.1*d$RNA] <- "#BC3C29FF"
d$col[d$RNA >1.1*d$Protein] <- "#2D6DB1"

label=paste("AUROC(protein) > 1.1*AUROC(RNA): ",sum(d$Protein>d$RNA*1.1,na.rm = TRUE),
            "\nAUROC(RNA) > 1.1*AUROC(protein): ",sum(d$RNA>d$Protein*1.1,na.rm = TRUE),sep = "")
ggplot(d,aes(x=RNA,y=Protein)) +
  geom_point(aes(colour=col)) +
  geom_line(data=data.frame(x=c(0.5,1),y=c(0.5,1)),aes(x=x,y=y),linetype=2)+
  geom_line(data=data.frame(x=1.1*c(0.5,1),y=c(0.5,1)),aes(x=x,y=y),linetype=2)+
  geom_line(data=data.frame(x=c(0.5,1),y=1.1*c(0.5,1)),aes(x=x,y=y),linetype=2)+
  coord_cartesian(xlim =c(NA,1.01), ylim = c(NA,1.01)) +
  ggpubr::theme_classic2()+
  annotate("text",x=1,y=0.52,label=label,hjust=1,size=3.2)+
  theme(legend.position = "none")+
  scale_color_manual(breaks = c("#BC3C29FF","#2D6DB1","black"),
                     values=c("#BC3C29FF","#2D6DB1","black"))

##################################
############FigureS9g#############
##################################
library(data.table)
library(tidyverse)
library(tableone)
library(forestploter)
require(openxlsx)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
## binary variable ##
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
seg_count=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/sample_seg_counts.txt")
clinical_anno=read_excel("E:/Projects/PREDIX_HER2/Clin/data_tidy/PREDIX_HER2_23Jan2023.xlsx")%>%dplyr::filter(TREAT%in%c("Experimentell","Standard"))%>%mutate(patientID=as.character(patientID))
seg_count=seg_count[order(seg_count$segment_count),]
seg_count=seg_count[!duplicated(substr(seg_count$sample,1,4)),]
seg_count=seg_count[substr(seg_count$sample,1,4)%in%clinical_anno$patientID,]
sampleid=seg_count$sample
gene=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/all_thresholded.by_genes.txt"))
amp=gene[gene$Cytoband%in%AMP,c("Gene Symbol","Locus ID","Cytoband",sampleid)]
amp=amp[amp$`Gene Symbol`%in%c("ERLIN2","RAB11FIP1","FADD","PPFIA1","CTTN"),]
row.names(amp)=amp$`Gene Symbol`
amp$`Gene Symbol`=NULL;amp$`Locus ID`=NULL;amp$Cytoband=NULL
amp=t(amp)%>%as.data.frame()
amp$patientID=substr(row.names(amp),1,4)%>%as.integer()
data=inner_join(clin,amp,by="patientID")%>%as.data.frame()
data$CNA_ERLIN2_Amp[data$ERLIN2%in%c(0,-1,-2)]="0"
data$CNA_ERLIN2_Amp[data$ERLIN2%in%c(1,2)]="1"
data$CNA_RAB11FIP1_Amp[data$RAB11FIP1%in%c(0,-1,-2)]="0"
data$CNA_RAB11FIP1_Amp[data$RAB11FIP1%in%c(1,2)]="1"
data$CNA_FADD_Amp[data$FADD%in%c(0,-1,-2)]="0"
data$CNA_FADD_Amp[data$FADD%in%c(1,2)]="1"
data$CNA_PPFIA1_Amp[data$PPFIA1%in%c(0,-1,-2)]="0"
data$CNA_PPFIA1_Amp[data$PPFIA1%in%c(1,2)]="1"
data$CNA_CTTN_Amp[data$CTTN%in%c(0,-1,-2)]="0"
data$CNA_CTTN_Amp[data$CTTN%in%c(1,2)]="1"
bin_var=c("CNA_ERLIN2_Amp","CNA_RAB11FIP1_Amp","CNA_FADD_Amp","CNA_PPFIA1_Amp","CNA_CTTN_Amp")
data[,bin_var]=lapply(as.data.frame(data[,bin_var]),function(x) as.factor(x))
data$Arm=factor(data$Arm,levels = c("DHP","T-DM1"))
res=Logistic_batch_bin_subgroup(data,bin_var)%>%as.data.frame()

whole<- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = data[data$CNA_CTTN_Amp=="0",])
ShowRegTable(whole)
whole<- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = data[data$CNA_CTTN_Amp=="1",])
ShowRegTable(whole)
interaction_1<- glm(as.numeric(pCR) ~ CNA_CTTN_Amp+ER+Arm, family = "binomial", data = data)
interaction_2<- glm(as.numeric(pCR) ~ CNA_CTTN_Amp+ER+Arm+Arm*CNA_CTTN_Amp, family = "binomial", data = data)
lrtest(interaction_1,interaction_2)


table(data$CNA_CTTN_Amp)
##################################
############FigureS9g#############
##################################
library(tableone);library(data.table);library(tidyverse);library(ggplot2);library(ggpubr);library(sm)
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Meta_adjust_V2.tsv')%>%as.data.frame()
meta$sampleID=meta$Sample_ID
MS=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Quant_adjust_V2.tsv')%>%as.data.frame()
row.names(MS)=MS$Gene_name
MS$Gene_name=NULL
all.equal(colnames(MS),meta$Sample_ID)
colnames(MS)=meta$patientID
MS=MS[c("ERBB2","RAB11FIP1"),]%>%t()%>%as.data.frame()
range(MS$ERBB2)
range(MS$RAB11FIP1)
MS$patientID=row.names(MS)%>%as.double()
data=left_join(MS,clin,by="patientID")
data$group=paste0(data$Arm,"_",data$Response)
table(data$group)

data=data[data$Arm=="T-DM1",] 
table(data$Response[data$ERBB2>median(data$ERBB2)&data$RAB11FIP1>median(data$RAB11FIP1)])
table(data$Response[data$ERBB2>median(data$ERBB2)&data$RAB11FIP1<median(data$RAB11FIP1)])
table(data$Response[data$ERBB2<median(data$ERBB2)&data$RAB11FIP1<median(data$RAB11FIP1)])
table(data$Response[data$ERBB2<median(data$ERBB2)&data$RAB11FIP1>median(data$RAB11FIP1)])

outputDir="E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS9/"
pdf(paste0(outputDir,"FigS9g_density.pdf"),width=(5.38*2)/2.54,height=(5.1*2)/2.54, useDingbats = F, onefile = T)

par(mfrow=c(2,2), oma = c(2, 2, 0, 0),
    mar=c(2,1.5,1.5,1), font.lab=2,
    mgp = c(2.5, 1, 0), cex.axis=1.2, 
    cex.main=1.2, cex.lab=1.5)

for (q in c("DHP_RD","DHP_pCR","T-DM1_RD","T-DM1_pCR")){
  ERBB2   <- data[data$group %in% q,"ERBB2"]
  RAB11FIP1 <- data[data$group %in% q,"RAB11FIP1"]
  y     <- cbind(`ERBB2 protein`=ERBB2,`RAB11FIP1 protein` = RAB11FIP1)
  sm.density(y,display="image",panel=F,ylim=c(-1,2),xlim=c(-3,3),xlab="",ylab="")
  title(q, line = 0.7)
  abline(v=mean(data$ERBB2),lty="dashed")
  abline(h=mean(data$RAB11FIP1),lty="dashed")
}
mtext("RAB11FIP1 (protein)",side=2,line=1,outer=TRUE,cex=1,las=0,font=1)
mtext("ERBB2 (protein)",side=1,line=1,outer=TRUE,cex=1,font=1)
dev.off()


library(tableone);library(data.table);library(tidyverse);library(ggplot2);library(ggpubr);library(sm)
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Meta_adjust_V2.tsv')%>%as.data.frame()
meta$sampleID=meta$Sample_ID
MS=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Quant_adjust_V2.tsv')%>%as.data.frame()
row.names(MS)=MS$Gene_name
MS$Gene_name=NULL
all.equal(colnames(MS),meta$Sample_ID)
colnames(MS)=meta$patientID
MS=MS[c("ERBB2","ERLIN2"),]%>%t()%>%as.data.frame()
range(MS$ERBB2)
range(MS$ERLIN2)
MS$patientID=row.names(MS)%>%as.double()
data=left_join(MS,clin,by="patientID")
data$group=paste0(data$Arm,"_",data$Response)
table(data$group)

data=data[data$Arm=="T-DM1",] 
table(data$Response[data$ERBB2>median(data$ERBB2)&data$ERLIN2>median(data$ERLIN2)])
table(data$Response[data$ERBB2>median(data$ERBB2)&data$ERLIN2<median(data$ERLIN2)])
table(data$Response[data$ERBB2<median(data$ERBB2)&data$ERLIN2<median(data$ERLIN2)])
table(data$Response[data$ERBB2<median(data$ERBB2)&data$ERLIN2>median(data$ERLIN2)])

outputDir="E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS8/"
pdf(paste0(outputDir,"FigS9g_density_ERLIN2.pdf"),width=(5.38*2)/2.54,height=(5.1*2)/2.54, useDingbats = F, onefile = T)

par(mfrow=c(2,2), oma = c(2, 2, 0, 0),
    mar=c(2,1.5,1.5,1), font.lab=2,
    mgp = c(2.5, 1, 0), cex.axis=1.2, 
    cex.main=1.2, cex.lab=1.5)

for (q in c("DHP_RD","DHP_pCR","T-DM1_RD","T-DM1_pCR")){
  ERBB2   <- data[data$group %in% q,"ERBB2"]
  ERLIN2 <- data[data$group %in% q,"ERLIN2"]
  y     <- cbind(`ERBB2 protein`=ERBB2,`ERLIN2 protein` = ERLIN2)
  sm.density(y,display="image",panel=F,ylim=c(-1,2),xlim=c(-3,3),xlab="",ylab="")
  title(q, line = 0.7)
  abline(v=mean(data$ERBB2),lty="dashed")
  abline(h=mean(data$ERLIN2),lty="dashed")
}
mtext("ERLIN2 (protein)",side=2,line=1,outer=TRUE,cex=1,las=0,font=1)
mtext("ERBB2 (protein)",side=1,line=1,outer=TRUE,cex=1,font=1)
dev.off()
##################################
############FigureS9g#############
##################################
library(tableone);library(data.table);library(tidyverse);library(ggplot2);library(ggpubr)
MS=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/ERBB2_MS.rds")
MS=MS[,c("ERBB2","patientID","HER2_low")]
colnames(MS)[1]="HER2_prot"
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
clin$Response[clin$ISH_HER2_copy<6&clin$Arm=="T-DM1"]
# CNA CUTseq
seg_count=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/sample_seg_counts.txt")
seg_count=seg_count[order(seg_count$segment_count),]
seg_count=seg_count[!duplicated(substr(seg_count$sample,1,4)),]
seg_count=seg_count[substr(seg_count$sample,1,4)%in%clin$patientID,]
sampleid=seg_count$sample
CNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/all_data_by_genes.txt")%>%as.data.frame()
amp=c("ERBB2")
CNA=CNA[CNA$`Gene Symbol`%in%amp,c("Gene Symbol",sampleid)]
row.names(CNA)=CNA$`Gene Symbol`
CNA$`Gene Symbol`=NULL
CNA=as.data.frame(t(CNA))
CNA$patientID=substr(row.names(CNA),1,4)%>%as.integer() 
colnames(CNA)[1]="ERBB2_CN"
# CNA WES
WES_CNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/GISTIC2_WES/all_data_by_genes.txt")%>%as.data.frame()
amp=c("ERBB2")
WES_CNA=WES_CNA[WES_CNA$`Gene Symbol`%in%amp,]
row.names(WES_CNA)=WES_CNA$`Gene Symbol`
WES_CNA$`Gene Symbol`=NULL;WES_CNA$`Gene ID`=NULL;WES_CNA$Cytoband=NULL
WES_CNA=as.data.frame(t(WES_CNA))
WES_CNA$patientID=substr(row.names(WES_CNA),9,12)%>%as.integer() 
colnames(WES_CNA)[1]="ERBB2_CN"
pid=setdiff(WES_CNA$patientID,CNA$patientID)
WES_CNA=WES_CNA[WES_CNA$patientID%in%pid,]
# merge CNA
CNA=rbind(CNA,WES_CNA)
# mRNA
tpm=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/TMM-normalized-TPM.rds")
tpm=tpm["ERBB2",]
tpm=t(tpm)%>%as.data.frame()
tpm$patientID=substr(row.names(tpm),9,12)%>%as.double()
colnames(tpm)[1]="ERBB2_mRNA"
# data
data=left_join(clin,CNA,by="patientID")%>%left_join(tpm,by="patientID")%>%left_join(MS,by="patientID")
data$CNA=scale(data$ERBB2_CN)
data$mRNA=scale(data$ERBB2_mRNA)
data$Protien=scale(data$HER2_prot)
data$group="HER2-positive"
data$group[data$ISH_HER2_copy<6]="HER2-low"
#########Plot ############
data=data[,c("group","CNA","mRNA","Protien")]
data=gather(data,key = "omic",value = "value",-group)
P=ggplot(data, aes(x=omic, y=value,fill=group)) + 
  geom_violin(trim=FALSE,color="white") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线图
  scale_fill_manual(values = c("#38618c","#ff5964"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  theme(panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank(),
        legend.position = "top")+
  ylab("ERBB2")+xlab("")+ ylim(-3,3)
P+stat_compare_means(aes(group =group), label = "p.format")
# 4.5X4
##################################
############FigureS9h#############
##################################
library(tableone);library(data.table);library(tidyverse);library(ggplot2);library(ggpubr)
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Meta_adjust_V2.tsv')%>%as.data.frame()
meta$sampleID=meta$Sample_ID
MS=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Quant_adjust_V2.tsv')%>%as.data.frame()
row.names(MS)=MS$Gene_name
MS$Gene_name=NULL
all.equal(colnames(MS),meta$Sample_ID)
colnames(MS)=meta$patientID
MS=MS[c("TOP2A"),]%>%t()%>%as.data.frame()
colnames(MS)[1]="TOP2A_prot"
MS$patientID=row.names(MS)%>%as.double()
meta=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/ERBB2_PG.rds")
# CNA CUTseq
seg_count=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/sample_seg_counts.txt")
seg_count=seg_count[order(seg_count$segment_count),]
seg_count=seg_count[!duplicated(substr(seg_count$sample,1,4)),]
seg_count=seg_count[substr(seg_count$sample,1,4)%in%clin$patientID,]
sampleid=seg_count$sample
CNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/all_data_by_genes.txt")%>%as.data.frame()
amp=c("TOP2A")
CNA=CNA[CNA$`Gene Symbol`%in%amp,c("Gene Symbol",sampleid)]
row.names(CNA)=CNA$`Gene Symbol`
CNA$`Gene Symbol`=NULL
CNA=as.data.frame(t(CNA))
CNA$patientID=substr(row.names(CNA),1,4)%>%as.integer() 
colnames(CNA)[1]="TOP2A_CN"
# CNA WES
WES_CNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/GISTIC2_WES/all_data_by_genes.txt")%>%as.data.frame()
amp=c("TOP2A")
WES_CNA=WES_CNA[WES_CNA$`Gene Symbol`%in%amp,]
row.names(WES_CNA)=WES_CNA$`Gene Symbol`
WES_CNA$`Gene Symbol`=NULL;WES_CNA$`Gene ID`=NULL;WES_CNA$Cytoband=NULL
WES_CNA=as.data.frame(t(WES_CNA))
WES_CNA$patientID=substr(row.names(WES_CNA),9,12)%>%as.integer() 
colnames(WES_CNA)[1]="TOP2A_CN"
pid=setdiff(WES_CNA$patientID,CNA$patientID)
WES_CNA=WES_CNA[WES_CNA$patientID%in%pid,]
# merge CNA
CNA=rbind(CNA,WES_CNA)
# mRNA
tpm=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/TMM-normalized-TPM.rds")
tpm=tpm["TOP2A",]
tpm=t(tpm)%>%as.data.frame()
tpm$patientID=substr(row.names(tpm),9,12)%>%as.double()
colnames(tpm)[1]="TOP2A_mRNA"
# data
data=left_join(meta,CNA,by="patientID")%>%left_join(tpm,by="patientID")%>%left_join(MS,by="patientID")
data$CNA=scale(data$TOP2A_CN)
data$mRNA=scale(data$TOP2A_mRNA)
data$Protien=scale(data$TOP2A_prot)
data$group="HER2-positive"
data$group[data$HER2_prot%in%c("Pseudo_HER2")]="HER2-low"
#########Plot ############
data=data[,c("group","CNA","mRNA","Protien")]
data=gather(data,key = "omic",value = "value",-group)
P=ggplot(data, aes(x=omic, y=value,fill=group)) + 
  geom_violin(trim=FALSE,color="white") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  geom_boxplot(width=0.2,position=position_dodge(0.9))+ #绘制箱线图
  scale_fill_manual(values = c("#38618c","#ff5964"))+ #设置填充的颜色
  theme_bw()+ #背景变为白色
  theme(panel.grid.major = element_blank(),   #不显示网格线
        panel.grid.minor = element_blank(),
        legend.position = "top")+
  ylab("ERBB2")+xlab("")+ ylim(-3,3)
P+stat_compare_means(aes(group =group), label = "p.format")


##################################
##################################
############FigureS9h#############
##################################
library(tidyverse);library(tableone);library(data.table);library(lmtest);library(forestploter)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/ERBB2_PG.rds")
table(data$ERBB2_PG,data$Arm)
whole<- glm(as.numeric(pCR) ~ ERBB2_PG+ER, family = "binomial", data = data)
ShowRegTable(whole)
whole<- glm(as.numeric(pCR) ~ ERBB2_PG+ER, family = "binomial", data = data[data$Arm=="DHP",])
ShowRegTable(whole)
whole<- glm(as.numeric(pCR) ~ ERBB2_PG+ER, family = "binomial", data = data[data$Arm=="T-DM1",])
ShowRegTable(whole)

library(forestploter)
require(openxlsx)
# forest plot #
library(data.table)
library(readxl)
df=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS9/Forest.xlsx")
df$HR=ifelse(is.na(df$OR), "", df$OR)
df$Low=ifelse(is.na(df$Low), "", df$Low)
df$High=ifelse(is.na(df$High), "", df$High)
df$`OR (95% CI)`=ifelse(is.na(df$`OR (95% CI)`), "", df$`OR (95% CI)`)
df$P=ifelse(is.na(df$P), "", df$P)
df$` ` <- paste(rep(" ", 10), collapse = " ")
a=df$`OR (95% CI)`;b=df$P
df$`OR (95% CI)`=NULL;df$`P`=NULL
df$`OR (95% CI)`=a
df$`P`=b
df$OR=as.numeric(df$OR)
df$Low=as.numeric(df$Low)
df$High=as.numeric(df$High)


tm <- forest_theme(base_size = 10,
                   refline_col = "black",
                   arrow_type = "closed",
                   footnote_col = "blue")

p <- forest(df[,c(1,6:8)],
            est = df$OR,
            lower = df$Low, 
            upper = df$High,
            ci_column = 2,
            ref_line = 1,
            arrow_lab = c("RD", "pCR"),
            xlim = c(0, 10),
            ticks_at = c(0, 1, 2, 4,6,8),
            theme = tm)

# Print plot
plot(p)
# 4X2 P 65%
##################################
############FigureS9i#############
##################################
library(data.table);library(ggplot2);library(irr);library(tidyverse)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/ERBB2_PG.rds")
data$ERBB2_PG_pos="No"
data$ERBB2_PG_pos[data$ERBB2_PG=="Positive"]="Yes"
data$Her2E="No"
data$Her2E[data$sspbc.subtype=="Her2"]="Yes"
table(data$Her2E,data$ERBB2_PG_pos)
kappa_result <- kappa2(data[,c("Her2E","ERBB2_PG_pos")], weight = "unweighted")
print(kappa_result)
# Cohen’s Kappa value (Altman 1999, Landis JR (1977)).  Kappa = 0.319 fair agreement
df=table(data$ERBB2_PG_pos,data$Her2E)%>%as.data.frame()
colnames(df)=c("ERBB2_PG_pos","Her2E","Freq")
f=ggplot(data =df, mapping = aes(x =ERBB2_PG_pos, y = Her2E)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
  scale_fill_gradient(low="white", high="#009194") +
  theme_bw()
f
#4X3

library(data.table);library(ggplot2);library(irr);library(tidyverse)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/ERBB2_PG.rds")
data$ERBB2_PG_pos="No"
data$ERBB2_PG_pos[data$ERBB2_PG=="Positive"]="Yes"
data$Her2_prolif="No"
data$Her2_prolif[data$Her2ADC=="S1"]="Yes"
table(data$Her2_prolif,data$ERBB2_PG_pos)
kappa_result <- kappa2(data[,c("Her2_prolif","ERBB2_PG_pos")], weight = "unweighted")
print(kappa_result)
# Cohen’s Kappa value (Altman 1999, Landis JR (1977)).  Kappa = 0.319 fair agreement
df=table(data$ERBB2_PG_pos,data$Her2_prolif)%>%as.data.frame()
colnames(df)=c("ERBB2_PG_pos","Her2_prolif","Freq")
f=ggplot(data =df, mapping = aes(x =ERBB2_PG_pos, y = Her2_prolif)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
  scale_fill_gradient(low="white", high="#009194") +
  theme_bw()
f
##################################
############FigureS9j#############
##################################
library(ggplot2);library(ggpubr);library(dplyr)
mydata=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/ERBB2_PG.rds")
chisq.test(mydata$Response,mydata$ERBB2_class) #0.006
table(mydata$EFS.status,mydata$ERBB2_class)
mydata=mydata[!is.na(mydata$ERBB2_class),]
mydata$ERBB2_class=factor(mydata$ERBB2_class,levels = c("Her2E ERBB2 PG+",
      "Her2E ERBB2 PG-","Other ERBB2 PG+","Other ERBB2 PG-"))
mytable=mydata %>%
  group_by(ERBB2_class,Response) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
ggbarplot(mytable, "ERBB2_class", "freq",
          fill = "Response", color = "Response", 
          palette = c("#fdb462","#00A087FF"))+coord_flip()

#6X4

library(ggplot2);library(ggpubr);library(dplyr)
mydata=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/ERBB2_PG.rds")
mydata$ERBB2_class=NA
mydata$ERBB2_class[mydata$ERBB2_PG=="Negative"&mydata$Her2ADC=="S1"]="Her2E ERBB2 PG-"
mydata$ERBB2_class[mydata$ERBB2_PG=="Negative"&mydata$Her2ADC!="S1"]="Other ERBB2 PG-"
mydata$ERBB2_class[mydata$ERBB2_PG=="Positive"&mydata$Her2ADC=="S1"]="Her2E ERBB2 PG+"
mydata$ERBB2_class[mydata$ERBB2_PG=="Positive"&mydata$Her2ADC!="S1"]="Other ERBB2 PG+"

chisq.test(mydata$Response,mydata$ERBB2_class) #0.006
table(mydata$EFS.status,mydata$ERBB2_class)
mydata=mydata[!is.na(mydata$ERBB2_class),]
mydata$ERBB2_class=factor(mydata$ERBB2_class,levels = c("Her2E ERBB2 PG+",
                                                        "Her2E ERBB2 PG-","Other ERBB2 PG+","Other ERBB2 PG-"))
mytable=mydata %>%
  group_by(ERBB2_class,Response) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))
ggbarplot(mytable, "ERBB2_class", "freq",
          fill = "Response", color = "Response", 
          palette = c("#fdb462","#00A087FF"))+coord_flip()


