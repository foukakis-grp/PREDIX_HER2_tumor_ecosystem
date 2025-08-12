############################
#######Figure1.A############
############################
############Upset plot###############
library(ComplexUpset);library(ggvenn);library(ggsci);library("scales");library(ggplot2);library(readxl);library(data.table)
show_col(pal_npg("nrc")(10))
show_col(pal_npg("nrc", alpha = 0.6)(10))
wes=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/Mutect2_tumor_normal_baseline_PREDIX_HER2.PASS.txt")
wes=data.frame(sampleID.wes=unique(wes$Tumor_Sample_Barcode),WES=1)%>%mutate(patientID=substr(sampleID.wes,9,12))
rna=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/TMM-normalized-TPM.rds")
rna=data.frame(sampleID.rna=colnames(rna),RNAseq=1)%>%mutate(patientID=substr(sampleID.rna,9,12))
swgs=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/PREDIX_HER2_CUTseq_gistics2_baseline_amp_peak_curated.txt")
swgs=data.frame(patientID=swgs$patientID,sWGS=1)
images=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/image_metrics_PREDIX_HER2.rds")
images$patientID=as.character(images$patientID)
images$Image=1
images=images[,c("patientID","Image")]
# MS
MS=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Meta_adjust_V2.tsv')%>%as.data.frame()
MS$sampleID.ms=MS$Code
MS$MS=1
MS=MS[,c("sampleID.ms","patientID","MS")]
MS=MS[!duplicated(MS$patientID),]
# Xenium
xenium=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/xenium_meta.csv")
xenium=xenium[xenium$tpt=="pre",]
xenium$xenium=1
xenium=xenium[,c("patientID","xenium")]
xenium$patientID=as.character(xenium$patientID)
# clin
meta=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta=meta[,c("patientID","Arm","ER")]
# merge
meta$patientID=as.character(meta$patientID)
wes$patientID=as.character(wes$patientID)
swgs$patientID=as.character(swgs$patientID)
rna$patientID=as.character(rna$patientID)
MS$patientID=as.character(MS$patientID)
images$patientID=as.character(images$patientID)
data=left_join(meta,wes,by="patientID")%>%left_join(swgs,by="patientID")%>%
    left_join(rna,by="patientID")%>%left_join(MS,by="patientID")%>%
    left_join(images,by="patientID")%>%left_join(xenium,by="patientID")
data$WES[is.na(data$WES)]=0
data$sWGS[is.na(data$sWGS)]=0
data$RNAseq[is.na(data$RNAseq)]=0
data$MS[is.na(data$MS)]=0
data$Image[is.na(data$Image)]=0
data$xenium[is.na(data$xenium)]=0
table(data$WES)
table(data$sWGS)
table(data$RNAseq)
table(data$MS)
table(data$Image)
table(data$xenium)
write.csv(data, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/Fig1a.csv")
id=c("WES","sWGS","RNAseq","MS","Image","xenium")
library(ComplexHeatmap)
library(circlize)
data=data[order(data$Arm),]
ha = HeatmapAnnotation(Arm=data$Arm,
                       ER = data$ER,
                       WES = data$WES,
                       sWGS=data$sWGS,
                       RNAseq=data$RNAseq,
                       MS=data$MS,
                       Image=data$Image,
                       Xenium=data$xenium,
                       col = list(Arm =c("DHP" =  "#8491B4FF", "T-DM1" = "#91D1C2FF"),
                                  ER=c("positive" =  "#616569", "negative" = "#eeeeee"),
                                  WES=c("1" =  "#916ba6", "0" = "#eeeeee"),
                                  sWGS=c("1" =  "#916ba6", "0" = "#eeeeee"),
                                  RNAseq=c("1" =  "#916ba6", "0" = "#eeeeee"),
                                  DNA_cycle2=c("1" = "#916ba6", "0" = "#eeeeee"),
                                  MS=c("1" =  "#916ba6", "0" = "#eeeeee"),
                                  Image=c("1" =  "#916ba6", "0" = "#eeeeee"),
                                  Xenium=c("1" =  "#916ba6", "0" = "#eeeeee")))


draw(ha)

set.seed(123)
mat = matrix(rnorm(80, 2), 8, 197)
Heatmap(mat, cluster_columns=FALSE,top_annotation = ha)
# 11X6 P

#Upset
#https://krassowski.github.io/complex-upset/articles/Examples_R.html
set_size = function(w, h, factor=1.5) {
  s = 1 * factor
  options(
    repr.plot.width=w * s,
    repr.plot.height=h * s,
    repr.plot.res=100 / factor,
    jupyter.plot_mimetypes='image/png',
    jupyter.plot_scale=1
  )
}

set_size(5, 3)
upset(
  data, id,
  min_size=1,
  width_ratio=0.5,
  set_sizes=(
    upset_set_size(
      geom=geom_bar(
        aes(fill=Arm, x=group),
        width=0.5),
      position='right',
    )
  ),
  guides='over'
)

# color Experimentall: #91D1C2FF Standard: #8491B4FF; landscape 6x4




######pCR barplot#####
library(ggpubr);library(data.table)
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
table(clin$Arm)
clin$pCR[clin$Arm=="DHP"]%>%table()
clin$pCR[clin$Arm=="T-DM1"]%>%table()
# Data
df <- data.frame(Arm=c("DHP", "T-DM1"),
                 pCR=c(45.5,43.9))
ggbarplot(df, "Arm", "pCR",
          fill = "Arm", color = "white",
          palette = c("#8491B4FF","#91D1C2FF"),
          label = TRUE, lab.pos = "in", lab.col = "white")
# Portrait 3X5
median(clin$OS.time)
library("survival");library(survminer)
sfit <- survfit(Surv(EFS.time,EFS.status)~Arm, data=clin)
sfit
summary(sfit)

fit <- survfit(Surv(EFS.time,EFS.status) ~ Arm, data = clin)
ggsurvplot(fit, palette = c("#8491B4FF","#91D1C2FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))

library(tableone)
cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(Arm)+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES)+as.factor(pCR), data=clin) ##DII_density_with_supp
ShowRegTable(cox.test)
(test.ph <- cox.zph(cox.test))
############################
##########Flowchart#########
##########Data Curation#####
############################
library(data.table);library(tidyverse);library(readxl)
blood=fread("E:/Projects/PREDIX_HER2/HLAtyping/recal.bam.TSV/samplesheet_hlatyping.csv")%>%mutate(patientID=as.character(substr(sample,9,12)))
tumor=as.data.frame(fread("E:/Projects/PREDIX_HER2/HLAtyping/recal.bam.TSV/WES/samplesheet_WES.csv"))
tumor=tumor[tumor$tpt=="0",]
tumor$patientID=as.character(substr(tumor$sampleID,9,12))
cutseq=as.data.frame(fread("E:/Projects/PREDIX_HER2/CUTseq/data/meta/meta_CUTseq.txt"))%>%filter(Sample_Cohort=="PREDIX HER2",tpt=='Baseline')
CNV=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/DifferentialCNA.txt")%>%as.data.frame()
CNV$patientID=substr(CNV$sampleID,1,4);CNV$sampleID=NULL
CNV=CNV[!duplicated(CNV$patientID),!duplicated(colnames(CNV))]
clinical_anno=read_excel("E:/Projects/PREDIX_HER2/Clin/data_tidy/PREDIX_HER2_23Jan2023.xlsx")%>%filter(TREAT%in%c("Experimentell","Standard"))%>%mutate(patientID=as.character(patientID))
# for normal blood
n=intersect(clinical_anno$patientID,blood$patientID)
# for tumor
n=intersect(clinical_anno$patientID,tumor$patientID)
# for cutseq data, and final curation 
n=intersect(clinical_anno$patientID,cutseq$patientID)
n=intersect(clinical_anno$patientID,CNV$patientID)
#################################
###clin/Genomic Data curation#####
################################
#clin
library(tidyverse);library(data.table);library(readxl)
clin=read_excel("E:/Projects/PREDIX_HER2/Clin/data_tidy/PREDIX_HER2_23Jan2023.xlsx")%>%filter(TREAT%in%c("Experimentell","Standard"))%>%mutate(patientID=as.character(patientID))
clin$ER="positive"
clin$ER[clin$ERPRdic=='ERandPR-']="negative"
clin$Arm="DHP"
clin$Arm[clin$TREAT=="Experimentell"]="T-DM1"
clin$Response="pCR"
clin$Response[clin$pCR2020=="No"]="RD"
clin$Response=factor(clin$Response,levels = c("RD","pCR"))
clin$patientID=as.integer(clin$patientID)

library(tidyverse)
load("E:/Projects/PREDIX_HER2/RNAseq/data/Mats220711/PREDHER2_PD11juli.rdata")
EFS=`PREDHER2_juni22HJ-PD11juli`
rm(`PREDHER2_juni22HJ-PD11juli`)
colnames(EFS)[1]="patientID"
table(EFS$lastfustatus)
EFS=EFS[,c("patientID","PDdat","opdat")]
#read follow-up data 2024.1
surv_data=read_excel("E:/Projects/PREDIX_HER2/Clin/data/Alex240201/Pred-HER2_feb24.xlsx")
surv_data$arm=NULL
EFS=left_join(EFS,surv_data,by="patientID")

EFS$PD.status=0
EFS$PD.status[!is.na(EFS$PDdat)]=1
EFS$PD.time=as.numeric(difftime(EFS$PDdat,EFS$regdat,units="days")/30) 

table(EFS$reltyp1)
EFS$reltyp1[is.na(EFS$reltyp1)]="No_event"
EFS$LRR.status=0
EFS$LRR.status[EFS$reltyp1=="locoregional"]=1
EFS$LRR.time[EFS$reltyp1=="locoregional"]=as.numeric(difftime(EFS$relapsedt1[EFS$reltyp1=="locoregional"],EFS$regdat[EFS$reltyp1=="locoregional"],units="days")/30) 

EFS$DM.status=0
EFS$DM.status[EFS$reltyp1=="distant"]=1
EFS$DM.time[EFS$reltyp1=="distant"]=as.numeric(difftime(EFS$relapsedt1[EFS$reltyp1=="distant"],EFS$regdat[EFS$reltyp1=="distant"],units="days")/30) 

EFS$CBC.status=0
EFS$CBC.status[EFS$reltyp1=="contralateral"]=1
EFS$CBC.time[EFS$reltyp1=="contralateral"]=as.numeric(difftime(EFS$relapsedt1[EFS$reltyp1=="contralateral"],EFS$regdat[EFS$reltyp1=="contralateral"],units="days")/30) 

EFS$dead[is.na(EFS$dead)]=0
EFS$OS.status=EFS$dead
EFS$OS.time[EFS$OS.status==1]=as.numeric(difftime(EFS$deathdt[EFS$OS.status==1],EFS$regdat[EFS$OS.status==1],units="days")/30)
EFS$FU.time=as.numeric(difftime(EFS$lastdat,EFS$regdat,units="days")/30)
median(EFS$FU.time)

EFS$PD.time[is.na(EFS$PD.time)]=EFS$FU.time[is.na(EFS$PD.time)]
EFS$LRR.time[is.na(EFS$LRR.time)]=EFS$FU.time[is.na(EFS$LRR.time)]
EFS$DM.time[is.na(EFS$DM.time)]=EFS$FU.time[is.na(EFS$DM.time)]
EFS$CBC.time[is.na(EFS$CBC.time)]=EFS$FU.time[is.na(EFS$CBC.time)]
EFS$OS.time[is.na(EFS$OS.time)]=EFS$FU.time[is.na(EFS$OS.time)]

EFS$EFS.status=0
EFS$EFS.status[EFS$PD.status==1|EFS$LRR.status==1|EFS$DM.status==1|EFS$CBC.status==1|EFS$OS.status==1]=1
table(EFS$EFS.status)  
EFS = EFS %>% mutate(EFS.time = pmin(PD.time,LRR.time,DM.time,CBC.time,OS.time,FU.time))
EFS$EFS.time[EFS$patientID=="1518"]=EFS$OS.time[EFS$patientID=="1518"]

var=intersect(colnames(EFS),colnames(clin))
var=var[2:14]
clin[,var]=NULL
clin=left_join(clin,EFS,by="patientID")

#generate RFS
RFS=data.frame(patientID=clin$patientID,status=clin$reltyp1)
RFS$LRR.status=0
RFS$LRR.status[RFS$status=="locoregional"]=1
RFS$LRR.time[RFS$status=="locoregional"]=as.numeric(difftime(clin$relapsedt1[clin$reltyp1=="locoregional"],clin$opdat[clin$reltyp1=="locoregional"],units="days")/30) 

RFS$DM.status=0
RFS$DM.status[RFS$status=="distant"]=1
RFS$DM.time[RFS$status=="distant"]=as.numeric(difftime(clin$relapsedt1[clin$reltyp1=="distant"],clin$opdat[clin$reltyp1=="distant"],units="days")/30) 

RFS$OS.status=clin$dead
RFS$OS.time[RFS$OS.status==1]=as.numeric(difftime(clin$deathdt[clin$OS.status==1],clin$opdat[clin$OS.status==1],units="days")/30)
RFS$FU.time=as.numeric(difftime(clin$lastdat,clin$opdat,units="days")/30)
median(EFS$FU.time)
RFS$LRR.time[is.na(RFS$LRR.time)]=RFS$FU.time[is.na(RFS$LRR.time)]
RFS$DM.time[is.na(RFS$DM.time)]=RFS$FU.time[is.na(RFS$DM.time)]
RFS$OS.time[is.na(RFS$OS.time)]=RFS$FU.time[is.na(RFS$OS.time)]

RFS$RFS.status=0
RFS$RFS.status[RFS$LRR.status==1|RFS$DM.status==1|RFS$OS.status==1]=1
RFS = RFS %>% mutate(RFS.time = pmin(LRR.time,DM.time,OS.time,FU.time))
RFS$RFS.time[RFS$patientID=="1518"]=RFS$OS.time[RFS$patientID=="1518"]
RFS$RFS.time[RFS$patientID=="1301"]=0
table(RFS$RFS.status)

clin=left_join(clin,RFS[,c("patientID","RFS.status","RFS.time")],by="patientID")

require("survival");library("survminer")
fit <- survfit(Surv(RFS.time, RFS.status) ~ Arm, data = clin)
ggsurvplot(fit, pval = TRUE,data = clin,xlim = c(0,90))

saveRDS(clin,file="E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
write.table(clin,file="E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt",quote = F,row.names =F,sep="\t")
#wes
library(data.table)
wes=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/Mutect2_tumor_normal_baseline_PREDIX_HER2.PASS.oncokb_annotated.txt")
pid=Reduce(intersect, list(clinical_anno$patientID,tumor$patientID,blood$patientID,substr(wes$Tumor_Sample_Barcode,9,12)%>%unique())) 
sampleid=tumor$sampleID[tumor$patientID%in%pid]
wes=wes[wes$Tumor_Sample_Barcode%in%sampleid,]
purity=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/PREDIX_HER2_PureCN_baseline_purity.csv"))%>%filter(sampleID%in%sampleid)
allele.cnv=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/low_purity0.1/PREDIX_HER2_PureCN_baseline_LOH.csv"))%>%filter(sampleID%in%sampleid)
seg=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/PREDIX_HER2_PureCN_baseline_copy.csv"))%>%filter(sampleID%in%sampleid)
lohhla=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PREDIX_HER2_baseline_lohhla.csv"))
lohhla=lohhla[lohhla$sampleID%in%sampleid,]
#loh=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/PREDIX_HER2_PureCN_baseline_LOH.csv"))%>%filter(sampleID%in%sampleid)
loh_gene=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/PREDIX_HER2_PureCN_baseline_gene.csv"))%>%filter(Sampleid%in%sampleid)
colnames(loh_gene)[1]="sampleID"
# wes_gistic
#peaks
peaks=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/GISTIC2_WES/all_lesions.conf_99.txt")%>%as.data.frame()
peaks$V203=NULL
peaks=peaks[grepl("CN values",peaks$`Unique Name`),] 
peaks$Descriptor[grepl("Amplification",peaks$`Unique Name`)]=paste0("Amp",peaks$Descriptor[grepl("Amplification",peaks$`Unique Name`)])
peaks$Descriptor[grepl("Deletion",peaks$`Unique Name`)]=paste0("Del",peaks$Descriptor[grepl("Deletion",peaks$`Unique Name`)])
sampleid=colnames(peaks)[10:ncol(peaks)]
peaks=peaks[,c("Descriptor",sampleid)]
colnames(peaks)[1]="gene"
#gene-level
gene=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/GISTIC2_WES/all_data_by_genes.txt")%>%as.data.frame()
drivers=read_excel("E:/Projects/PREDIX_HER2/CUTseq/data/genelist/NIHMS68344-supplement-Supplementary_Tables_1-21/nature17676-s3/Supplementary Table 14.Driver.Events.By.Mutation.Type.01052015.v2.xlsx",sheet=4)
drivers_cna=drivers$Gene[drivers$Mutation_Type%in%c("CopyNumber")]
gene=gene[gene$`Gene Symbol`%in%drivers_cna,c("Gene Symbol",sampleid)]
colnames(gene)[1]="gene"
#merge
df=rbind(gene,peaks)
df=df[!duplicated(df$gene),]
rownames(df)=df$gene
df$gene=NULL
df=t(df)%>%as.data.frame()
df$patientID=substr(rownames(df),9,12)%>%as.character()
pid=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/PREDIX_HER2_PureCN_baseline_purity_curated.csv")
df=df[pid$sampleID,]
write.table(df,file="E:/Projects/PREDIX_HER2/Multimodal/Data/GISTIC2_WES/PREDIX_HER2_wes_gistics2_baseline_curated.txt",quote = F,row.names =F,sep="\t")
###
write.table(wes,file="E:/Projects/PREDIX_HER2/Multimodal/Data/WES/Mutect2_tumor_normal_baseline_PREDIX_HER2.PASS.oncokb_annotated_curated.txt",quote = F,row.names =F,sep="\t")
write.table(purity,file="E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/PREDIX_HER2_PureCN_baseline_purity_curated.csv",quote = F,row.names =F,sep="\t")
write.table(allele.cnv,file="E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/PREDIX_HER2_PureCN_baseline_LOH_curated.csv",quote = F,row.names =F,sep="\t")
write.table(seg,file="E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/PREDIX_HER2_PureCN_baseline_seg_curated.csv",quote = F,row.names =F,sep="\t")
write.table(lohhla,file="E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PREDIX_HER2_baseline_lohhla_curated.csv",quote = F,row.names =F,sep="\t")
write.table(loh_gene,file="E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/PREDIX_HER2_PureCN_baseline_gene_curated.csv",quote = F,row.names =F,sep="\t")
####
#cutseq_gistic
library(readxl)
clinical_anno=read_excel("E:/Projects/PREDIX_HER2/Clin/data_tidy/PREDIX_HER2_23Jan2023.xlsx")%>%filter(TREAT%in%c("Experimentell","Standard"))%>%mutate(patientID=as.character(patientID))
seg_count=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/sample_seg_counts.txt")
seg_count=seg_count[order(seg_count$segment_count),]
seg_count=seg_count[!duplicated(substr(seg_count$sample,1,4)),]
seg_count=seg_count[substr(seg_count$sample,1,4)%in%clinical_anno$patientID,]
sampleid=seg_count$sample
#peaks
peaks=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/all_lesions.conf_99.txt")%>%as.data.frame()
peaks$V203=NULL
peaks=peaks[grepl("CN values",peaks$`Unique Name`),] 
peaks=peaks[peaks$`q values`<0.1,]
peaks$Descriptor[grepl("Amplification",peaks$`Unique Name`)]=paste0("Amp",peaks$Descriptor[grepl("Amplification",peaks$`Unique Name`)])
peaks$Descriptor[grepl("Deletion",peaks$`Unique Name`)]=paste0("Del",peaks$Descriptor[grepl("Deletion",peaks$`Unique Name`)])
peaks=peaks[,c("Descriptor",sampleid)]
colnames(peaks)[1]="peak"
peaks=peaks[!duplicated(peaks$peak),]
rownames(peaks)=peaks$peak
peaks$peak=NULL
peaks=t(peaks)%>%as.data.frame()
peaks$patientID=substr(rownames(peaks),1,4)%>%as.character()
amp_peak=peaks[,grepl("Amp", colnames(peaks), ignore.case = TRUE)]%>%mutate(patientID=peaks$patientID)
del_peak=peaks[,grepl("Del", colnames(peaks), ignore.case = TRUE)]%>%mutate(patientID=peaks$patientID)

write.table(amp_peak,file="E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/PREDIX_HER2_CUTseq_gistics2_baseline_amp_peak_curated.txt",quote = F,row.names =F,sep="\t")
write.table(del_peak,file="E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/PREDIX_HER2_CUTseq_gistics2_baseline_del_peak_curated.txt",quote = F,row.names =F,sep="\t")
#gene-level
gene=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/all_data_by_genes.txt")%>%as.data.frame()
drivers=fread("E:/Projects/Collaboration/BEVPAC/CUTseq/genelist_nik-zainal-etal.tsv") 
colnames(gene)[1]="gene"
gene=gene[gene$gene%in%unique(c(drivers$Gene,"RAB11FIP1","FADD","PPFIA1","RPL19","MED1","CDK12","PPP1R1B","MIEN1","GRB7")),
          c("gene",sampleid)]
rownames(gene)=gene$gene
gene$gene=NULL
gene=t(gene)%>%as.data.frame()
gene$patientID=substr(rownames(gene),1,4)%>%as.character()
write.table(gene,file="E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/PREDIX_HER2_CUTseq_gistics2_baseline_gene_curated.txt",quote = F,row.names =F,sep="\t")
# CUTseq gistic2 and WES as complementary (threshold)
CNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/GISTIC2_WES/all_thresholded.by_genes.txt")%>%as.data.frame()
row.names(CNA)=CNA$`Gene Symbol`
CNA$`Locus ID`=NULL;CNA$Cytoband=NULL;CNA$`Gene Symbol`=NULL;
CNA=t(CNA)%>%as.data.frame()
CNA$sampleID=rownames(CNA)
row.names(CNA)=substr(row.names(CNA),9,12)
CNA$sampleID=NULL
CUTseq=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/all_thresholded.by_genes.txt")%>%as.data.frame()
row.names(CUTseq)=CUTseq$`Gene Symbol`
CUTseq$`Locus ID`=NULL;CUTseq$Cytoband=NULL;CUTseq$`Gene Symbol`=NULL;
CUTseq=t(CUTseq)%>%as.data.frame()
library(readxl)
clinical_anno=read_excel("E:/Projects/PREDIX_HER2/Clin/data_tidy/PREDIX_HER2_23Jan2023.xlsx")%>%filter(TREAT%in%c("Experimentell","Standard"))%>%mutate(patientID=as.character(patientID))
seg_count=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/sample_seg_counts.txt")
seg_count=seg_count[order(seg_count$segment_count),]
seg_count=seg_count[!duplicated(substr(seg_count$sample,1,4)),]
seg_count=seg_count[substr(seg_count$sample,1,4)%in%clinical_anno$patientID,]
sampleid=seg_count$sample
CUTseq=CUTseq[sampleid,]
row.names(CUTseq)=substr(row.names(CUTseq),1,4)
CUTseq=CUTseq[,colnames(CNA)]
all.equal(colnames(CNA),colnames(CUTseq))
pid_add=setdiff(rownames(CNA),rownames(CUTseq))
CNA=CNA[pid_add,]
CNA=rbind(CUTseq,CNA)
saveRDS(CNA,file="E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/CUTseq_gene_baseline_complemental.rds")

# CUTseq gistic2 and WES as complementary (GISTIC score)
CNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/GISTIC2_WES/all_data_by_genes.txt")%>%as.data.frame()
row.names(CNA)=CNA$`Gene Symbol`
CNA$`Locus ID`=NULL;CNA$Cytoband=NULL;CNA$`Gene Symbol`=NULL;
CNA=t(CNA)%>%as.data.frame()
CNA$sampleID=rownames(CNA)
row.names(CNA)=substr(row.names(CNA),9,12)
CNA$sampleID=NULL
CUTseq=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/all_data_by_genes.txt")%>%as.data.frame()
row.names(CUTseq)=CUTseq$`Gene Symbol`
CUTseq$`Locus ID`=NULL;CUTseq$Cytoband=NULL;CUTseq$`Gene Symbol`=NULL;
CUTseq=t(CUTseq)%>%as.data.frame()
library(readxl)
clinical_anno=read_excel("E:/Projects/PREDIX_HER2/Clin/data_tidy/PREDIX_HER2_23Jan2023.xlsx")%>%filter(TREAT%in%c("Experimentell","Standard"))%>%mutate(patientID=as.character(patientID))
seg_count=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/sample_seg_counts.txt")
seg_count=seg_count[order(seg_count$segment_count),]
seg_count=seg_count[!duplicated(substr(seg_count$sample,1,4)),]
seg_count=seg_count[substr(seg_count$sample,1,4)%in%clinical_anno$patientID,]
sampleid=seg_count$sample
CUTseq=CUTseq[sampleid,]
row.names(CUTseq)=substr(row.names(CUTseq),1,4)
CUTseq=CUTseq[,colnames(CNA)]
all.equal(colnames(CNA),colnames(CUTseq))
pid_add=setdiff(rownames(CNA),rownames(CUTseq))
CNA=CNA[pid_add,]
CNA=rbind(CUTseq,CNA)
CNA=na.omit(CNA)
saveRDS(CNA,file="E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/CUTseq_gene_baseline_GISTIC_score_complemental.rds")
#######image
meta=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
image=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Digital_image/handcraft_features_FULL.csv")%>%as.data.frame()
image$V1=NULL
image$TREAT=NULL
image$ERPRdic=NULL
image$pCR=NULL
image=image[,c("patientID","CellProps__ALL__ImmuneCells","MinDist__ALL__Tumor__ImmuneCells__Mean","median_FeD_Cent_mst__edgelength_min_div_max__ALL")]
data=left_join(image,meta,by="patientID")%>%as.data.frame()
varibale=c("CellProps__ALL__ImmuneCells","MinDist__ALL__Tumor__ImmuneCells__Mean","median_FeD_Cent_mst__edgelength_min_div_max__ALL")
data=data[,c("patientID","pCR","Arm","ER",varibale)]
colnames(data)=c("patientID","pCR","Arm","ER","Immune_Cell_prop","Distance_tumor_immune","Cell_Interaction")
varibale=c("Immune_Cell_prop","Distance_tumor_immune","Cell_Interaction")
data[,varibale]=data[,varibale] %>% mutate(across(where(is.numeric), scale))
image=data[,c("patientID","Immune_Cell_prop","Distance_tumor_immune","Cell_Interaction")]
image=image[order(image$Immune_Cell_prop,decreasing = T),]
image=image[!duplicated(image$patientID),]
saveRDS(image,file="E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/image_metrics_PREDIX_HER2.rds")



# export data for proteomics
library(data.table);library(tidyverse)
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
clin$HER2_low=0
clin$HER2_low[clin$HER2neu1=="2+"&clin$ISH_HER2_copy<6]=1 
clin$ISH_HER2_CEP17_ratio[clin$HER2neu1=="2+"&clin$ISH_HER2_copy<6]
clin$ISH_HER2_CEP17_ratio[clin$ISH_HER2_copy<6]

table(clin$HER2_low)
clin=select(clin,c("patientID","attAGE","Arm","Response","ER","HER2neu1","HER2_low","HISTGR","TUMSIZE","ANYNODES"))
purity=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/PREDIX_HER2_PureCN_baseline_purity.csv")
dna=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.txt")
dna=left_join(dna,purity,by="sampleID")
dna=select(dna,c("patientID","Purity","TMB_uniform","CNV_burden","LOH_Del_burden"))
rna=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.txt")%>%select(c("patientID","HER2DX","sspbc.subtype"))
data=left_join(clin,dna,by="patientID")%>%left_join(rna,by="patientID")
write.table(data,file="E:/Projects/PREDIX_HER2/Multimodal/Data/metadata_proteomic.txt",quote = F,row.names =F,sep="\t")

maf=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/Mutect2_tumor_normal_baseline_PREDIX_HER2.PASS.oncokb_annotated_curated.txt")
vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                 "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del",
                 "In_Frame_Ins", "Missense_Mutation")
maf$vaf=maf$t_alt_count/maf$t_depth
freq=maf%>%group_by(Tumor_Sample_Barcode)%>%summarise(n = n())
maf=maf%>%filter(Hugo_Symbol=="ERBB2",vaf>0.01,Variant_Classification%in%vc.nonSilent,n_depth>25,t_depth>25,is.na(gnomAD_AF)|gnomAD_AF<0.01)
write.table(maf,file="E:/Projects/PREDIX_HER2/Multimodal/Data/ERBB2_mut_proteomic.txt",quote = F,row.names =F,sep="\t")


# tumor size after two cycles
library(readxl);library(data.table);library(ggpubr)
data=read_excel("E:/Projects/PREDIX_HER2/Clin/data/Mats231215/Kang_PredixHer2_dec23.xlsx")
data$size_C2=NA
data$size_C2=data$radsize_C2
data$size_C2[is.na(data$size_C2)]=data$clinsize_C2
table(is.na(data$size_C2))
WES=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PREDIX_HER2_vaf_purity.txt")
cycle2=WES[substr(WES$sampleID,14,14)==2,] 
data$rna_C2="Missing"
data$rna_C2[data$patientID%in%substr(cycle2$sampleID,9,12)]="Available"
table(data$rna_C2)

ggboxplot(data, x = "rna_C2", y = "size_C2",fill = "rna_C2",
          width = 0.8,palette = c("#00AFBB", "#E7B800"))+stat_summary(fun.y=mean)

ggdensity(data, x = "size_C2",
          add = "mean", rug = TRUE,
          color = "rna_C2", palette = c("#00AFBB", "#E7B800"))

data$size_C2[data$rna_C2=="Available"]%>%na.omit()%>%median()
data$size_C2[data$rna_C2=="Missing"]%>%na.omit()%>%median()
##################


image=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Digital_image/handcraft_features.csv") 
image=image[,c("V1","patientID")]





