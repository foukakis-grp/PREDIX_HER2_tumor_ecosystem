#############################################
#######################Fcr###################
#############################################
tpm=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/TMM-normalized-TPM.rds")
tpm=t(tpm[c("FCGR3A","FCGR3B"),])%>%as.data.frame()
FC=data.frame(patientID=rownames(tpm),FCGR3=rowMeans(tpm),
              FCGR3A=tpm$FCGR3A,
              FCGR3B=tpm$FCGR3B)
FC$patientID=substr(FC$patientID,9,12)
#############################################
###################Neoantigen################
#############################################
library(data.table);library(stringr);library(tidyverse);library(readxl)
Genotype=fread("E:/Projects/PREDIX_HER2/pVACseq/data/Optitype_HLA_PREDIX_HER2.txt")%>%as.data.frame()
table(Genotype$A1)
table(Genotype$A2)
Genotype$`HLA-A*02:01`="No"
Genotype$`HLA-A*02:01`[Genotype$A1=="A*02:01"|Genotype$A2=="A*02:01"]="Yes"
Genotype$`HLA-A*03:01`="No"
Genotype$`HLA-A*03:01`[Genotype$A1=="A*03:01"|Genotype$A2=="A*03:01"]="Yes"
Genotype$`HLA-B*07:02`="No"
Genotype$`HLA-B*07:02`[Genotype$B1=="B*07:02"|Genotype$B2=="B*07:02"]="Yes"
Genotype$`HLA-B*08:01`="No"
Genotype$`HLA-B*08:01`[Genotype$B1=="B*08:01"|Genotype$B2=="B*08:01"]="Yes"
Genotype$`HLA-B*40:01`="No"
Genotype$`HLA-B*40:01`[Genotype$B1=="B*40:01"|Genotype$B2=="B*40:01"]="Yes"
Genotype$`HLA-C*06:02`="No"
Genotype$`HLA-C*06:02`[Genotype$C1=="C*06:02"|Genotype$C2=="C*06:02"]="Yes"
Genotype$patientID=substr(Genotype$sampleID,9,12)%>%as.character()
ref=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Resource/Supertype_HAL.xls")
ref$Allele <- sub("(A\\*\\d{2})(\\d{2})", "\\1:\\2", ref$Allele)
ref$Allele <- sub("(B\\*\\d{2})(\\d{2})", "\\1:\\2", ref$Allele)
Genotype$Allele=Genotype$A1
Genotype=left_join(Genotype,ref,by="Allele")
Genotype$SupertypeA1=Genotype$Supertype;Genotype$Supertype=NULL
Genotype$Allele=Genotype$A2
Genotype=left_join(Genotype,ref,by="Allele")
Genotype$SupertypeA2=Genotype$Supertype;Genotype$Supertype=NULL
table(Genotype$SupertypeA1)
table(Genotype$SupertypeA2) # A01 A03

Genotype$Allele=Genotype$B1
Genotype=left_join(Genotype,ref,by="Allele")
Genotype$SupertypeB1=Genotype$Supertype;Genotype$Supertype=NULL
Genotype$Allele=Genotype$B2
Genotype=left_join(Genotype,ref,by="Allele")
Genotype$SupertypeB2=Genotype$Supertype;Genotype$Supertype=NULL
table(Genotype$SupertypeB1)
table(Genotype$SupertypeB2) # A01 A03

table(Genotype$SupertypeA1,Genotype$SupertypeA2)
table(Genotype$SupertypeB1,Genotype$SupertypeB2)

Genotype$A01="No"
Genotype$A01[Genotype$SupertypeA1%in%c("A01")&Genotype$SupertypeA2%in%c("A01")]="Yes"
table(Genotype$A01)
Genotype$A02="No"
Genotype$A02[Genotype$SupertypeA1=="A02"&Genotype$SupertypeA2=="A02"]="Yes"
Genotype$A03="No"
Genotype$A03[Genotype$SupertypeA1%in%c("A03")&Genotype$SupertypeA2%in%c("A03")]="Yes"
Genotype$A24="No"
Genotype$A24[Genotype$SupertypeA1%in%c("A24")&Genotype$SupertypeA2%in%c("A24")]="Yes"   

Genotype$B07="No"
Genotype$B07[Genotype$SupertypeB1=="B07"|Genotype$SupertypeB2=="B07"]="Yes"
Genotype$B08="No"
Genotype$B08[Genotype$SupertypeB1=="B08"|Genotype$SupertypeB2=="B08"]="Yes"
Genotype$B27="No"
Genotype$B27[Genotype$SupertypeB1=="B27"|Genotype$SupertypeB2=="B27"]="Yes"
Genotype$B44="No"
Genotype$B44[Genotype$SupertypeB1=="B44"|Genotype$SupertypeB2=="B44"]="Yes"
Genotype=Genotype[,c("patientID","A01","A02","A03","A24","B07","B08","B27","B44")]

#############################################
###################Neoantigen################
#############################################
DNA1=fread("E:/Projects/PREDIX_HER2/pVACseq/result/PREDIX_HER2_baseline_DNA_neoantigen_batch1.txt")%>%as.data.frame()
DNA2=fread("E:/Projects/PREDIX_HER2/pVACseq/result/PREDIX_HER2_baseline_DNA_neoantigen_batch2.txt")%>%as.data.frame()
DNA1=DNA1[DNA1$Chromosome!="Chromosome",]
DNA2=DNA2[DNA2$Chromosome!="Chromosome",]
DNA1$sampleID=str_extract(DNA1$`baseline/batch1/UE-2971-1102-0/MHC_Class_I/UE-2971-1102-0.filtered.tsv`, "UE-\\d{4}-\\d{4}-\\d{1}")
DNA2$sampleID=str_extract(DNA2$`baseline/batch2/UE-2971-1101-0/MHC_Class_I/UE-2971-1101-0.filtered.tsv`, "UE-\\d{4}-\\d{4}-\\d{1}")
colnames(DNA1)
var=c("sampleID","Chromosome","Start","Stop","Reference","Variant","Transcript","Transcript Support Level",                                              
      "Ensembl Gene ID","Variant Type","Mutation","Protein Position","Gene Name","HGVSc","HGVSp","HLA Allele",                                                            
      "Peptide Length","Sub-peptide Position","Mutation Position","MT Epitope Seq","WT Epitope Seq")
DNA=rbind(DNA1[,var],DNA2[,var])
freq=DNA%>%group_by(sampleID)%>%summarise(n = n())%>%as.data.frame()
colnames(freq)=c("sampleID","Neoantigen_DNA")
pid_DNA=fread("E:/Projects/PREDIX_HER2/pVACseq/data/pVACseq_meta_baseline.txt")
pid_DNA=pid_DNA[,c("patientID","sampleID_WES")];colnames(pid_DNA)=c("patientID","sampleID")
DNA=left_join(pid_DNA,freq,by="sampleID")
DNA$Neoantigen_DNA[is.na(DNA$Neoantigen_DNA)]=0
##
RNA=fread("E:/Projects/PREDIX_HER2/pVACseq/result/PREDIX_HER2_baseline_RNA_neoantigen.txt")%>%as.data.frame()
RNA=RNA[RNA$Chromosome!="Chromosome",]
RNA$sampleID=str_extract(RNA$`baseline_RNAfusion/UF-3016-1101-0/MHC_Class_I/UF-3016-1101-0.filtered.tsv`, "UF-\\d{4}-\\d{4}-\\d{1}")
freq=RNA%>%group_by(sampleID)%>%summarise(n = n())%>%as.data.frame()
freq$patientID=substr(freq$sampleID,9,12)
colnames(freq)=c("sampleID","Neoantigen_RNA","patientID")
pid_RNA=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.rds")%>%as.data.frame()
pid_RNA=data.frame(patientID=pid_RNA[,c("patientID")]) 
RNA=left_join(pid_RNA,freq,by="patientID")
RNA$Neoantigen_RNA[is.na(RNA$Neoantigen_RNA)]=0
RNA$sampleID=NULL
RNA$patientID=as.integer(RNA$patientID)
# merged Neoantigen
neoantigen=left_join(DNA,RNA,by="patientID")
#neoantigen$Neoantigen=neoantigen$Neoantigen_DNA+neoantigen$Neoantigen_RNA
neoantigen$Neoantigen=neoantigen$Neoantigen_DNA
table(is.na(neoantigen$Neoantigen))
neoantigen$patientID=as.character(neoantigen$patientID)
neoantigen=neoantigen[,c("patientID","Neoantigen_DNA")]
##########################################
################HLA#######################
##########################################
library(HLAdivR)
library(tidyverse)
HLA=fread("E:/Projects/PREDIX_HER2/pVACseq/data/Optitype_HLA_PREDIX_HER2.txt")
HLADiversityScore_v <- Vectorize(HLADiversityScore)

# Calculate A, B and C grantham distance
hla.df <- HLA %>%
  mutate(A.grantham = HLADiversityScore_v(A1,A2, exons23 = TRUE)) %>%
  mutate(B.grantham = HLADiversityScore_v(B1,B2, exons23 = TRUE)) %>%
  mutate(C.grantham = HLADiversityScore_v(C1,C2, exons23 = TRUE)) %>%
  mutate(meanHED = (A.grantham + B.grantham + C.grantham)/3) %>%
  mutate(any.homozygous = ifelse(A.grantham == 0, TRUE,
                                 ifelse(B.grantham == 0, TRUE,
                                        ifelse(C.grantham == 0, TRUE,
                                               FALSE))))
hla.df=hla.df[,c("sampleID","meanHED","any.homozygous")]%>%as.data.frame()
hla.df$patientID=substr(hla.df$sampleID,9,12)%>%as.character()
hla.df$sampleID=NULL
#############################################
################LOHHLA#######################
#############################################
library(data.table);library(tidyverse)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
lohhla=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/HLA/PREDIX_HER2_baseline_lohhla.csv"))
lohhla=lohhla[lohhla$message=="analysis_completed",]
heterozygosity_id=unique(lohhla$sampleID) 
lohhla$ID=paste0(lohhla$sampleID,lohhla$HLA_A_type1)
lohhla=lohhla[lohhla$PVal<0.05,]
sampleID=unique(lohhla$sampleID[(lohhla$HLA_type1copyNum_withBAFBin<0.5&lohhla$HLA_A_type1==lohhla$LossAllele)|
                                  (lohhla$HLA_type2copyNum_withBAFBin<0.5&lohhla$HLA_A_type2==lohhla$LossAllele)])
lohhla_data=data.frame(sampleID=data$sampleID,lohhla=NA)
lohhla_data$HLA_heterozygosity="No"
lohhla_data$HLA_heterozygosity[lohhla_data$sampleID%in%heterozygosity_id]="Yes"
table(lohhla_data$HLA_heterozygosity)
lohhla_data$lohhla=0
lohhla_data$lohhla[lohhla_data$sampleID%in%sampleID]=1
lohhla_data$patientID=as.character(substr(lohhla_data$sampleID,9,12))
lohhla_data$sampleID=NULL
#############################################
########Immune cell and Immune Exclusion#####
#############################################
rna=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.rds")
#Danaher <- grep("^Danaher", colnames(rna), value = TRUE)
Danaher=c("Cytolysis_score","Danaher-B-cells","Danaher-DC","Danaher-Macrophages","Danaher-T-cells","Danaher-CD8-T-cells",
          "Danaher-Neutrophils","Danaher-Cytotoxic-cells","Danaher-Treg","Danaher-NK-CD56dim-cells",
          "Danaher-Mast-cells","Danaher-NK-cells","Danaher-CD45","Danaher-TILs" )
#TIDE <- grep("^TIDE", colnames(rna), value = TRUE)
TIDE=c("TIDE_Dysfunction","TIDE_Exclusion","TIDE_CAF","TIDE_TAM_M2")
immune_cell=rna[,c("patientID",Danaher,TIDE)]
#############################################
###################TCR#######################
#############################################
dna=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")%>%as.data.frame()
dna$patientID=as.character(dna$patientID)
TCRA=fread("E:/Projects/PREDIX_HER2/Multimodal/Analyses/TcellExTRECT/TcellExTRECT_tumor_PREDIX_HER2_baseline_purity_CN_adjusted.txt")
TCRA$patientID=substr(TCRA$sample,9,12)%>%as.character()
TCRA=TCRA[,c("patientID","TCRA.tcell.fraction.adj")]
TCRA=TCRA[TCRA$patientID%in%dna$patientID,]
# RNA
rna=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.rds")%>%as.data.frame()
Tumor_Tcellclone=as.data.frame(fread('E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/TRUST4/PREDIX_HER2_Baseline_TRUST4_TCR_clonality.txt'))
colnames(Tumor_Tcellclone)=c("sampleID",'TCR_clonality')
Tumor_Bcellclone=as.data.frame(fread('E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/TRUST4/PREDIX_HER2_Baseline_TRUST4_BCR_heavy_clonality.txt'))
colnames(Tumor_Bcellclone)=c("sampleID",'BCR_clonality')
Tumor_Bcellclone[,3]=NULL
Trust4_tumor=left_join(Tumor_Tcellclone,Tumor_Bcellclone,by="sampleID")
Trust4_tumor$patientID=substr(Trust4_tumor$sampleID,9,12)
Immune_Repertoire=Trust4_tumor[Trust4_tumor$patientID%in%rna$patientID,]
Immune_Repertoire$sampleID=NULL
##############################################
###################Merge######################
##############################################
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
clin$patientID=as.character(clin$patientID)
neoantigen=as.data.frame(neoantigen)
TCRA=as.data.frame(TCRA)
data=left_join(clin,Genotype,by="patientID")%>%
  left_join(neoantigen,by="patientID")%>%
  left_join(hla.df,by="patientID")%>%  
  left_join(lohhla_data,by="patientID")%>%
  left_join(immune_cell,by="patientID")%>%
  left_join(TCRA,by="patientID")%>%
  left_join(Immune_Repertoire,by="patientID")%>%
  left_join(FC,by="patientID")
  
data$patientID=as.integer(data$patientID)
saveRDS(data,file="E:/Projects/PREDIX_HER2/Multimodal/Data/TME_subtype/pheno_immune_multiomics.rds")
