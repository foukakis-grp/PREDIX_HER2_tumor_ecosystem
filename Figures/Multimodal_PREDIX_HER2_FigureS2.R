#############################################
##########Figure2A,B TMB ~ purity,ER#########
#############################################
library(ggpubr);library(ggplot2);library(data.table);library(tidyverse)
genomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
purity=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/PREDIX_HER2_PureCN_baseline_purity_curated.csv"))
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic$patientID=as.integer(genomic$patientID)
data=left_join(genomic,purity,by="sampleID")%>%left_join(clin,by="patientID")
data=data[data$totalTMB*41.2>5,]
ggscatter(data, x = "Purity", y = "totalTMB",
                add = "reg.line",  # Add regressin line
                xlim=c(0.1,0.75),ylim=c(0,20),
                add.params = list(color = "#0E4C92", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE # Add confidence interval
) + stat_cor(method = "pearson", label.x = 0.3, label.y = 20) 

#6X6 50%
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source(paste0(baseDir,"/Code/theme.R"))
d=data%>%select(c("ER","totalTMB"))
d <- reshape2::melt(d,id.vars=c("ER"))
figure_font_size=13
Fig <-
  ggplot(d,aes(x=ER,y=value,fill=ER))+
  geom_boxplot(outlier.size = 0.5, width=0.7)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_ER(name="ER status")+
  stat_compare_means(aes(group=ER),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.30)+
  labs(y="TMB(mut/MB)",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))+ 
  coord_cartesian(ylim = c(0,20))
Fig

#4X5 p 
#############################################
######FigureS2c ERBB2/TP53_lollipop########
#############################################
library(data.table);library(tidyverse)
maf=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/Mutect2_tumor_normal_baseline_PREDIX_HER2.PASS.oncokb_annotated_vaf0.05_curated.rds")
maf$vaf=maf$t_alt_count/maf$t_depth
vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                 "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del",
                 "In_Frame_Ins", "Missense_Mutation")
maf=maf%>%filter(Variant_Classification%in%vc.nonSilent)
table(maf$ONCOGENIC)

a=maf[maf$Hugo_Symbol=="TP53"&maf$VARIANT_IN_ONCOKB==T,]
a=maf[maf$Hugo_Symbol=="ERBB2"&maf$VARIANT_IN_ONCOKB==T&maf$ONCOGENIC%in%c("Oncogenic","Likely Oncogenic"),]
table(maf$VARIANT_IN_ONCOKB,maf$ONCOGENIC)

maf=select(maf,c("Tumor_Sample_Barcode","Hugo_Symbol","Chromosome","Start_Position","End_Position","Variant_Classification",'Reference_Allele',"Tumor_Seq_Allele2","Amino_acids"))
colnames(maf)=c("Sample_ID","Hugo_Symbol","Chromosome","Start_Position","End_Position","Mutation_Type",'Reference_Allele',"Variant_Allele","Protein_Change")
df=maf[maf$Hugo_Symbol%in%c("ERBB2"),]

write.table(df,file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Appeal/FigureS2/ERBB2_TP53.txt",quote = F,row.names =F,sep="\t")

#HER2 mutation "UE-2971-1119-0" "UE-2971-1121-0" "UE-2971-1208-0" "UE-2971-1215-0" "UE-2971-1514-0" "UE-2971-1516-0" "UE-2971-1520-0"
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
genomic=left_join(genomic,clin,by="patientID")

genomic$pCR2020[genomic$sampleID%in%c("UE-2971-1514-0","UE-2971-1516-0","UE-2971-1119-0",
                                      "UE-2971-1121-0","UE-2971-1208-0")]
table(genomic$coding_mutation_ERBB2)
#############################################
############FigureS2d forest plot############
#############################################
library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
library(forestplot)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
genomic=left_join(genomic,clin,by="patientID")
variable=c("coding_mutation_TP53","coding_mutation_PIK3CA","coding_mutation_ERBB2",
           "coding_mutation_TP53_oncokb","coding_mutation_PIK3CA_oncokb")
results=Logistic_batch_adjER(genomic,"pCR","Arm",variable,"ER")%>%as.data.frame()
colnames(results)
Total=results[,c("biomarker","whole_OR","Whole_LCI","Whole_UCI","whole_lr_p")]
Total$group="All"
DHP=results[,c("biomarker","DHP_OR","DHP_LCI","DHP_UCI","DHP_lr_p")]
DHP$group="DHP"
TDM1=results[,c("biomarker","TDM1_OR","TDM1_LCI","TDM1_UCI","TDM1_lr_p")]
TDM1$group="T-DM1"
colname=c("biomarker","OR","LCI","UCI","p","group")
colnames(Total)=colname;colnames(DHP)=colname;colnames(TDM1)=colname
df=rbind(Total,DHP,TDM1)
df$biomarker=gsub("coding_mutation_", "",df$biomarker)
#df$biomarker=gsub("_oncokb", "",df$biomarker)
biomarker=df$biomarker
df$OR=as.numeric(df$OR);df$LCI=as.numeric(df$LCI);df$UCI=as.numeric(df$UCI)
#
df |>
  group_by(group) |>
  forestplot(labeltext=biomarker,
             mean=OR,lower=LCI,upper=UCI,
             zero = 1,
             boxsize = .25, # We set the box size to better visualize the type
             line.margin = .1, # We need to add this to avoid crowding
             lty.ci = c(3),
             clip = c(0,4),xlab = " OR with 95% CI") |> 
  fp_add_lines("black") |> 
  fp_add_header("Mutation") |> 
  fp_set_style(box = c("#e5c06e", "#8491B4FF","#91D1C2FF") |> lapply(function(x) gpar(fill = x, col = "black")),
               default = gpar(vertices = TRUE)) |> 
  fp_set_zebra_style("#F5F9F9")


#TransNeo
meta=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/transneo/41586_2021_4278_MOESM4_ESM.xlsx",sheet=1)%>%as.data.frame()
meta=meta[meta$HER2.status=="POS"&meta$pCR.RD!="NA",]
maf=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/transneo/41586_2021_4278_MOESM4_ESM.xlsx",sheet=2)%>%as.data.frame()
maf$Tumor_Sample_Barcode=maf$Donor.ID;maf$Chromosome=maf$Chr;maf$Start_Position=maf$Start;maf$End_Position=maf$Start
maf$Reference_Allele=maf$Ref_Allele;maf$Tumor_Seq_Allele2=maf$Tumor_Alt_Allele
maf=maf[maf$Donor.ID%in%meta$Donor.ID,] 
write.table(maf,file="E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/transneo/TransNeo_HER2_maf.txt",
            quote = F,sep = "\t",row.names = F)

maf=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/transneo/TransNeo_HER2_maf.oncokb_annotated.txt")
maf=maf[maf$VARIANT_IN_ONCOKB==T&maf$ONCOGENIC%in%c("Likely Oncogenic","Oncogenic")&maf$Hugo_Symbol=="TP53",]
meta=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/transneo/41586_2021_4278_MOESM4_ESM.xlsx",sheet=1)%>%as.data.frame()
meta=meta[meta$HER2.status=="POS"&meta$pCR.RD!="NA",]
meta$TP53="No"
meta$TP53[meta$Donor.ID%in%maf$Donor.ID]="Yes"
table(meta$TP53,meta$pCR.RD)
library(tableone)
meta$pCR=1
meta$pCR[meta$pCR.RD=="RD"]=0
res<- glm(pCR ~ TP53+ER.status, family = "binomial", data =meta)
ShowRegTable(res)
library(vcd);library("ggsci")
meta$pCR.RD=factor(meta$pCR.RD,levels = c("RD","pCR"))
meta$TP53=factor(meta$TP53,levels = c("Yes","No"))
mosaic( ~ ER.status + pCR.RD + TP53, data = meta,
        highlighting = "TP53", highlighting_fill = c("#BC3C29FF","#6F99ADFF"),
        direction = c("h","v","h"))


#############################################
############FigureS2e/f   KM plot############
#############################################
#DHP KM 
library("survival");library(survminer);library(tableone)
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
genomic=left_join(genomic,clin,by="patientID")
table(genomic$coding_mutation_TP53_oncokb)
dhp=genomic[genomic$Arm=="DHP",]
dhp$TP53="Wild type"
dhp$TP53[dhp$coding_mutation_TP53_oncokb==1]="Mut"
dhp$TP53=factor(dhp$TP53,levels = c("Wild type","Mut"))
fit <- survfit(Surv(EFS.time,EFS.status) ~ TP53, data = dhp)
ggsurvplot(fit, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))

table(dhp$coding_mutation_TP53_oncokb)
cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(TP53)+as.factor(Response)+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=dhp) ##DII_density_with_supp
(test.ph <- cox.zph(cox.test))
ShowRegTable(cox.test)

#T-DM1 KM 
tdm1=genomic[genomic$Arm=="T-DM1",]
tdm1$TP53="Wild type"
tdm1$TP53[tdm1$coding_mutation_TP53_oncokb==1]="Mut"
tdm1$TP53=factor(tdm1$TP53,levels = c("Wild type","Mut"))
fit <- survfit(Surv(EFS.time,EFS.status) ~ TP53, data = tdm1)
ggsurvplot(fit, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))


table(tdm1$coding_mutation_TP53_oncokb)
cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(TP53)+as.factor(Response)+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=tdm1) ##DII_density_with_supp
(test.ph <- cox.zph(cox.test))
ShowRegTable(cox.test)

# 6X4.5
########################################
#######TP53 mutation validation#########
# https://doi.org/10.1002/cam4.4652
########################################
#MSKCC -DHP
library("survival");library(survminer);library(tableone);library(data.table);library(tidyverse)
tumor=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/brca_mapk_hp_msk_2021/brca_mapk_hp_msk_2021/data_clinical_sample.txt")
#tumor=tumor[tumor$SAMPLE_TYPE=="Primary",]
meta=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/brca_mapk_hp_msk_2021/brca_mapk_hp_msk_2021/data_clinical_patient.txt")
meta=left_join(meta,tumor,by="PATIENT_ID")
mut=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/brca_mapk_hp_msk_2021/brca_mapk_hp_msk_2021/data_mutations.oncokb_annotated.txt")
mut$patientID=substr(mut$Tumor_Sample_Barcode,1,9)
meta=meta[meta$PATIENT_ID%in%mut$patientID,]
mut=mut[mut$Hugo_Symbol=="TP53",]
mut=mut[mut$VARIANT_IN_ONCOKB==TRUE,]
mut=mut[mut$ONCOGENIC%in%c("Likely Oncogenic","Oncogenic")]
mut$vaf=mut$t_alt_count/(mut$t_alt_count+mut$t_ref_count)
range(mut$vaf)
mut=mut[mut$vaf>0.05,]
length(unique(mut$patientID))
meta$TP53="Wild type"
meta$TP53[meta$PATIENT_ID%in%mut$patientID]="Mut"
table(meta$TP53)

meta$TP53=factor(meta$TP53,levels = c("Wild type","Mut"))
meta$PFS.STATUS=substr(meta$PFS_STATUS,1,1)
meta$PFS.STATUS=as.numeric(meta$PFS.STATUS)
meta=meta%>%filter(!is.na(meta$PFS.STATUS)) 
#meta=meta[meta$DX_PSTAGE%in%c("IV"),]

fit <- survfit(Surv(PFS_MONTHS,PFS.STATUS) ~ TP53, data = meta)
ggsurvplot(fit, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))

table(meta$TP53)
cox.test <- coxph(Surv(PFS_MONTHS,PFS.STATUS)~as.factor(TP53)+as.factor(DX_PSTAGE), data=meta) ##DII_density_with_supp
(test.ph <- cox.zph(cox.test))
ShowRegTable(cox.test)

# TCGA 
#tumor=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/brca_tcga_pan_can_atlas_2018/brca_tcga_pan_can_atlas_2018/data_clinical_sample.txt")
meta=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/brca_tcga_pan_can_atlas_2018/brca_tcga_pan_can_atlas_2018/data_clinical_patient_HER2.txt")
meta$ER_STATUS_BY_IHC
meta=meta[meta$HER2_FISH_STATUS=="Positive"|meta$IHC_HER2=="Positive",]
mut=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/brca_tcga_pan_can_atlas_2018/brca_tcga_pan_can_atlas_2018/data_mutations.oncokb_annotated.txt")
mut=mut%>%filter(VARIANT_IN_ONCOKB==TRUE,
                 ONCOGENIC%in%c("Likely Oncogenic","Oncogenic"),Hugo_Symbol=="TP53")
mut$patientID=substr(mut$Tumor_Sample_Barcode,1,12)
meta$TP53="Wild type"
meta$TP53[meta$PATIENT_ID%in%mut$patientID]="Mut"
meta$TP53=factor(meta$TP53,levels = c("Wild type","Mut"))
table(meta$TP53)
meta$DFS.STATUS=substr(meta$DFS_STATUS,1,1)%>%as.numeric()
meta$DFS_MONTHS=as.numeric(meta$DFS_MONTHS)
meta=meta[!is.na(meta$DFS_MONTHS),]
trt=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/brca_tcga_pan_can_atlas_2018/brca_tcga_pan_can_atlas_2018/data_timeline_treatment.txt")
meta=meta[meta$PATIENT_ID%in%trt$PATIENT_ID[trt$AGENT%in%c("Trastuzumab","Paclitaxel + Doxorubicin + Cyclophosphamide + Trastuzumab",
                                                           "Docetaxel + Carboplatin + Trastuzumab")],]
meta$stage[meta$AJCC_PATHOLOGIC_TUMOR_STAGE%in%c("Stage I","Stage IA")]="I"
meta$stage[meta$AJCC_PATHOLOGIC_TUMOR_STAGE%in%c("Stage IIA","Stage IIB")]="II"
meta$stage[meta$AJCC_PATHOLOGIC_TUMOR_STAGE%in%c("Stage IIIA","Stage IIIC")]="III"

fit <- survfit(Surv(DFS_MONTHS,DFS.STATUS) ~ TP53, data = meta)
ggsurvplot(fit, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
cox.test <- coxph(Surv(DFS_MONTHS,DFS.STATUS)~as.factor(TP53)+as.factor(ER_STATUS_BY_IHC)+stage, data=meta) ##DII_density_with_supp
(test.ph <- cox.zph(cox.test))
ShowRegTable(cox.test)


# SUMMIT trial
library("survival");library(survminer);library(tableone);library(data.table);library(tidyverse)
meta=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/summit_2018/summit_2018/data_clinical_patient.txt")
meta=meta[meta$COHORT=="Breast",]
mut=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/summit_2018/summit_2018/data_mutations.oncokb_annotated.txt")
mut$patientID=mut$Tumor_Sample_Barcode
#meta=meta[meta$PATIENT_ID%in%mut$patientID,]
mut=mut[mut$VARIANT_IN_ONCOKB==TRUE,]
mut=mut[mut$ONCOGENIC%in%c("Likely Oncogenic","Oncogenic")]
mut=mut[mut$Hugo_Symbol=="TP53",]
meta$TP53="Wild type"
meta$TP53[meta$PATIENT_ID%in%mut$patientID]="Mut"
table(meta$TP53)

meta$TP53=factor(meta$TP53,levels = c("Wild type","Mut"))
meta$PFS.STATUS=0
meta$PFS.STATUS[meta$BEST_ORR%in%c("PD","SD")]=1

fit <- survfit(Surv(PFS_MONTHS,PFS.STATUS) ~ TP53, data = meta)
ggsurvplot(fit, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 5,
           xlim = c(0, 20))

cox.test <- coxph(Surv(PFS_MONTHS,PFS.STATUS)~as.factor(TP53), data=meta) ##DII_density_with_supp
(test.ph <- cox.zph(cox.test))
ShowRegTable(cox.test)



# METABRIC
tumor=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/brca_metabric/brca_metabric/data_clinical_sample.txt")
tumor=tumor[tumor$HER2_STATUS=="Positive",]
meta=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/brca_metabric/brca_metabric/data_clinical_patient.txt")
meta=left_join(tumor,meta,by="PATIENT_ID")
mut=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/brca_metabric/brca_metabric/data_mutations.oncokb_annotated.txt")
mut=mut[mut$Tumor_Sample_Barcode%in%tumor$SAMPLE_ID,]
mut=mut[mut$VARIANT_IN_ONCOKB==TRUE,]
mut=mut[mut$ONCOGENIC%in%c("Likely Oncogenic","Oncogenic")]
mut=mut[mut$Hugo_Symbol=="TP53",]
meta$TP53="Wild type"
meta$TP53[meta$PATIENT_ID%in%mut$Tumor_Sample_Barcode]="Mut"
meta$TP53=factor(meta$TP53,levels = c("Wild type","Mut"))
table(meta$TP53)
meta$HR="Negative"
meta$HR[meta$ER_STATUS=="Positive"|meta$PR_STATUS=="Positive"]="Positive"
meta=meta%>%filter(!is.na(meta$TUMOR_STAGE),TMB_NONSYNONYMOUS>0) # 

meta$RFS.STATUS=substr(meta$RFS_STATUS,1,1)%>%as.numeric()
meta=meta[meta$RFS_MONTHS>12,]


fit <- survfit(Surv(RFS_MONTHS,RFS.STATUS) ~ TP53, data = meta)
ggsurvplot(fit, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
cox.test <- coxph(Surv(RFS_MONTHS,RFS.STATUS)~as.factor(TP53)+as.factor(HR)+TUMOR_STAGE+
                    HORMONE_THERAPY+RADIO_THERAPY+CHEMOTHERAPY, data=meta) ##DII_density_with_supp
(test.ph <- cox.zph(cox.test))
ShowRegTable(cox.test)


#MSKCC T-DM1
tumor=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/msk_chord_2024/msk_chord_2024/data_clinical_sample.txt")
tumor=tumor[tumor$SAMPLE_TYPE=="Primary"&tumor$PRIMARY_SITE=="Breast",]
meta=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/msk_chord_2024/msk_chord_2024/data_clinical_patient.txt")
meta=meta[meta$PATIENT_ID%in%tumor$PATIENT_ID,]
meta=meta[meta$GENDER=="Female"&meta$HER2=="Yes"&meta$STAGE_HIGHEST_RECORDED%in%c("Stage 1-3","Stage 4"),]
tumor=tumor[tumor$PATIENT_ID%in%meta$PATIENT_ID,]
meta=left_join(meta,tumor,by="PATIENT_ID")
trt=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/msk_chord_2024/msk_chord_2024/data_timeline_treatment.txt")
trt=trt[trt$PATIENT_ID%in%meta$PATIENT_ID,]
trt$trt_duration=trt$STOP_DATE-trt$START_DATE
trt=trt[trt$trt_duration>100,]
dhp=intersect(trt$PATIENT_ID[trt$AGENT=="PERTUZUMAB"],trt$PATIENT_ID[trt$AGENT=="TRASTUZUMAB"])%>%unique()
dhp=c(dhp,unique(trt$PATIENT_ID[trt$AGENT=="HYALURONIDASE/PERTUZUMAB/TRASTUZUMAB"]))

meta=meta[meta$PATIENT_ID%in%dhp,]
tumor=tumor[tumor$PATIENT_ID%in%dhp,]
mut=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/msk_chord_2024/msk_chord_2024/data_mutations.oncokb_annotated.txt")
mut=mut[mut$VARIANT_IN_ONCOKB==TRUE&mut$ONCOGENIC%in%c("Likely Oncogenic","Oncogenic"),]
#mut=mut[mut$ONCOGENIC%in%c("Likely Oncogenic","Oncogenic")]
mut=mut[mut$Tumor_Sample_Barcode%in%tumor$SAMPLE_ID&mut$Hugo_Symbol=="TP53",]
tumor=tumor[tumor$SAMPLE_ID%in%mut$Tumor_Sample_Barcode,]

meta$TP53="Wild type"
meta$TP53[meta$PATIENT_ID%in%tumor$PATIENT_ID]="Mut" 
table(meta$TP53)

meta$TP53=factor(meta$TP53,levels = c("Wild type","Mut"))
meta$OS.STATUS=substr(meta$OS_STATUS,1,1)%>%as.numeric()

fit <- survfit(Surv(OS_MONTHS,OS.STATUS) ~ TP53, data = meta)
ggsurvplot(fit, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
table(meta$TP53)
cox.test <- coxph(Surv(OS_MONTHS,OS.STATUS)~as.factor(TP53)+as.factor(HR)+CURRENT_AGE_DEID+CLINICAL_GROUP, data=meta) ##DII_density_with_supp
(test.ph <- cox.zph(cox.test))
ShowRegTable(cox.test)
#############################################
##########Figure2C  volcano plot#############
#############################################
library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
library(forestplot);library(ggpubr)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")%>%as.data.frame()
bin_var=c("coding_mutation_ERBB2","coding_mutation_PIK3CA","coding_mutation_TP53","coding_mutation_ABCA13",          
"coding_mutation_AHNAK","coding_mutation_ANK2","coding_mutation_APOB","coding_mutation_BIRC6",
"coding_mutation_CACNA1E","coding_mutation_CCDC168","coding_mutation_CDK12","coding_mutation_CELSR3",
"coding_mutation_DNAH10","coding_mutation_DNAH14","coding_mutation_DNAH6","coding_mutation_FCGBP",
"coding_mutation_FLG","coding_mutation_FSIP2","coding_mutation_IGFN1","coding_mutation_KMT2C",
"coding_mutation_MED1","coding_mutation_MED24","coding_mutation_MUC16","coding_mutation_MUC17",
"coding_mutation_MUC5AC","coding_mutation_MYO15B","coding_mutation_OBSCN","coding_mutation_PKHD1L1",
"coding_mutation_RNF213","coding_mutation_RYR1","coding_mutation_SCN4A","coding_mutation_SHROOM2",
"coding_mutation_STARD9","coding_mutation_SYNE1","coding_mutation_TAF1L","coding_mutation_TTN",
"coding_mutation_USH2A","coding_mutation_ZFHX4","coding_mutation_ZNF469")

bin_var=c("coding_mutation_HER_pathway","coding_mutation_PIK3_AKT_pathway",
          "coding_mutation_MAPK_ERK_pathway","coding_mutation_CDK_RB_pathway")
genomic$patientID=as.integer(genomic$patientID)
genomic[,bin_var]=lapply(as.data.frame(genomic[,bin_var]),function(x) as.factor(x))
genomic=left_join(genomic,clin,by="patientID")%>%as.data.frame()
str(genomic)
colnames(genomic)
results=Logistic_batch_adjER(genomic,"pCR","Arm",bin_var,"ER")%>%as.data.frame()
results$whole_OR=as.numeric(results$whole_OR)
results$whole_lr_p=as.numeric(results$whole_lr_p)
results$lnOR=results$whole_OR%>%log(base=exp(1))
results$log10p=-log10(results$whole_lr_p)
results$biomarker=gsub("coding_mutation_","",results$biomarker)
library(EnhancedVolcano)
plotVolcano=function(res,title){
  EnhancedVolcano(res,
                  lab = res$biomarker,
                  x = 'lnOR',
                  y = 'whole_lr_p',
                  xlab="lnOR",
                  xlim = c(-2.5,2),
                  ylim = c(0,max(res$log10p)),
                  pCutoff = 0.05,
                  FCcutoff = 1,
                  pointSize = 2.5,
                  labSize = 3.0,
                  labCol = 'black',
                  labFace = 'bold',
                  boxedLabels = TRUE,
                  colAlpha = 4/5,
                  legendPosition = 'top',
                  legendLabSize = 12,
                  legendIconSize = 3.0,
                  drawConnectors = T,
                  widthConnectors = 0.01,
                  colConnectors="black",
                  title=title)}

plotVolcano(results,"Non-synonymous Mutation (N=39)")
# 5X6 
#############################################
##########Figure2D  K-M plot    #############
#############################################
library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
library(forestplot);library(ggpubr)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")%>%as.data.frame()
bin_var=c("coding_mutation_ERBB2","coding_mutation_PIK3CA","coding_mutation_TP53",
          "coding_mutation_ERBB2_oncokb","coding_mutation_PIK3CA_oncokb","coding_mutation_TP53_oncokb")

genomic[,bin_var]=lapply(as.data.frame(genomic[,bin_var]),function(x) as.factor(x))
genomic$patientID=as.integer(genomic$patientID)
genomic=left_join(genomic,clin,by="patientID")%>%as.data.frame()
library(Blasso);library(tidyverse)
Mut=genomic[,c("sampleID",bin_var)]
colnames(Mut)[1]="ID"
clin=genomic[,c("sampleID","EFS.status",'EFS.time')] 
colnames(clin)=c("ID","status","time")
res<-best_predictor_cox(target_data = clin, 
                        features = Mut, 
                        status = "status",
                        time = "time",
                        nfolds = 5,
                        permutation = 1000)
res
library("survival");library(survminer)
genomic$TP53="Mut"
genomic$TP53[genomic$coding_mutation_TP53_oncokb=="0"]="Wild-type"
table(genomic$TP53)
fit <- survfit(Surv(EFS.time,EFS.status) ~ TP53, data = genomic[genomic$Arm=="DHP",])
ggsurvplot(fit, palette = c("#BC3C29FF","#6F99ADFF"),
           pval = TRUE,
           ggtheme = theme_bw(), # Change ggplot2 theme
           break.time.by = 10,
           xlim = c(0, 60))
table(genomic$IGFN1)
# 4X4
library(tableone)
genomic$TP53=factor(genomic$TP53,levels = c("Wild-type","Mut"))
table(genomic$TP53)
cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(TP53)+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=genomic[genomic$Arm=="T-DM1",]) ##DII_density_with_supp
ShowRegTable(cox.test)

cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(TP53)+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=genomic[genomic$Arm=="DHP",]) ##DII_density_with_supp
ShowRegTable(cox.test)

table(genomic$coding_mutation_ERBB2_oncokb,genomic$Response,genomic$Arm)
library(tableone)
cox.test <- coxph(Surv(EFS.time,EFS.status)~TMB_uniform+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES)+Arm, data=genomic) ##DII_density_with_supp
ShowRegTable(cox.test)

##################################
#############Figure2F#############
##################################
library(vcd);library("ggsci");library(tidyverse);library(data.table)
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.txt")
genomic=left_join(genomic,clin,by="patientID")
genomic$Response=factor(genomic$Response,levels = c("RD","pCR"))
genomic$TP53[genomic$coding_mutation_TP53_oncokb==1]="Mutation"
genomic$TP53[genomic$coding_mutation_TP53_oncokb==0]="Wild type"


ftable(genomic$Arm,genomic$ER,genomic$TP53,genomic$Response)
ftable(genomic$Arm)
mosaic( ~ Arm + ER + Response + TP53, data = genomic,
        highlighting = "TP53", highlighting_fill = c("#BC3C29FF","#6F99ADFF"),
        direction = c("v","h","v","h"))

# 5X5 55%
DHP=genomic[genomic$Arm=="DHP",]
TDM1=genomic[genomic$Arm=="T-DM1",]

table(DHP$coding_mutation_TP53_oncokb,DHP$Response)
table(TDM1$coding_mutation_TP53_oncokb,TDM1$Response)

interaction_2<- glm(pCR ~ coding_mutation_TP53_oncokb*Arm+coding_mutation_TP53_oncokb+Arm+ER, family = "binomial", data =genomic)
library(tableone)
ShowRegTable(interaction_2)


interaction_2<- glm(pCR ~ coding_mutation_TP53_oncokb+ER, family = "binomial", data =genomic[genomic$Arm=="T-DM1",])
library(tableone)
ShowRegTable(interaction_2)

##################################
#############FigureS2C############
##################################
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
genomic=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.txt")
#genomic$patientID=as.character(genomic$patientID)
genomic=left_join(genomic,clin,by="patientID")
genomic%>%group_by(Arm,Response)%>%summarise(mean=median(TMB_uniform))
d=genomic%>%select(c("Arm","Response","totalTMB","TMB_uniform","TMB_clone","TMB_oncogenic"))
d <- reshape2::melt(d,id.vars=c("Arm","Response"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source(paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
Fig <-
  ggplot(d,aes(x=Arm,y=value,fill=Response))+
  geom_boxplot(outlier.size = 0.5, width=0.7)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Treatment Response")+
  stat_compare_means(aes(group=Response),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.98)+
  labs(y="mut/MB",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))
Fig


# 10X5
##################################
#############FigureS2D############
##################################
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
genomic=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.txt")
#genomic$patientID=as.character(genomic$patientID)
genomic=left_join(genomic,clin,by="patientID")
d=genomic%>%select(c("Arm","Response","ER","subclone_per"))
d <- reshape2::melt(d,id.vars=c("Arm","Response","ER"))
d$value=100*d$value
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
Fig <-
  ggplot(d,aes(x=Arm,y=value,fill=Response))+
  geom_boxplot(outlier.size = 0.5, width=0.7)+
  facet_wrap(variable~ER,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Treatment Response")+
  stat_compare_means(aes(group=Response),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.98)+
  labs(y="Subclone Percentage (%)",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))
Fig

# 5X5
##################################
#############FigureS2e############
##################################
library(deconstructSigs);library(BSgenome.Hsapiens.UCSC.hg38);library(data.table);library(tidyverse);library(circlize);library(ComplexHeatmap)
# https://www.nature.com/articles/s41586-019-1056-z#Sec2
mutect=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/Mutect2_tumor_normal_baseline_PREDIX_HER2.PASS.oncokb_annotated_curated.txt")
maf=mutect
maf$vaf=maf$t_alt_count/maf$t_depth
#maf=maf%>%filter(Variant_Type=="SNP",vaf>=0.01,n_depth>25,t_depth>25,is.na(gnomAD_AF)|gnomAD_AF<0.01)
#freq=maf%>%group_by(Tumor_Sample_Barcode)%>%summarise(n = n())
#tumor_barcode=freq$Tumor_Sample_Barcode[freq$n>10]
#maf=maf%>%filter(Tumor_Sample_Barcode%in%tumor_barcode)
# Convert to deconstructSigs input
sigs.input <- mut.to.sigs.input(mut.ref = maf, 
                                sample.id = "Tumor_Sample_Barcode", 
                                chr = "Chromosome", 
                                pos = "Start_Position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele2",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)


w=lapply(row.names(sigs.input), function(i){
  ## signatures.cosmic signatures.nature2013
  sample_1 = whichSignatures(tumor.ref = sigs.input[,], 
                             signatures.ref = signatures.cosmic, 
                             sample.id =  i, 
                             contexts.needed = TRUE,
                             tri.counts.method = 'exome')
  print(i)
  return(sample_1$weights)
})
w=do.call(rbind,w)
w$patientID=substr(row.names(w),9,12)%>%as.integer()
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
sig=left_join(w,clin,by="patientID")
sig$mut_sig=rowSums(sig[,c("Signature.2","Signature.3","Signature.6","Signature.7","Signature.10","Signature.13","Signature.15")])
sig=sig[order(-sig$mut_sig),]
given_order=sig$patientID%>%as.integer()

long_data <- sig[,1:31]%>%gather(key = "signature", value = "value", -patientID)
table(long_data$patientID)
long_data=left_join(long_data,clin,by='patientID')
long_data$patientID=paste0("PID",long_data$patientID)
long_data=long_data%>%filter(signature%in%c("Signature.2","Signature.3","Signature.6","Signature.7","Signature.10","Signature.13","Signature.15"))
long_data=long_data%>%filter(value!=0)
# 计算每个 patientID 的总和
sum_values <- aggregate(value ~ patientID, data = long_data, sum)

# 根据总和值对 patientID 进行排序
sorted_patientIDs <- sum_values[order(sum_values$value, decreasing = TRUE), ]$patientID

# 将 patientID 转换为 factor，并按照排序后的顺序重新赋值
long_data$patientID <- factor(long_data$patientID, levels = sorted_patientIDs)

library("scales");library(ggsci)
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
ggplot(long_data, aes(fill = signature, y = value, x = patientID)) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_x_discrete(limits = sorted_patientIDs)+scale_fill_brewer(palette="Set1")+
  theme(panel.background = element_rect(fill='transparent'),
        plot.background = element_rect(fill='transparent', color=NA),
        strip.background = element_blank(),
        axis.text.x=element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))

# 6X5 70%
##################################
#############FigureS2f############
##################################
library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
library(forestplot)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.txt")
variable=c("COSMIC.Signature.2","COSMIC.Signature.3","COSMIC.Signature.6","COSMIC.Signature.7","COSMIC.Signature.10",
           "COSMIC.Signature.13","COSMIC.Signature.15")
genomic=left_join(genomic,clin,by="patientID")%>%as.data.frame()
colnames(genomic)
results=Logistic_batch_adjER(genomic,"pCR","Arm",variable,"ER")%>%as.data.frame()
colnames(results)
Total=results[,c("biomarker","whole_OR","Whole_LCI","Whole_UCI","whole_lr_p")]
Total$group="All"
DHP=results[,c("biomarker","DHP_OR","DHP_LCI","DHP_UCI","DHP_lr_p")]
DHP$group="DHP"
TDM1=results[,c("biomarker","TDM1_OR","TDM1_LCI","TDM1_UCI","TDM1_lr_p")]
TDM1$group="T-DM1"
colname=c("biomarker","OR","LCI","UCI","p","group")
colnames(Total)=colname;colnames(DHP)=colname;colnames(TDM1)=colname
df=rbind(Total,DHP,TDM1)
df$biomarker=gsub("coding_mutation_", "",df$biomarker)
biomarker=df$biomarker
df$OR=as.numeric(df$OR);df$LCI=as.numeric(df$LCI);df$UCI=as.numeric(df$UCI)

df$biomarker=gsub("COSMIC.", "", df$biomarker)
#
df |>
  group_by(group) |>
  forestplot(labeltext=biomarker,
             mean=OR,lower=LCI,upper=UCI,
             zero = 1,
             boxsize = .25, # We set the box size to better visualize the type
             line.margin = .1, # We need to add this to avoid crowding
             lty.ci = c(3),
             clip = c(0,2.5),xlab = " OR with 95% CI") |> 
  fp_add_lines("black") |> 
  fp_add_header("COSMIC Signature") |> 
  fp_set_style(box = c("#DF8F44FF", "#8491B4FF","#91D1C2FF") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) |> 
  fp_set_zebra_style("#F5F9F9")

#7X5 landscape 65%
##################################
########ProfileExtract############
##################################
library(data.table);library(rtracklayer);library(tidyverse)
library(deconstructSigs);library(BSgenome.Hsapiens.UCSC.hg38);library(data.table)
# https://www.nature.com/articles/s41586-019-1056-z#Sec2
maf=fread('E:/Projects/PREDIX_HER2/Validation/brca_tcga/brca_tcga/data_mutations.txt')
table(maf$Chromosome)
maf=filter(maf,Variant_Type=='SNP',Chromosome!="MT",Chromosome!="GL000209.1",Chromosome!='X',Chromosome!="Y")
maf$Chromosome=paste0("chr",maf$Chromosome)
maf$patientID=substr(maf$Tumor_Sample_Barcode,1,12)
her2=fread("E:/Projects/PREDIX_HER2/Validation/brca_tcga/brca_tcga/data_clinical_patient.txt")
her2=her2%>%filter(HER2_FISH_STATUS=="Positive"|IHC_HER2=="Positive")
maf=maf%>%filter(patientID%in%her2$PATIENT_ID)
maf=select(maf,c("Chromosome","Start_Position","End_Position","Tumor_Sample_Barcode","Reference_Allele","Tumor_Seq_Allele2"))
#write.table(maf,file="E:/Projects/PREDIX_HER2/Multimodal/Analyses/SigProfilerExtractor/tcga_her2.bed",quote = F,row.names =F,sep="\t")
maf=fread("E:/Projects/PREDIX_HER2/Multimodal/Analyses/SigProfilerExtractor/tcga_her2_hg38.bed")
str(maf)
mutect=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/Mutect2_tumor_normal_baseline_PREDIX_HER2.PASS.oncokb_annotated_curated.txt")
mutect=mutect%>%filter(Variant_Type=='SNP')
mutect=select(mutect,c("Chromosome","Start_Position","End_Position","Tumor_Sample_Barcode","Reference_Allele","Tumor_Seq_Allele2"))
str(mutect)
df=rbind(mutect,maf)
table(df$Tumor_Seq_Allele2)
sigs.input <- mut.to.sigs.input(mut.ref = df, 
                                sample.id = "Tumor_Sample_Barcode", 
                                chr = "Chromosome", 
                                pos = "Start_Position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele2",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)
cosine=t(sigs.input)%>%as.data.frame()
cosine$'Mutation Types'=row.names(cosine)
cosine <- cosine[, c(ncol(cosine), 1:(ncol(cosine)-1))]
write.table(cosine,file='E:/Projects/PREDIX_HER2/Multimodal/Analyses/SigProfilerExtractor/PREDIX_HER2_tcga_cosine96.txt',quote = F,row.names =F,sep="\t")
##################################
#############FigureS2D############
##################################
library(GGally)
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.txt")
genomic=left_join(genomic,clin,by="patientID")
colnames(genomic)
variable=c("totalTMB","TMB_uniform","TMB_clone","TMB_subclone","TMB_oncogenic")
ggpairs(genomic, columns =variable, aes(color = Response, alpha = 0.5),
        lower = list(continuous = "smooth"))

#6X6 60%

##################################
#############FigureS2E############
##################################
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
genomic=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.txt")
#genomic$patientID=as.character(genomic$patientID)
genomic=left_join(genomic,clin,by="patientID")
genomic$manuTILS
d=genomic%>%select(c("coding_mutation_TP53_oncokb","ER","manuTILS"))
d <- reshape2::melt(d,id.vars=c("coding_mutation_TP53_oncokb","ER"))
d$value=as.numeric(d$value)
d$coding_mutation_TP53_oncokb[d$coding_mutation_TP53_oncokb=="1"]='Mutation'
d$coding_mutation_TP53_oncokb[d$coding_mutation_TP53_oncokb=="0"]='Wild-type'
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
Fig <-
  ggplot(d,aes(x=coding_mutation_TP53_oncokb,y=value,fill=coding_mutation_TP53_oncokb))+
  geom_boxplot(outlier.size = 0.5, width=0.7)+
  facet_wrap(variable~ER,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Treatment Response")+
  stat_compare_means(aes(group=coding_mutation_TP53_oncokb),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.98)+
  labs(y="TILs (%)",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))
Fig
#6X6 60%
##################################
#############FigureS2E############
#########TMB vs TILs by pCR#######
##################################
library(ggpubr)
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
genomic=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.txt")
#genomic$patientID=as.character(genomic$patientID)
genomic=left_join(genomic,clin,by="patientID")
# Extending the regression line --> fullrange = TRUE
# Add marginal rug (marginal density) ---> rug = TRUE
genomic$TMB_clone=as.numeric(genomic$TMB_clone)
genomic$manuTILS=as.numeric(genomic$manuTILS)
ggscatter(genomic, x ="TMB_clone", y ="manuTILS",
          add = "reg.line",                         # Add regression line
          color = "Response", palette = "jco",           # Color by groups "cyl"
          shape = "Response",                            # Change point shape by groups "cyl"
          fullrange = TRUE,                         # Extending the regression line
          rug = TRUE                                # Add marginal rug
)+stat_cor(aes(color = Response))           # Add correlation coefficient





