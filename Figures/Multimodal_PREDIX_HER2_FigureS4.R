##########FigureS4a#########
library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
library(forestplot)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
genomic=left_join(genomic,clin,by="patientID")
variable=c("coding_mutation_MAPK_ERK_pathway","coding_mutation_PIK3_AKT_pathway","coding_mutation_CDK_RB_pathway")
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
df$biomarker=gsub("_oncokb", "",df$biomarker)
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
#6X6 landscape 65%
##########FigureS4b Mut Sig#########
library(deconstructSigs);library(BSgenome.Hsapiens.UCSC.hg38);library(data.table);library(tidyverse)
# https://www.nature.com/articles/s41586-019-1056-z#Sec2
maf=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/Mutect2_tumor_normal_baseline_PREDIX_HER2.PASS.oncokb_annotated_vaf0.5_curated.rds")
maf$vaf=maf$t_alt_count/maf$t_depth
maf=maf%>%filter(Variant_Type=="SNP")
freq=maf%>%group_by(Tumor_Sample_Barcode)%>%summarise(n = n())
maf=maf[maf$Tumor_Sample_Barcode%in%freq$Tumor_Sample_Barcode[freq$n>10],]
sigs.input <- mut.to.sigs.input(mut.ref = maf, 
                                sample.id = "Tumor_Sample_Barcode", 
                                chr = "Chromosome", 
                                pos = "Start_Position", 
                                ref = "Reference_Allele", 
                                alt = "Tumor_Seq_Allele2",
                                bsg = BSgenome.Hsapiens.UCSC.hg38)
cosine=t(sigs.input)%>%as.data.frame()
cosine$'Mutation Types'=row.names(cosine)
cosine <- cosine[, c(ncol(cosine), 1:(ncol(cosine)-1))]
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
w$sampleID=row.names(w)
genomic=w[,c("sampleID","Signature.2","Signature.3","Signature.6","Signature.7","Signature.10","Signature.13")]
colnames(genomic)[2:7] <- paste0("COSMIC.",colnames(genomic)[2:8])
genomic$mut_sig=rowSums(genomic[,2:7])
genomic=genomic[order(-genomic$mut_sig),]
given_order=genomic$sampleID
sig=genomic[,c("sampleID","COSMIC.Signature.2","COSMIC.Signature.3","COSMIC.Signature.6","COSMIC.Signature.7",              
               "COSMIC.Signature.10","COSMIC.Signature.13")]
long_data <- sig%>%gather(key = "signature", value = "value", -sampleID)
long_data=long_data%>%filter(value!=0)
# 计算每个 patientID 的总和
sum_values <- aggregate(value ~ sampleID, data = long_data, sum)

# 根据总和值对 patientID 进行排序
sorted_patientIDs <- sum_values[order(sum_values$value, decreasing = TRUE), ]$sampleID

# 将 patientID 转换为 factor，并按照排序后的顺序重新赋值
long_data$sampleID <- factor(long_data$sampleID, levels = sorted_patientIDs)

library("scales");library(ggsci)
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
ggplot(long_data, aes(fill = signature, y = value, x =sampleID)) + 
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

##########FigureS4c Mut Sig#########
library(tidyverse);library(data.table);library(pheatmap)
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
genomic=left_join(genomic,clin,by="patientID")
var=c("CX1","CX2","CX3","CX4","CX5","CX9")
df=genomic[,var];row.names(df)=genomic$sampleID;df=t(df)
ann_df <-genomic[,c("sampleID","Arm","ER")]
row.names(ann_df) <- ann_df$sampleID
pheatmap(df)

##########FigureS4c cor matrix#########
library(corrplot)
genomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
var=c("COSMIC.Signature.2","COSMIC.Signature.3","COSMIC.Signature.6","COSMIC.Signature.7",              
      "COSMIC.Signature.10","COSMIC.Signature.13","LOH_Del_burden","CNV_burden","TMB_uniform",
      "TMB_clone","HRD")
df=genomic[,var]
#M=cor(df,na.rm=T)
M <- cor(df, use = "pairwise.complete.obs")
testRes = cor.mtest(df, conf.level = 0.95)

corrplot(M, p.mat = testRes$p, tl.col = "black",method = 'circle', type = 'lower', insig='blank',
         order = 'AOE', diag = FALSE,col =c("#053061","#246BAE","#549EC9","#A7CFE4", "#E5F0F6","#FDEBDF","#F7B698","#DC6F58","#B51F2E","#67001F"))$corrPos -> p1
text(p1$x, p1$y, round(p1$corr, 2))
# 25X25 20%
##########FigureS4d#########
library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
library(forestplot);library(ggpubr)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
data=left_join(genomic,clin,by="patientID")
variable0=c("COSMIC.Signature.2","COSMIC.Signature.3","COSMIC.Signature.6",
            "COSMIC.Signature.7","COSMIC.Signature.10","COSMIC.Signature.13")
data[,variable0]=data[,variable0] %>% mutate(across(where(is.numeric), scale))
df=data[!is.na(data$COSMIC.Signature.2),]
results0=Logistic_batch_adjER(df,"pCR","Arm",variable0,"ER")%>%as.data.frame()

variable=c("TMB_uniform","TMB_clone","CNV_burden","LOH_Del_burden","HRD")
data[,variable]=data[,variable] %>% mutate(across(where(is.numeric), scale))
df=data
results=Logistic_batch_adjER(df,"pCR","Arm",variable,"ER")%>%as.data.frame()
results=rbind(results,results0)

results$whole_lr_p=as.numeric(results$whole_lr_p)
results$log10P=-log10(results$whole_lr_p)
results$biomarker=factor(results$biomarker,levels =c("HRD","LOH_Del_burden","CNV_burden","COSMIC.Signature.13","COSMIC.Signature.10","COSMIC.Signature.7","COSMIC.Signature.6",
                                                     "COSMIC.Signature.3","COSMIC.Signature.2","TMB_clone","TMB_uniform"))
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
ggplot(data=results,aes(biomarker,log10P))+
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.7,size=0.25,fill="#916ba6")+
  labs(x="Genomic Metrics",y="-log10p")+
  theme_manuscript(base_size = figure_font_size)+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))+
  coord_flip() 
################################################
#################COSMIC 6/13####################
################################################
library(readxl);library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
library(forestplot);library(ggpubr)
library(deconstructSigs);library(BSgenome.Hsapiens.UCSC.hg19);library(data.table)
meta=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/transneo/41586_2021_4278_MOESM4_ESM.xlsx",sheet=1)%>%as.data.frame()
meta=meta[meta$HER2.status=="POS",]
maf=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/transneo/41586_2021_4278_MOESM4_ESM.xlsx",sheet=2)%>%as.data.frame()
maf=maf[maf$Donor.ID%in%meta$Donor.ID,]
freq=maf%>%group_by(Donor.ID)%>%summarise(n = n())
maf=maf[maf$Donor.ID%in%freq$Donor.ID[freq$n>10],]
meta=meta[meta$Donor.ID%in%unique(maf$Donor.ID),]
maf$Chr <- ifelse(grepl("^chr", maf$Chr), maf$Chr, paste0("chr", maf$Chr))
sigs.input <- mut.to.sigs.input(mut.ref = maf, 
                                sample.id = "Donor.ID", 
                                chr = "Chr", 
                                pos = "Start", 
                                ref = "Ref_Allele", 
                                alt = "Tumor_Alt_Allele",
                                bsg = BSgenome.Hsapiens.UCSC.hg19)
cosine=t(sigs.input)%>%as.data.frame()
cosine$'Mutation Types'=row.names(cosine)
cosine <- cosine[, c(ncol(cosine), 1:(ncol(cosine)-1))]
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
df <- w[, colMeans(w == 0) < 0.90]
ncol(df)
df[df==0]=0.001
# normalize columns using the first column (log2) using apply()
df_norm <- apply(df[,2:8], 2, function(x) log2(x/df$Signature.1))%>%as.data.frame()
df_norm$Donor.ID=row.names(df_norm)
d=left_join(meta,df_norm,by="Donor.ID")
d=d[d$pCR.RD!="NA",]

d=d%>%select(c("pCR.RD","Signature.6","Signature.13"))
d <- reshape2::melt(d,id.vars=c("pCR.RD"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
Figs4e <-
  ggplot(d,aes(x=pCR.RD,y=value,fill=pCR.RD))+
  geom_boxplot(outlier.size = 0.2, width=0.8)+
  facet_wrap(variable~.,scales="free_y",nrow=2)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Treatment Response")+
  stat_compare_means(aes(group=pCR.RD),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.98)+
  labs(y="",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.2,0.2,0.5), "lines"))
Figs4e

d=left_join(meta,df_norm,by="Donor.ID")
d=d[d$pCR.RD!="NA",]
d$pCR=0
d$pCR[d$pCR.RD=="pCR"]=1
res <- glm(as.numeric(pCR) ~ Signature.6 + ER.status, family = "binomial", data = d)
ShowRegTable(res)
res <- glm(as.numeric(pCR) ~ Signature.13 + ER.status, family = "binomial", data = d)
ShowRegTable(res)
################################################
#################TMB clonal ####################
################################################
library(readxl);library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
library(forestplot);library(ggpubr);library(data.table)
meta=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/transneo/41586_2021_4278_MOESM4_ESM.xlsx",sheet=1)%>%as.data.frame()
meta=meta[meta$HER2.status=="POS",]
maf=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/transneo/41586_2021_4278_MOESM4_ESM.xlsx",sheet=2)%>%as.data.frame()
maf=maf[maf$Donor.ID%in%meta$Donor.ID,]
table(maf$MAF_Variant)
vc.nonSilent = c("Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "Translation_Start_Site",
                 "Nonsense_Mutation", "Nonstop_Mutation", "In_Frame_Del",
                 "In_Frame_Ins", "Missense_Mutation")
maf_clone=maf%>%filter(MAF_Variant%in%vc.nonSilent,ClonalStatus=="CLONAL",tumour.vaf>0.05)
freq=maf_clone%>%group_by(Donor.ID)%>%summarise(n = n())%>%as.data.frame()
tmb_clone=data.frame(Donor.ID=freq$Donor.ID,TMB_clone=freq$n/45.54)
meta=left_join(meta,tmb_clone,by="Donor.ID")
meta$TMB_clone[is.na(meta$TMB_clone)]=0

meta=meta[meta$pCR.RD!="NA",]
meta$pCR=0
meta$pCR[meta$pCR.RD=="pCR"]=1
res <- glm(as.numeric(pCR) ~ TMB_clone + ER.status, family = "binomial", data = meta)
ShowRegTable(res)

##########FigureS4f#########
library(tidyverse);library(data.table);library("survival");library(survminer)
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
genomic=left_join(genomic,clin,by="patientID")
genomic=genomic[genomic$Arm=="T-DM1",]
cutoff=quantile(genomic$LOH_Del_burden, probs = c(0.75))
genomic$group="low"
genomic$group[genomic$LOH_Del_burden>0.21]="high"
# 4. Fit survival curves and visualize
fit <- survfit(Surv(EFS.time,EFS.status) ~ group, data = genomic)
ggsurvplot(fit,data = genomic, palette = c("#BC3C29FF","#6F99ADFF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
# 5X5 50%
library(tableone)
table(genomic$group)
cox.test <- coxph(Surv(EFS.time,EFS.status)~group+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=genomic) ##DII_density_with_supp
ShowRegTable(cox.test)
test.ph <- cox.zph(cox.test)
test.ph
##########FigureS4g,h#########
library(Blasso);library(tidyverse)
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
genomic=left_join(genomic,clin,by="patientID")
genomic=genomic[genomic$Arm=="DHP",]
genomic=genomic[genomic$Arm=="T-DM1",] ###
library("survival");library(survminer)
# 1. Determine the optimal cutpoint of variables
res.cut <- surv_cutpoint(genomic, time = "EFS.time", event = "EFS.status",variables = c("COSMIC.Signature.13"))
summary(res.cut)
# 2. Plot cutpoint for DEPDC1
# palette = "npg" (nature publishing group), see ?ggpubr::ggpar
plot(res.cut, "COSMIC.Signature.13", palette = "npg")
# 3. Categorize variables
res.cat <- surv_categorize(res.cut)
head(res.cat)
res.cat$ER=genomic$ER
res.cat$TUMSIZE=genomic$TUMSIZE
res.cat$ANYNODES=genomic$ANYNODES
# 4. Fit survival curves and visualize
fit <- survfit(Surv(EFS.time,EFS.status) ~COSMIC.Signature.13, data = res.cat)
ggsurvplot(fit,data = res.cat, palette = c("#BC3C29FF","#6F99ADFF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
library(tableone)
table(res.cat$COSMIC.Signature.13)
res.cat$COSMIC.Signature.13=factor(res.cat$COSMIC.Signature.13,levels = c("low","high"))
cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(COSMIC.Signature.13)+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=res.cat) ##DII_density_with_supp
ShowRegTable(cox.test)
table(res.cat$COSMIC.Signature.13)

##########FigureS4i#########
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.txt")%>%as.data.frame()
genomic=inner_join(genomic,clin,by="patientID")%>%as.data.frame()
variable=c("LOH_Del_burden")
genomic=slice_metric(genomic,variable,20,80,5)
## re-run the batch logistic
selected_columns <- grep("_per_", colnames(genomic), value = TRUE)
Interact_result=Logistic_batch_adjER(genomic,"pCR","Arm",selected_columns ,"ER")%>%as.data.frame()
Interact_result=separate(Interact_result,biomarker, into = c("group", "per"), sep = "_per_", remove = FALSE)
Interact_result$per=as.numeric(Interact_result$per)
Interact_result$interaction_coefficient=as.numeric(Interact_result$interaction_coefficient)
Interact_result$interaction_coefficient_LCI=as.numeric(Interact_result$interaction_coefficient_LCI)
Interact_result$interaction_coefficient_UCI=as.numeric(Interact_result$interaction_coefficient_UCI)


pd <- position_dodge(3)
ggplot(Interact_result, 
       aes(x = per, 
           y = interaction_coefficient, group=group, color=group)) +
  geom_point(position=pd, 
             size=3) +
  geom_line(position=pd, 
            size = 1) +
  geom_errorbar(aes(ymin = interaction_coefficient_LCI, 
                    ymax = interaction_coefficient_UCI), 
                width = .1, 
                position=pd, 
                size=1) +
  scale_color_brewer(palette="Set1") +
  theme_minimal() +
  labs(title = "",
       subtitle = "(mean +/- standard error)",
       x = "", 
       y = "",
       color = "Genomic Biomarker")
#9X3
##########FigureS4j#########
library(ggpubr)
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.txt")%>%as.data.frame()
genomic=inner_join(genomic,clin,by="patientID")%>%as.data.frame()
cna=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/PREDIX_HER2_CUTseq_gistics2_baseline_gene_curated.txt")%>%as.data.frame()
data=left_join(genomic,cna,by="patientID")
cor.test(data$BRCA2,data$COSMIC.Signature.13)
data$COSMIC.Signature.13=scale(data$COSMIC.Signature.13)

ggscatter(data, x = "BRCA2", y = "COSMIC.Signature.13",
          add = "reg.line",  # Add regressin line
          xlim=c(-1,1),ylim=c(-10,10),size = 0.5,
          add.params = list(color = "#0E4C92", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE # Add confidence interval
) + stat_cor(method = "pearson", label.x = 0.3, label.y = 8) 

ggscatter(data, x = "BRCA2", y = "LOH_Del_burden",
          add = "reg.line",  # Add regressin line
          xlim=c(-1,1),ylim=c(0,0.5),size = 0.5,
          add.params = list(color = "#0E4C92", fill = "lightgray"), # Customize reg. line
          conf.int = TRUE # Add confidence interval
) + stat_cor(method = "pearson", label.x = 0.3, label.y = 0.5) 

#3.5X3.5
##################################
###########FigureS4A,B############
##################################
library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
library(forestplot)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.txt")
variable=colnames(genomic)[3:ncol(genomic)]
genomic=left_join(genomic,clin,by="patientID")%>%as.data.frame()
colnames(genomic)
#######ER-negative
results=Logistic_batch(genomic[genomic$ER=="negative",],"pCR","Arm",variable)%>%as.data.frame()
colnames(results)
biomarker=union(results$biomarker[results$TDM1_lr_p<0.106],results$biomarker[results$DHP_lr_p<0.106])
TDM1=results[results$biomarker%in%biomarker,c("biomarker","TDM1_OR","TDM1_lr_p")]
TDM1$group="TDM1"
DHP=results[results$biomarker%in%biomarker,c("biomarker","DHP_OR","DHP_lr_p")]
DHP$group="DHP"
colnames(TDM1)=c("Genomic","OR","Pvalue","group")
colnames(DHP)=c("Genomic","OR","Pvalue","group")
df=rbind(TDM1,DHP)
df$Pvalue=as.numeric(df$Pvalue)
df$log10P=-log10(df$Pvalue)
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
ggplot(data=df,aes(Genomic,log10P,fill=group))+
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.7,size=0.25)+
  scale_fill_Arm(name="Treatment Arm")+
  labs(x="Genomic Metrics",y="-log10(p value)")+
  theme_manuscript(base_size = figure_font_size)+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))+
  coord_flip()

#######ER-positive
results=Logistic_batch(genomic[genomic$ER=="positive",],"pCR","Arm",variable)%>%as.data.frame()
colnames(results)
biomarker=union(results$biomarker[results$TDM1_lr_p<0.106],results$biomarker[results$DHP_lr_p<0.106])
TDM1=results[results$biomarker%in%biomarker,c("biomarker","TDM1_OR","TDM1_lr_p")]
TDM1$group="TDM1"
DHP=results[results$biomarker%in%biomarker,c("biomarker","DHP_OR","DHP_lr_p")]
DHP$group="DHP"
colnames(TDM1)=c("Genomic","OR","Pvalue","group")
colnames(DHP)=c("Genomic","OR","Pvalue","group")
df=rbind(TDM1,DHP)
df$Pvalue=as.numeric(df$Pvalue)
df$log10P=-log10(df$Pvalue)
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
ggplot(data=df,aes(Genomic,log10P,fill=group))+
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.7,size=0.25)+
  scale_fill_Arm(name="Treatment Arm")+
  labs(x="Genomic Metrics",y="-log10(p value)")+
  theme_manuscript(base_size = figure_font_size)+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))+
  coord_flip()








