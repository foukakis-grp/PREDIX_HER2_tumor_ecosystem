#FigureS7a
library(data.table);library(ggplot2)
d=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.rds")
d=d[,c("Her2ADC","G0scores")]
d <- reshape2::melt(d,id.vars=c("Her2ADC"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
eFig7a <-
  ggplot(d,aes(x=Her2ADC,y=value,fill=Her2ADC))+
  geom_boxplot(outlier.size = 0.5, width=0.7)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_HER2ADC(name="Her2ADC")+
  stat_compare_means(aes(group=Her2ADC),label = "p.format", hide.ns = F,size=5,
                     color="black", label.y.npc = 0.94)+
  labs(y="Normalized G0 arrest score",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))
eFig7a 

#FigureS7b
library(tidyverse);library(data.table);library(ggpubr);library(ggplot2);library(ggridges);library(coin)
abundance=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/TMM-normalized-TPM.rds")
colnames(abundance)=substr(colnames(abundance),9,12)
rna=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.rds")
abundance=abundance[,rna$patientID]
res=cbind(rna,t(abundance[c("NR2F1","NFE2L2","GAS6","CD8A","MKI67"),]))
d=res[,c("Her2ADC","NR2F1","NFE2L2","GAS6")]
d <- reshape2::melt(d,id.vars=c("Her2ADC"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
eFig7b <-
  ggplot(d,aes(x=Her2ADC,y=value,fill=Her2ADC))+
  geom_boxplot(outlier.size = 0.5, width=0.7)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_HER2ADC(name="Her2ADC")+
  stat_compare_means(aes(group=Her2ADC),label = "p.format", hide.ns = F,size=5,
                     color="black", label.y.npc = 0.94)+
  labs(y="TMM-normalized TPM (log2)",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))
eFig7b 

ggarrange(eFig7a,eFig7b,widths = c(1,2.5),nrow = 1,common.legend = TRUE)
# 10X5
#FigureS7c forest
library(tidyverse)
library(tableone)
library(forestploter)
require(openxlsx)
library(data.table)
d=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.rds")
d$patientID=as.integer(d$patientID)
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
d=left_join(d,clin,by="patientID")
table(d$Arm,d$Her2ADC)
res<- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = d[d$Her2ADC=="S3",])
ShowRegTable(res)

interaction_1<- glm(as.numeric(pCR) ~ Her2ADC+Arm+ER, family = "binomial", data = d)
interaction_2<- glm(as.numeric(pCR) ~ Her2ADC*Arm+Her2ADC+Arm+ER, family = "binomial", data = d)
interaction_lr <- lrtest(interaction_1,interaction_2)
interaction_lr

df=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS7/Subgroup_Forestplot.csv")
df$DHP <- ifelse(is.na(df$DHP), "", df$DHP)
df$`T-DM1` <- ifelse(is.na(df$`T-DM1`), "", df$`T-DM1`)
df$`P for interaction` <- ifelse(is.na(df$`P for interaction`), "", df$`P for interaction`)
df$` ` <- paste(rep(" ", 20), collapse = " ")
a=df$`OR (95% CI)`;b=df$`P for interaction`
df$`OR (95% CI)`=NULL;df$`P for interaction`=NULL
df$`OR (95% CI)`=a
df$`P for interaction`=b

tm <- forest_theme(base_size = 10,
                   refline_col = "black",
                   arrow_type = "closed",
                   footnote_col = "blue")

p <- forest(df[,c(1:3,7:9)],
            est = df$OR,
            lower = df$Low, 
            upper = df$High,
            ci_column = 4,
            ref_line = 1,
            arrow_lab = c("DHP Better", "T-DM1 Better"),
            xlim = c(0, 4),
            ticks_at = c(0,0.5, 1, 2, 3),
            theme = tm)

# Print plot
plot(p)

#7X5
ggplot2::ggsave(filename = "E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS7/FigureS7c.pdf", plot = p,
                dpi = 300,
                width = 7.5, height = 5, units = "in")

#FigureS7d,e
library("survival");library(survminer)
d=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.rds")
d$patientID=as.integer(d$patientID)
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
df=left_join(d,clin,by="patientID")
fit <- survfit(Surv(EFS.time,EFS.status) ~Her2ADC, data = df)
ggsurvplot(fit,data = df, palette = c("#4DE66699","#F39B7FFF","#4D1A6699"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
# 5X5 
fit <- survfit(Surv(EFS.time,EFS.status) ~Her2ADC+Arm, data = df)
ggsurvplot(fit,data = df, palette = "npg",
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))

# 5X5.5 P

# 5X5 
fit <- survfit(Surv(EFS.time,EFS.status) ~Arm, data = df[df$Her2ADC=="S3",])
ggsurvplot(fit,data = df, palette = "npg",
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))

# 5X5.5 P
res.cox <- coxph(Surv(EFS.time,EFS.status) ~ Arm+ strata(as.factor(ER))+as.factor(ANYNODES)+as.factor(TUMSIZE), data = df[df$Her2ADC=="S3",])
cox.zph(res.cox)
ShowRegTable(res.cox)

#TCGA f,g,h
library(ggplot2)
d=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Validation/TCGA/TCGA_transcriptomic.rds")
df=table(d$Her2ADC,d$sspbc.subtype)%>%as.data.frame()
colnames(df)=c("Her2ADC","Intrinsic_subtype","Freq")
f=ggplot(data =  df, mapping = aes(x = Her2ADC, y = Intrinsic_subtype)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
  scale_fill_gradient(low="white", high="#009194") +
  theme_bw()
f

df=read_excel('E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS7/HER2ADC_GSEA_TCGA.xlsx')
df$PVAL=-log10(df$padj)
S1_pathway=df$pathway[df$HER2ADC=='S1'&df$PVAL>1.3&df$NES>0]%>%as.character()
S2_pathway=df$pathway[df$HER2ADC=='S2'&df$PVAL>1.3&df$NES>0]%>%as.character()
S3_pathway=df$pathway[df$HER2ADC=='S3'&df$PVAL>1.3&df$NES>0]%>%as.character()
pathway=c(S1_pathway,S2_pathway,S3_pathway)
df$pathway <- gsub("HALLMARK_","",df$pathway)
df$pathway <- gsub("_"," ",df$pathway)
df$pathway <- gsub("PI3K AKT MTOR","PI3K/AKT/mTOR",df$pathway)
pathway=c('MTORC1 SIGNALING',"G2M CHECKPOINT",
          "MYC TARGETS V1",'PI3K/AKT/mTOR SIGNALING',
          "OXIDATIVE PHOSPHORYLATION","GLYCOLYSIS","FATTY ACID METABOLISM","XENOBIOTIC METABOLISM",  ## S1
          "ESTROGEN RESPONSE EARLY","ESTROGEN RESPONSE LATE","DNA REPAIR",  ##S2
          "MYOGENESIS","KRAS SIGNALING UP","APICAL JUNCTION") ##S3
df=df[df$pathway%in%pathway,]
head(df)
df$pathway=factor(df$pathway,levels =c(pathway))
# #buble plot
g <- ggplot(df, aes(x = HER2ADC, y = pathway, size = PVAL, color = NES)) +
  geom_point(shape=19) +
  #ggtitle("pathway heterogeneity") +
  labs(x = NULL, y = NULL,
       size = "-log10 pvalue", color = "NES") +
  scale_size(range = c(0,4)) +
  scale_color_gradient2( low = "#375E97",mid = "white", high = "#FB6542") +
  #scale_color_gradient2(low="red",mid="white",high="blue",midpoint = 1) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(colour="black",size=9),
        axis.line = element_line(size=0.3, colour = "black"),
        #panel.grid.major = element_line(colour = "#d3d3d3"),
        #panel.grid.minor = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.3),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 9,hjust=1,vjust=0.5),
        axis.text.y=element_text(colour="black", size = 9)) +
  theme(plot.margin = unit(rep(1,4),"lines"))
g

TCGA_transcriptomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Validation/TCGA/TCGA_transcriptomic.rds")
d=TCGA_transcriptomic[,c("Her2ADC","G0scores")]
median_data <- mean(d$G0scores)
sd_data <- sd(d$G0scores)
lower_bound <- median_data - 1.5 * sd_data
upper_bound <- median_data + 1.5 * sd_data
d=d[d$G0scores>lower_bound,]
d <- reshape2::melt(d,id.vars=c("Her2ADC"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
h <-
  ggplot(d,aes(x=Her2ADC,y=value,fill=Her2ADC))+
  geom_boxplot(outlier.size = 0.5, width=0.7)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_HER2ADC(name="Her2ADC")+
  stat_compare_means(aes(group=Her2ADC),label = "p.format", hide.ns = F,size=5,
                     color="black", label.y.npc = 0.94)+
  labs(y="Normalized G0 arrest score",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"),
        legend.position = "none")
h

ggarrange(f,g,h,nrow = 1,widths = c(1.2,1.2,0.7))
# 16X5 50%   f 85%


#SCANB i,j,k
library(ggplot2)
d=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Validation/SCAN-B/SCANB_transcriptomic.rds")
d$SSP.Subtype=factor(d$SSP.Subtype,levels = c("LumA","LumB","Her2","Basal"))
df=table(d$Her2ADC,d$SSP.Subtype)%>%as.data.frame()
colnames(df)=c("Her2ADC","Intrinsic_subtype","Freq")
f=ggplot(data =  df, mapping = aes(x = Her2ADC, y = Intrinsic_subtype)) +
  geom_tile(aes(fill = Freq), colour = "white") +
  geom_text(aes(label = sprintf("%1.0f", Freq)), vjust = 1) +
  scale_fill_gradient(low="white", high="#009194") +
  theme_bw()
f

df=read_excel('E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS7/HER2ADC_GSEA_SCANB.xlsx')
df$PVAL=-log10(df$padj)
S1_pathway=df$pathway[df$HER2ADC=='S1'&df$PVAL>1.3&df$NES>0]%>%as.character()
S2_pathway=df$pathway[df$HER2ADC=='S2'&df$PVAL>1.3&df$NES>0]%>%as.character()
S3_pathway=df$pathway[df$HER2ADC=='S3'&df$PVAL>1.3&df$NES>0]%>%as.character()
pathway=c(S1_pathway,S2_pathway,S3_pathway)
df$pathway <- gsub("HALLMARK_","",df$pathway)
df$pathway <- gsub("_"," ",df$pathway)
df$pathway <- gsub("PI3K AKT MTOR","PI3K/AKT/mTOR",df$pathway)
pathway=c("MITOTIC SPINDLE",'MTORC1 SIGNALING','E2F TARGETS',"G2M CHECKPOINT",
          "MYC TARGETS V1",'PI3K/AKT/mTOR SIGNALING',
          "OXIDATIVE PHOSPHORYLATION","GLYCOLYSIS",  ## S1
          "ESTROGEN RESPONSE EARLY","ESTROGEN RESPONSE LATE","DNA REPAIR",  ##S2
          "MYOGENESIS","UV RESPONSE DN","KRAS SIGNALING UP",
          "ADIPOGENESIS","APICAL JUNCTION") ##S3
df=df[df$pathway%in%pathway,]
head(df)
df$pathway=factor(df$pathway,levels =c(pathway))
# #buble plot
g <- ggplot(df, aes(x = HER2ADC, y = pathway, size = PVAL, color = NES)) +
  geom_point(shape=19) +
  #ggtitle("pathway heterogeneity") +
  labs(x = NULL, y = NULL,
       size = "-log10 pvalue", color = "NES") +
  scale_size(range = c(0,4)) +
  scale_color_gradient2( low = "#375E97",mid = "white", high = "#FB6542") +
  #scale_color_gradient2(low="red",mid="white",high="blue",midpoint = 1) +
  theme(legend.position = "top", legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key.size = unit(0.5, "cm"),
        legend.text = element_text(colour="black",size=9),
        axis.line = element_line(size=0.3, colour = "black"),
        #panel.grid.major = element_line(colour = "#d3d3d3"),
        #panel.grid.minor = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.3),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 9,hjust=1,vjust=0.5),
        axis.text.y=element_text(colour="black", size = 9)) +
  theme(plot.margin = unit(rep(1,4),"lines"))
g

TCGA_transcriptomic=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Validation/SCAN-B/SCANB_transcriptomic.rds")
d=TCGA_transcriptomic[,c("Her2ADC","G0scores")]
d <- reshape2::melt(d,id.vars=c("Her2ADC"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
h <-
  ggplot(d,aes(x=Her2ADC,y=value,fill=Her2ADC))+
  geom_boxplot(outlier.size = 0.5, width=0.7)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_HER2ADC(name="Her2ADC")+
  stat_compare_means(aes(group=Her2ADC),label = "p.format", hide.ns = F,size=5,
                     color="black", label.y.npc = 0.94)+
  labs(y="Normalized G0 arrest score",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"),
        legend.position = "none")
h

ggarrange(f,g,h,nrow = 1,widths = c(1.2,1.2,0.7))
# 16X5 50%   f 85%

#SCANB-survival
library(survminer);library(survival);library(tableone)
res=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Validation/SCAN-B/SCANB_transcriptomic.rds")
# Overall survival (OS), recurrence-free interval (RFi), and breast cancer-free interval (BCFi)
res$RFi_days=res$RFi_days/30
res$DRFi_days=res$DRFi_days/30
res$BCFi_days=res$BCFi_days/30
res$OS_days=res$OS_days/30

fit <- survfit(Surv(RFi_days, RFi_event) ~Her2ADC,data = res)
ggsurvplot(fit,data = res, palette = c("#4DE66699","#F39B7FFF","#4D1A6699"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 100))


fit <- survfit(Surv(DRFi_days, DRFi_event) ~Her2ADC,data = res)
ggsurvplot(fit,data = res, palette = c("#4DE66699","#F39B7FFF","#4D1A6699"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 100))

fit <- survfit(Surv(BCFi_days, BCFi_event) ~Her2ADC,data = res)
ggsurvplot(fit,data = res, palette = c("#4DE66699","#F39B7FFF","#4D1A6699"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 90))

fit <- survfit(Surv(OS_days, OS_event) ~Her2ADC,data = res)
ggsurvplot(fit,data = res, palette = c("#4DE66699","#F39B7FFF","#4D1A6699"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 120))
# 4X4 50%
res$Her2ADC=factor(res$Her2ADC,levels = c("S1","S2","S3"))

res.cox <- coxph(Surv(RFi_days, RFi_event) ~ as.factor(Her2ADC)+as.factor(LN)+as.factor(T.size)+as.factor(ER), data = res)
cox.zph(res.cox)
ShowRegTable(res.cox)

res.cox <- coxph(Surv(DRFi_days, DRFi_event) ~ as.factor(Her2ADC)+as.factor(LN)+as.factor(T.size)+as.factor(ER), data = res)
cox.zph(res.cox)
ShowRegTable(res.cox)

res.cox <- coxph(Surv(BCFi_days, BCFi_event) ~ as.factor(Her2ADC)+as.factor(LN)+as.factor(T.size)+as.factor(ER), data = res)
cox.zph(res.cox)
ShowRegTable(res.cox)

res.cox <- coxph(Surv(OS_days, OS_event) ~ as.factor(Her2ADC)+as.factor(LN)+as.factor(T.size)+as.factor(ER), data = res)
cox.zph(res.cox)
ShowRegTable(res.cox)
