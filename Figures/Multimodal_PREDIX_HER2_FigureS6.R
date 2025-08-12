# figS6a
library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
library(forestplot)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.txt")
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
df=left_join(df,clin,by="patientID")
df$sspbc.subtype=factor(df$sspbc.subtype,levels = c("Her2","LumA","LumB","Basal"))
res<- glm(as.numeric(pCR) ~ sspbc.subtype, family = "binomial", data = df)
ShowRegTable(res)
Total=data.frame(biomarker=c("LumA","LumB","Basal"),OR=c(0.28,0.16,0.16),
                 LCI=c(0.12,0.06,0.02),UCI=c(0.61,0.37,0.66),
                 p=c(0.002,0.001,0.02),group=c("All","All","All"))

res<- glm(as.numeric(pCR) ~ sspbc.subtype, family = "binomial", data = df[df$Arm=="DHP",])
ShowRegTable(res)
DHP=data.frame(biomarker=c("LumA","LumB","Basal"),OR=c(0.23,0.20,0.14),
                 LCI=c(0.07,0.06,0.01),UCI=c(0.66,0.62,1.04),
                 p=c(0.008,0.007,0.09),group=c("DHP","DHP","DHP"))

res<- glm(as.numeric(pCR) ~ sspbc.subtype, family = "binomial", data = df[df$Arm=="T-DM1",])
ShowRegTable(res)
TDM1=data.frame(biomarker=c("LumA","LumB","Basal"),OR=c(0.34,0.10,0.17),
               LCI=c(0.10,0.01,0.01),UCI=c(1.11,0.40,1.26),
               p=c(0.08,0.004,0.13),group=c("T-DM1","T-DM1","T-DM1"))

df=rbind(Total,DHP,TDM1)
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
             clip = c(0,2.5),xlab = " OR with 95% CI") |> 
  fp_add_lines("black") |> 
  fp_add_header("COSMIC Signature") |> 
  fp_set_style(box = c("#DF8F44FF", "#8491B4FF","#91D1C2FF") |> lapply(function(x) gpar(fill = x, col = "#555555")),
               default = gpar(vertices = TRUE)) |> 
  fp_set_zebra_style("#F5F9F9")

# figS6b
library("survival");library(survminer);library(data.table);library(tidyverse);library(tableone)
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.txt")
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
df=left_join(df,clin,by="patientID")
df$sspbc.subtype=factor(df$sspbc.subtype,levels = c("Basal","LumA","LumB","Her2"))
fit <- survfit(Surv(EFS.time,EFS.status) ~ sspbc.subtype, data = df)
summary(fit)
ggsurvplot(fit, palette = c("#e31a1c","#1f78b4","#a6cee3","#fb9a99"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))

cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(sspbc.subtype)+Arm+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=df) ##DII_density_with_supp
ShowRegTable(cox.test)

# figS6c
library("survival");library(survminer);library(data.table);library(tidyverse);library(tableone)
fusion=readRDS('E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/RNA_fusion_filtered.rds')
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
df=left_join(fusion,clin,by="patientID")
# freq of chromsome peak
# 1. 计算左端和右端染色体的频数
common_fusion=fusion[fusion$left_chr==fusion$right_chr,]
fusion=fusion[fusion$left_chr!=fusion$right_chr,]
left_freq <- table(fusion$left_chr)
right_freq <- table(fusion$right_chr)
common_freq=table(common_fusion$left_chr)
# 2. 合并频数表格（以确保两者具有相同的染色体名称）
all_chromosomes <- union(union(names(left_freq), names(right_freq)),names(common_freq))
left_freq <- as.numeric(left_freq[all_chromosomes])
left_freq[is.na(left_freq)]=0
right_freq <- as.numeric(right_freq[all_chromosomes])
right_freq[is.na(right_freq)]=0
common_freq <- as.numeric(common_freq[all_chromosomes])
common_freq[is.na(common_freq)]=0
# 3. 计算总频数
total_freq <- left_freq + right_freq + common_freq
names(total_freq) <- all_chromosomes
total_freq=as.data.frame(total_freq)
total_freq$chr=row.names(total_freq)
total_freq$chr=paste0("chr",total_freq$chr)
total_freq=total_freq[total_freq$chr!="chrXp",]
total_freq=total_freq[total_freq$chr!="chrXq",]
total_freq$chr=factor(total_freq$chr,levels =c('chr1',"chr2","chr3","chr4","chr5","chr6","chr7","chr8",
                                           "chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16",
                                           "chr17","chr18","chr19","chr20","chr21","chr22") )

library(ggpubr)
ggdotchart(total_freq, x = "chr", y = "total_freq",
           sorting = "ascending",                        
           add = "segments",                             
           xlab="", 
           ylab="Frequency",
           rotate = TRUE,
           dot.size = 3)
# 4X4
# figS6d chord plot
library(circlize)
library(viridis)
library(reshape2)
library(data.table)
library(tidyverse)
fusion=readRDS('E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/RNA_fusion_filtered.rds')
fusion=fusion[fusion$RNA_fusion%in%c("CCDC32--CBX3","ERBB2--CDK12","ERBB2--GRB7","KANSL1--ARL17A","KANSL1--ARL17B",
                                     "LHFPL5--CLPSL1","TFG--ADGRG7","TTC6--MIPOL1"),]
df=data.frame(from=fusion$left_gene,to=fusion$right_gene,value=1)
colnames(df)<-c('from','to','value')
#排序
# 颜色主题方案
mycolor <- viridis(14, alpha = 1, begin = 0, end = 1, option = "C")# fill N of gene_shared
names(mycolor) <-unique(union(df$from,df$to))
circos.clear()
circos.par(start.degree = 90, gap.degree = 0.5, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))
chordDiagram(
  x = df,
  grid.col = mycolor,
  transparency = 0.25,
  directional = 1,
  direction.type = c("arrows", "diffHeight"),
  diffHeight = -0.04,
  annotationTrack = "grid",
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow",
  link.sort = TRUE,
  link.largest.ontop = TRUE)
# 添加数据标签和坐标轴
circos.trackPlotRegion(
  track.index = 1,
  bg.border = NA,
  panel.fun = function(x, y) {
    xlim = get.cell.meta.data("xlim")
    ylim = get.cell.meta.data("ylim")
    sector.index = get.cell.meta.data("sector.index")
    # 添加数据标签
    circos.text(
      x = mean(xlim),
      y = ylim[1] + 1,
      labels = sector.index,
      facing = "clockwise",
      cex = 0.5,niceFacing = TRUE, adj = c(0,0.5)
    )
  })

# figS6g
library("survival");library(survminer);library(data.table);library(tidyverse);library(tableone)
rna=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.txt')
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
df=left_join(rna,clin,by="patientID")

fit <- survfit(Surv(EFS.time,EFS.status) ~ ERBB2_fusion, data = df)
summary(fit)
ggsurvplot(fit,data = df, palette = c("#BC3C29FF","#6F99ADFF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))

fit <- survfit(Surv(EFS.time,EFS.status) ~ chr17q12_fusion, data = df)
summary(fit)
ggsurvplot(fit,data = df, palette = c("#BC3C29FF","#6F99ADFF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))

# figS6h
library("survival");library(survminer);library(data.table);library(tidyverse);library(tableone);
library(ComplexHeatmap);library(circlize);library(RColorBrewer);library(ggstatsplot)
rna=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.txt')
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
df=left_join(rna,clin,by="patientID")
df$patientID=paste0("P",df$patientID)
df$sspbc.subtype=factor(df$sspbc.subtype,levels = c("LumA","LumB","Her2","Basal"))
mydata=df[order(df$sspbc.subtype),]%>%as.data.frame()
row.names(mydata)=mydata$patientID

cellular_process=c("ABC_transporter","Apoptosis","Lysosome","Endocytosis","EMT","Exosome")
metabolism=c("Oxidative_phosphorylation","Purine_metabolism","Citrate_cycles","Glutathione_metabolism",
             "Fatty_acid_metabolism","Glycolysis","Hypoxia")
name=c(cellular_process,metabolism)
mydata[,name]=scale(mydata[,name])
TME_score=as.data.frame(t(mydata[,name]))

mat1 = as.matrix(TME_score[,mydata$patientID[mydata$sspbc.subtype=="LumA"]])   # baseline samples
TME_cluster1=mydata[mydata$patientID%in%colnames(mat1),]
mat2 = as.matrix(TME_score[,mydata$patientID[mydata$sspbc.subtype=="LumB"]])   # baseline samples
TME_cluster2=mydata[mydata$patientID%in%colnames(mat2),]
mat3 = as.matrix(TME_score[,mydata$patientID[mydata$sspbc.subtype=="Her2"]])   # baseline samples
TME_cluster3=mydata[mydata$patientID%in%colnames(mat3),]
mat4 = as.matrix(TME_score[,mydata$patientID[mydata$sspbc.subtype=="Basal"]])   # baseline samples
TME_cluster4=mydata[mydata$patientID%in%colnames(mat4),]

TME_function=c(rep("Cellular_process",times=6),
               rep("Metabolism",times=7))

pheno=TME_cluster1
ha1=HeatmapAnnotation(foo=anno_block(gp= gpar(fill = "#1f78b4"),labels = c("LumA (n = 36)")),
                      Taxane_response=anno_barplot(pheno$Taxane_response,ylim =c(-4,1),border =F,show_name = FALSE,bar_width = 0.4),
                      HER2DX=anno_barplot(pheno$HER2DX_pCR_likelihood_score,ylim =c(0,12),border =F,show_name = FALSE,bar_width = 0.4),
                      PIK3CA=anno_barplot(pheno$pik3ca_sig,ylim =c(0.5,2),border =F,show_name = FALSE,bar_width = 0.4),
                      Arm=pheno$Arm,
                      Response=pheno$Response,
                      Recurrence=pheno$RFS.status,
                      chr17q12_fusion=pheno$chr17q12_fusion,
                      ERBB2_fusion=pheno$ERBB2_fusion,
                      col=list(Arm=c("DHP"="#8491B4FF","T-DM1"="#91D1C2FF"),
                               Response=c("RD"="#525252","pCR"="#f0f0f0"),
                               Recurrence=c("1"="#525252","0"="#f0f0f0"),
                               chr17q12_fusion=c("1"="#525252","0"="#f0f0f0"),
                               ERBB2_fusion=c("1"="#525252","0"="#f0f0f0")),
                      na_col="grey",show_annotation_name = FALSE,
                      border=TRUE)

pheno=TME_cluster2
ha2=HeatmapAnnotation(foo=anno_block(gp= gpar(fill = "#a6cee3"),labels = c("LumB (n = 35)")),
                      Taxane_response=anno_barplot(pheno$Taxane_response,ylim =c(-4,1),border =F,show_name = FALSE,bar_width = 0.4),
                      HER2DX=anno_barplot(pheno$HER2DX_pCR_likelihood_score,ylim =c(0,12),border =F,show_name = FALSE,bar_width = 0.4),
                      PIK3CA=anno_barplot(pheno$pik3ca_sig,ylim =c(0.5,2),border =F,show_name = FALSE,bar_width = 0.4),
                      Arm=pheno$Arm,
                      Response=pheno$Response,
                      Recurrence=pheno$RFS.status,
                      chr17q12_fusion=pheno$chr17q12_fusion,
                      ERBB2_fusion=pheno$ERBB2_fusion,
                      col=list(Arm=c("DHP"="#8491B4FF","T-DM1"="#91D1C2FF"),
                               Response=c("RD"="#525252","pCR"="#f0f0f0"),
                               Recurrence=c("1"="#525252","0"="#f0f0f0"),
                               chr17q12_fusion=c("1"="#525252","0"="#f0f0f0"),
                               ERBB2_fusion=c("1"="#525252","0"="#f0f0f0")),
                      na_col="grey",show_annotation_name = FALSE,
                      border=TRUE)


pheno=TME_cluster3
ha3=HeatmapAnnotation(foo=anno_block(gp= gpar(fill = "#fb9a99"),labels = c("HER2-enriched (n = 104)")),
                      Taxane_response=anno_barplot(pheno$Taxane_response,ylim =c(-4,1),border =F,show_name = FALSE,bar_width = 0.4),
                      HER2DX=anno_barplot(pheno$HER2DX_pCR_likelihood_score,ylim =c(0,12),border =F,show_name = FALSE,bar_width = 0.4),
                      PIK3CA=anno_barplot(pheno$pik3ca_sig,ylim =c(0.5,2),border =F,show_name = FALSE,bar_width = 0.4),
                      Arm=pheno$Arm,
                      Response=pheno$Response,
                      Recurrence=pheno$RFS.status,
                      chr17q12_fusion=pheno$chr17q12_fusion,
                      ERBB2_fusion=pheno$ERBB2_fusion,
                      col=list(Arm=c("DHP"="#8491B4FF","T-DM1"="#91D1C2FF"),
                               Response=c("RD"="#525252","pCR"="#f0f0f0"),
                               Recurrence=c("1"="#525252","0"="#f0f0f0"),
                               chr17q12_fusion=c("1"="#525252","0"="#f0f0f0"),
                               ERBB2_fusion=c("1"="#525252","0"="#f0f0f0")),
                      na_col="grey",show_annotation_name = FALSE,
                      border=TRUE)
table(mydata$TME_cluster)

pheno=TME_cluster4
ha4=HeatmapAnnotation(foo=anno_block(gp= gpar(fill = "#e31a1c"),labels = c("Basal-like (n = 10)")),
                      Taxane_response=anno_barplot(pheno$Taxane_response,ylim =c(-4,1),border =F,show_name = FALSE,bar_width = 0.4),
                      HER2DX=anno_barplot(pheno$HER2DX_pCR_likelihood_score,ylim =c(0,12),border =F,show_name = FALSE,bar_width = 0.4),
                      PIK3CA=anno_barplot(pheno$pik3ca_sig,ylim =c(0.5,2),border =F,show_name = FALSE,bar_width = 0.4),
                      Arm=pheno$Arm,
                      Response=pheno$Response,
                      Recurrence=pheno$RFS.status,
                      chr17q12_fusion=pheno$chr17q12_fusion,
                      ERBB2_fusion=pheno$ERBB2_fusion,
                      col=list(Arm=c("DHP"="#8491B4FF","T-DM1"="#91D1C2FF"),
                               Response=c("RD"="#525252","pCR"="#f0f0f0"),
                               Recurrence=c("1"="#525252","0"="#f0f0f0"),
                               chr17q12_fusion=c("1"="#525252","0"="#f0f0f0"),
                               ERBB2_fusion=c("1"="#525252","0"="#f0f0f0")),
                      na_col="grey",show_annotation_name =T,
                      border=TRUE)
col_fun = colorRamp2(c(-2, 0, 2), c("#377EB8", "white", "#E41A1C"))
ht_list =Heatmap(TME_function, name = "signature", show_row_names = FALSE, width = unit(10, "mm"),
                 col = structure(names = c("Cellular_process","Metabolism"),c("#A6761D","#666666")),
                                 #brewer.pal(6,"Dark2")[]),
                 row_split = factor(TME_function, levels = c("Cellular_process","Metabolism")))+ 
  Heatmap(mat1, col = col_fun, name = "Scaled Signature score",
          clustering_distance_columns = "spearman",
          show_row_dend = FALSE, show_column_dend = T,
          show_column_names = FALSE,
          show_row_names = FALSE,
          top_annotation = ha1,
          row_split = factor(TME_function, levels = c("Cellular_process","Metabolism")), 
          row_title_gp = gpar(col = "#FFFFFF00"), width = unit(3, "cm"))+
  Heatmap(mat2,col=col_fun, show_column_names = FALSE, show_row_names = FALSE,
          show_column_dend = T,top_annotation=ha2,
          show_heatmap_legend = FALSE, width = unit(3, "cm"))+
  Heatmap(mat3,col=col_fun, show_column_names = FALSE, show_row_names = FALSE,
          show_column_dend = T,top_annotation=ha3,
          show_heatmap_legend = FALSE, width = unit(9, "cm"))+
  Heatmap(mat4,col=col_fun, show_column_names = FALSE, 
          show_column_dend = T,top_annotation=ha4,
          show_heatmap_legend = FALSE, width = unit(1, "cm"))

p0=draw(ht_list,annotation_legend_side = "left", heatmap_legend_side = "left") 
#12X6

# figS6i
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
rna=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.txt')
data=inner_join(rna,clin,by="patientID")%>%as.data.frame()
variable=c("Taxane_response","HER2DX_pCR_likelihood_score","pik3ca_sig","ABC_transporter","Exosome","Hypoxia")
genomic=slice_metric(data,variable,20,80,5)
## re-run the batch logistic
selected_columns <- grep("_per_", colnames(genomic), value = TRUE)
Interact_result=Logistic_batch_adjER(genomic,"pCR","Arm",selected_columns ,"ER")%>%as.data.frame()
Interact_result=separate(Interact_result,biomarker, into = c("group", "per"), sep = "_per_", remove = FALSE)
Interact_result$per=as.numeric(Interact_result$per)
Interact_result$interaction_coefficient=as.numeric(Interact_result$interaction_coefficient)
Interact_result$interaction_coefficient_LCI=as.numeric(Interact_result$interaction_coefficient_LCI)
Interact_result$interaction_coefficient_UCI=as.numeric(Interact_result$interaction_coefficient_UCI)
# ABC_transporter_per_70, Exosome_per_45, Hypoxia_per_80, pik3ca_sig_per_80

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
  scale_color_brewer(palette="Dark2") +
  theme_minimal() +
  labs(title = "",
       subtitle = "(mean +/- standard error)",
       x = "", 
       y = "",
       color = "Genomic Biomarker")
#12X3.5

# figS6j
library(data.table)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
rna=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.txt')
df=inner_join(rna,clin,by="patientID")%>%as.data.frame()
res<- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = df[df$sspbc.subtype=="LumA",])
ShowRegTable(res)
res<- glm(as.numeric(pCR) ~ Arm, family = "binomial", data = df[df$sspbc.subtype=="LumB",])
ShowRegTable(res)
res<- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = df[df$sspbc.subtype=="Her2",])
ShowRegTable(res)
res<- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = df[df$sspbc.subtype=="Basal",])
ShowRegTable(res)

res<- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = df[df$sspbc.subtype=="LumA",])
ShowRegTable(res)

interaction_1<- glm(as.numeric(pCR) ~ sspbc.subtype+Arm+ER, family = "binomial", data = df)
interaction_2<- glm(as.numeric(pCR) ~ sspbc.subtype*Arm+sspbc.subtype+Arm+ER, family = "binomial", data =df)
interaction_lr <- lrtest(interaction_1,interaction_2)

# her2dx
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
rna=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.txt')
data=inner_join(rna,clin,by="patientID")%>%as.data.frame()
variable=c("Taxane_response","HER2DX_pCR_likelihood_score","pik3ca_sig","ABC_transporter","Exosome","Hypoxia")
genomic=slice_metric(data,variable,20,80,5)
## re-run the batch logistic
selected_columns <- grep("_per_", colnames(genomic), value = TRUE)
Interact_result=Logistic_batch_adjER(genomic,"pCR","Arm",selected_columns ,"ER")%>%as.data.frame()
Interact_result=separate(Interact_result,biomarker, into = c("group", "per"), sep = "_per_", remove = FALSE)
Interact_result$per=as.numeric(Interact_result$per)
Interact_result$interaction_coefficient=as.numeric(Interact_result$interaction_coefficient)
Interact_result$interaction_coefficient_LCI=as.numeric(Interact_result$interaction_coefficient_LCI)
Interact_result$interaction_coefficient_UCI=as.numeric(Interact_result$interaction_coefficient_UCI)
#HER2DX_pCR_likelihood_score_per_40
genomic$Arm=factor(genomic$Arm,levels = c("DHP","T-DM1"))
table(genomic$HER2DX_pCR_likelihood_score_per_40,genomic$Arm)
res=Logistic_batch_continuous_subgroup(genomic,"HER2DX_pCR_likelihood_score_per_40")%>%as.data.frame()
library(forestploter)
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS6/Subgroup_Forestplot.csv")
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

#7X5 65%

