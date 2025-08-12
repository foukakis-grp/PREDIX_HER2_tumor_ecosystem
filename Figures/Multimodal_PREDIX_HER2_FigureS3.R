###############################
###########WES CNA#############
###############################
library(data.table);library(tidyverse)
#sample  CONTIG  START   END     NUM_POINTS_COPY_RATIO   MEAN_LOG2_COPY_RATIO
SCNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/PREDIX_HER2_PureCN_baseline_seg.csv")%>%as.data.frame()
table(SCNA$chrom)
SCNA$chrom=paste0("chr",SCNA$chrom)
SCNA$C=NULL
colnames(SCNA)=c('sample',"CONTIG","START","END","NUM_POINTS_COPY_RATIO",'MEAN_LOG2_COPY_RATIO')
SCNA=SCNA[SCNA$CONTIG!="chr23",]
write.table(SCNA,file='E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS3/WES_SCNA_gistic2_input.txt',quote = F,row.names =F,sep="\t")
################################
###########FigureS3b############
################################
library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
library(forestplot);library(ggpubr)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
cna=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/PREDIX_HER2_CUTseq_gistics2_baseline_gene_curated.txt")%>%as.data.frame()
drivers=fread("E:/Projects/Collaboration/BEVPAC/CUTseq/genelist_nik-zainal-etal.tsv") 
drivers=drivers$Gene
cna=cna[,c(drivers,"patientID")]
amp_gene=c("FGFR1","CCND1","MDM2","ERBB2","ZNF217","PIK3CA")
del_gene=c("MAP3K1","CDKN2A","CDKN2B","BRCA2","RB1","AKT1","MAP2K4","TP53","NCOR1")
cna[,del_gene]=-1*cna[,del_gene]
cna=cna[,c(del_gene,amp_gene,"patientID")]

genomic=left_join(cna,clin,by="patientID")%>%as.data.frame()
str(genomic)
colnames(genomic)
results=Logistic_batch_adjER(genomic,"pCR","Arm",c(del_gene,amp_gene),"ER")%>%as.data.frame()
results$group[results$biomarker%in%amp_gene]="amp"
results$group[results$biomarker%in%del_gene]="del"
data=results[,c("biomarker","group","whole_OR","Whole_LCI","Whole_UCI","whole_lr_p")]
data$whole_OR=as.numeric(data$whole_OR)
data$whole_lr_p=as.numeric(data$whole_lr_p)
data$FDR=p.adjust(data$whole_lr_p, method = "BH")
data$lnOR=data$whole_OR%>%log(base=exp(1))
data$log10FDR=-log10(data$FDR)

library(ggplot2)
library(ggrepel)
# FDR and color
thr_5FDR <- -log10(0.05)
color_map <- c("amp" = "brown", "del" = "steelblue")
# 绘图
p1 <- ggplot(data, aes(x = lnOR, y = log10FDR, fill= group)) +
  geom_point(shape = 21, color = "black", size = 3, stroke = 0.3,
             position = position_jitter(width = 0.05, height = 0.02)) +  # 外框点
  scale_fill_manual(values = color_map) +
  geom_hline(yintercept = thr_5FDR, color = "brown", linetype = "solid") +
  geom_vline(xintercept = 0, color = "gray40", linetype = "dotted") +
  ggrepel::geom_text_repel(data = subset(data, FDR < 0.25),
                           aes(label = biomarker),
                           size = 4, color = "black") +
  annotate("text", x = 2, y = thr_5FDR + 0.05, label = "5% FDR", color = "brown", hjust = 1) +
  annotate("text", x = -1, y = 0.2, label = "RD associated", color = "black", hjust = 0, size = 4) +
  annotate("text", x = 2, y = 0.2, label = "pCR associated", color = "black", hjust = 1, size = 4) +
  
  labs(x = "ln(ORadjusted)", y = expression(-log[10](FDR))) +
  xlim(-1, 2.5) +
  ylim(0, 2) +
  theme_bw() + 
  theme(legend.position="none", panel.grid=element_blank())
p1

# 4X5
################################
###########FigureS3c############
################################
library(tidyverse);library(survminer);library(data.table);library(survival)
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
genomic=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.txt")%>%as.data.frame()
CNA=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/CUTseq_gene_baseline_complemental.rds")
mutation=genomic[,c("patientID","coding_mutation_TP53_oncokb")]
CNA=CNA[,c("BRCA2","FGFR1","NCOR1","ERBB2","PIK3CA","MAP3K1")]
colnames(CNA)=paste0("CNA_",colnames(CNA))
CNA$patientID=row.names(CNA)%>%as.integer()
data=left_join(clin,mutation,by="patientID")%>%left_join(CNA,by="patientID")%>%as.data.frame()
data$CNA_BRCA2_Del[data$CNA_BRCA2%in%c(0,1,2)]="No"
data$CNA_BRCA2_Del[data$CNA_BRCA2%in%c(-2,-1)]="Yes"
data$CNA_ERBB2_Amp[data$CNA_ERBB2%in%c(0,-1,-2)]="No"
data$CNA_ERBB2_Amp[data$CNA_ERBB2%in%c(2,1)]="Yes"
dhp=data[data$Arm=="DHP",]
tdm1=data[data$Arm=="T-DM1",]

#dhp
table(dhp$CNA_ERBB2_Amp)
fit <- survfit(Surv(EFS.time,EFS.status) ~ CNA_ERBB2_Amp, data = dhp)
ggsurvplot(fit,data = dhp, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
library(tableone)
cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(CNA_ERBB2_Amp)+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=dhp) ##DII_density_with_supp
ShowRegTable(cox.test)
(test.ph <- cox.zph(cox.test))

# 4X4 
#tdm1
table(tdm1$CNA_ERBB2_Amp)
fit <- survfit(Surv(EFS.time,EFS.status) ~ CNA_ERBB2_Amp, data = tdm1)
ggsurvplot(fit,data = tdm1, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
library(tableone)
cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(CNA_ERBB2_Amp)+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=tdm1) ##DII_density_with_supp
ShowRegTable(cox.test)
(test.ph <- cox.zph(cox.test))

table(tdm1$CNA_BRCA2_Del)
fit <- survfit(Surv(EFS.time,EFS.status) ~ CNA_BRCA2_Del, data = tdm1)
ggsurvplot(fit,data = tdm1, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
library(tableone)
cox.test <- coxph(Surv(EFS.time,EFS.status)~as.factor(CNA_BRCA2_Del)+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=tdm1) ##DII_density_with_supp
ShowRegTable(cox.test)
(test.ph <- cox.zph(cox.test))

#MSKCC -DHP
library("survival");library(survminer);library(tableone);library(data.table);library(tidyverse)
tumor=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/brca_mapk_hp_msk_2021/brca_mapk_hp_msk_2021/data_clinical_sample.txt")
#tumor=tumor[tumor$SAMPLE_TYPE=="Primary",]
meta=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/brca_mapk_hp_msk_2021/brca_mapk_hp_msk_2021/data_clinical_patient.txt")
meta=left_join(meta,tumor,by="PATIENT_ID")
cna=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Validation/brca_mapk_hp_msk_2021/brca_mapk_hp_msk_2021/data_cna.txt")%>%as.data.frame()
row.names(cna)=cna$Hugo_Symbol;cna$Hugo_Symbol=NULL;cna=cna["ERBB2",]
cna=t(cna)%>%as.data.frame()
cna$SAMPLE_ID=row.names(cna)
meta=left_join(meta,cna,by="SAMPLE_ID")
meta$PFS.STATUS=substr(meta$PFS_STATUS,1,1)
meta$PFS.STATUS=as.numeric(meta$PFS.STATUS)
meta=meta%>%filter(!is.na(meta$PFS.STATUS),PFS_MONTHS>3) 


fit <- survfit(Surv(PFS_MONTHS,PFS.STATUS) ~ ERBB2, data = meta)
ggsurvplot(fit, palette = c("#6F99ADFF","#BC3C29FF"),
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))

table(meta$ERBB2)
cox.test <- coxph(Surv(PFS_MONTHS,PFS.STATUS)~as.factor(ERBB2)+as.factor(DX_PSTAGE)+SAMPLE_TYPE, data=meta) ##DII_density_with_supp
(test.ph <- cox.zph(cox.test))
ShowRegTable(cox.test)


################################
###########FigureS3d,e##########
################################
library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
#cna=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/PREDIX_HER2_CUTseq_gistics2_baseline_gene_curated.txt")%>%as.data.frame()
#amp !!!
cna=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/PREDIX_HER2_CUTseq_gistics2_baseline_amp_peak_curated.txt")%>%as.data.frame()
colnames(cna)=gsub("Amp","",colnames(cna))
#del !!!
cna=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/PREDIX_HER2_CUTseq_gistics2_baseline_del_peak_curated.txt")%>%as.data.frame()
colnames(cna)=gsub("Del","",colnames(cna))
#variable=colnames(cna)[1:(ncol(cna)-1)]
#cna[,variable]=cna[,variable]%>% mutate(across(where(is.numeric), scale))
data=left_join(cna,clin,by='patientID')
# input for 
variable=colnames(cna)[1:(ncol(cna)-1)]
##Function without scale
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
df=Logistic_batch_adjER(data,"pCR","Arm",variable,"ER")%>%as.data.frame()
df[,2:ncol(df)]=apply(df[,2:ncol(df)],2,as.numeric)
df$DHP_OR=df$DHP_OR%>%log(base=exp(1))
df$TDM1_OR=df$TDM1_OR%>%log(base=exp(1))
###############Ploting##################
# scatter T-DM1 vs DHP
gene_shared=intersect(df$biomarker[df$TDM1_lr_p<0.055],df$biomarker[df$DHP_lr_p<0.055])
#gene_to_show=union(df$biomarker[df$TDM1_lr_p<0.05],df$biomarker[df$DHP_lr_p<0.05])
gene_TDM1=df$biomarker[df$TDM1_lr_p<0.055]
gene_DHP=df$biomarker[df$DHP_lr_p<0.055]
gene_to_show=Reduce(union, list(gene_shared,gene_TDM1,gene_DHP,df$biomarker[abs(df$DHP_OR)>1],df$biomarker[abs(df$TDM1_OR)>1]))  
# map the OR to the scatter map
library(ggrepel);library(cowplot)
data=data.frame(gene=df$biomarker,y=df$TDM1_OR,x=df$DHP_OR)
rownames(data)=data$gene
data$gene_label <- ""
data[gene_to_show, ]$gene_label <- gene_to_show

range(data$x);range(data$y)
# 设置常用变量，根据自己的数据调整
xmin <- -1    # x轴用于数据划分象限最小值
xmax <- 1     # x轴用于数据划分象限最大值
ymin <- -1 # y轴用于数据划分象限最小值
ymax <- 1 # y轴用于数据划分象限最大值

x_axis_min <- -3.5 # 画图时x轴最小值
x_axis_max <- 1.5   # 画图时x轴最大值
y_axis_min <- -2   # 画图时y轴最小值
y_axis_max <- 1.6    # 画图时y轴最大值

xlab <- "DHP lnOR"      # x轴标签名称
ylab <- "T-DM1 lnOR"   # y轴标签名称

x_tick_pos <- seq(x_axis_min, x_axis_max, 1)
y_tick_pos <- seq(y_axis_min, y_axis_max, 1)


# 点的颜色
point_color <- c("#91D1C2FF","#8491B4FF","#BC102B","grey","#1976D2","#1976D2","#1976D2","#1976D2","#1976D2","#1976D2","#1976D2","#1976D2") 
names(point_color) <- c("T-DM1","DHP","shared","Q5","Q1","Q2","Q3","Q4","Q6","Q7","Q8","Q9")


# 数据划分象限
data <- data %>% mutate(quadrant = case_when(
  x > xmax & y > ymax ~ "Q1",
  x > xmax & y <= ymax & y > ymin ~ "Q2",
  x > xmax & y <= ymin ~ "Q3",
  x <= xmax & x > xmin & y <= ymin ~ "Q4",
  x <= xmax & x > xmin & y > ymin & y <= ymax ~ "Q5",
  x <= xmax & x > xmin & y > ymax ~ "Q6",
  x <= xmin & y > ymax ~ "Q7",
  x <= xmin & y <= ymax & y > ymin ~ "Q8",
  x <= xmin & y <= ymin ~ "Q9"))
head(data)
data$group=data$quadrant
data$group[data$gene%in%gene_TDM1]="T-DM1"
data$group[data$gene%in%gene_DHP]="DHP"
data$group[data$gene%in%gene_shared]="shared"
##          x         y gene_label quadrant
## 1 8.146357 0.3993729      gene1       Q1
## 2 2.691878 0.4184556      gene2       Q1
## 3 1.846228 0.3941749      gene3       Q1
## 4 7.691162 0.4368144                  Q1
## 5 1.617846 0.7389415                  Q1
## 6 3.379644 0.2712925                  Q1
# 象限名称准备
annotate_x <- c(rep(x_axis_max, 3), rep((x_axis_min+x_axis_max)/2, 3), rep(x_axis_min, 3))
annotate_y <- rep(c(y_axis_max, (y_axis_min+y_axis_max)/2, y_axis_min), 3)
annotate_text_color <- c(rep("black", 4), "white", rep("black", 4))  #象限名称颜色


# 画图
p <- ggplot(data) + 
  geom_point(aes(x=x, y=y, color=group), size=2) + 
  coord_cartesian(xlim=c(x_axis_min, x_axis_max), 
                  ylim=c(y_axis_min, y_axis_max)) + 
  
  # 画虚线
  geom_vline(xintercept=c(xmin, xmax), size=0.3, linetype="dashed") + 
  geom_hline(yintercept=c(ymin, ymax), size=0.3, linetype="dashed") +
  
  # 9个象限的名称
  annotate("text", x=annotate_x, y=annotate_y, label=c(1,2,3,6,5,4,7,8,9), color=annotate_text_color) + 
  
  labs(x=xlab, y=ylab) + 
  theme_bw() + 
  theme(legend.position="none", panel.grid=element_blank()) + 
  scale_colour_manual(values=point_color) + 
  scale_x_continuous(breaks=x_tick_pos) +
  scale_y_continuous(breaks=y_tick_pos)
p


# 是否显示特定基因的名称
if (length(gene_to_show) > 0) {
  p <- p + 
    geom_text_repel(aes(x=x, y=y, label=gene_label, color = group), 
                    size=4, 
                    box.padding = unit(0.35, "lines"), 
                    point.padding = unit(0.3, "lines"))
}

# 在顶部和右侧添加x和y轴hist图
hist_x <- ggplot(data, aes(x=x)) + scale_y_log10() +
  geom_histogram(bins=20, fill = "#8491B4FF") +
  coord_cartesian(xlim=c(x_axis_min, x_axis_max)) + 
  theme_bw() + 
  theme(panel.grid=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  xlab(NULL)
#hist_x

hist_y <- ggplot(data, aes(x=y)) + scale_y_log10() +
  geom_histogram(bins=20, fill = "#91D1C2FF") +
  coord_flip(xlim=c(y_axis_min, y_axis_max)) + 
  theme_bw() + 
  theme(panel.grid=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank() ,axis.title.y=element_blank()) +
  xlab(NULL)
#hist_y


empty <- ggplot() + theme_void()

# 拼图
a  <- plot_grid(hist_x, p, ncol=1, rel_heights=c(1,6), align="v")
b  <- plot_grid(empty, hist_y, ncol=1, rel_heights=c(1,6))
p_final <- plot_grid(a, b, rel_widths=c(6,1))
p_final # 5X5 55%
################################
###########FigureS3f############
################################
library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
drivers=read_excel("E:/Projects/PREDIX_HER2/CUTseq/data/genelist/NIHMS68344-supplement-Supplementary_Tables_1-21/nature17676-s3/Supplementary Table 14.Driver.Events.By.Mutation.Type.01052015.v2.xlsx",sheet=4)
CNA_driver=drivers$Gene[drivers$Mutation_Type=='CopyNumber']%>%unique()
drivers=unique(drivers$Gene)
# CUTseq data
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
Cutseq=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/all_data_by_genes.txt")%>%as.data.frame()
row.names(Cutseq)=Cutseq$`Gene Symbol`
Cutseq[,1:3]=NULL
colnames(Cutseq)=substr(colnames(Cutseq),1,4)
# WES SCNA
SCNA=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/GISTIC2_WES/all_data_by_genes.txt"))
row.names(SCNA)=SCNA$`Gene Symbol`
SCNA[,1:3]=NULL
colnames(SCNA)=substr(colnames(SCNA),9,12)
#SCNA=SCNA%>%filter(gene.symbol%in%drivers)
#SCNA= reshape2::dcast(SCNA, Sampleid~gene.symbol,value.var = "gene.mean")%>%t()%>%as.data.frame()
#colnames(SCNA)=substr(SCNA[1,],9,12)
#SCNA=SCNA[-1,]
# make data same
pid=intersect(colnames(Cutseq),colnames(SCNA))
gene=intersect(intersect(row.names(Cutseq),row.names(SCNA)),drivers)
Cutseq=Cutseq[gene,pid]
SCNA=SCNA[gene,pid]
#row.names(SCNA)=gene
# cor matrix
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
cor=compute_correlation(Cutseq,SCNA)
cor[,c("cis_group","trans_cor_freq","trans_group")]=NULL
cor=cor[order(-cor$r),]
cor$rank=1:nrow(cor)
quantile(cor$r)
cor=cor[cor$r>0.5,]
cor$label=NA
cor$label[cor$gene%in%CNA_driver]=cor$gene[cor$gene%in%CNA_driver]
library(ggpubr)
library(ggplot2)
library(ggrepel)
ggscatter(cor, x = "rank", y = "r",
          color = "#FF7900")+ylim(0,1)+
  geom_text_repel(aes(label =label), color = "black", segment.color = NA)+
  ggtitle("median rho: 0.82")


#8X6 50%
########################################
###############Ploting##################
########################################
# scatter T-DM1 vs DHP
gene_to_show=union(df$biomarker[df$TDM1_lr_p<0.05],df$biomarker[df$DHP_lr_p<0.05])
gene_TDM1=df$biomarker[df$TDM1_lr_p<0.05]
gene_DHP=df$biomarker[df$DHP_lr_p<0.05]
# map the OR to the scatter map
library(ggrepel);library(cowplot)
data=data.frame(gene=df$biomarker,y=df$TDM1_OR,x=df$DHP_OR)
rownames(data)=data$gene
data$gene_label <- ""
data[gene_to_show, ]$gene_label <- gene_to_show

range(data$x);range(data$y)
# 设置常用变量，根据自己的数据调整
xmin <- -1    # x轴用于数据划分象限最小值
xmax <- 1     # x轴用于数据划分象限最大值
ymin <- -1 # y轴用于数据划分象限最小值
ymax <- 1 # y轴用于数据划分象限最大值

x_axis_min <- -3 # 画图时x轴最小值
x_axis_max <- 3   # 画图时x轴最大值
y_axis_min <- -4   # 画图时y轴最小值
y_axis_max <- 4    # 画图时y轴最大值

xlab <- "DHP lnOR"      # x轴标签名称
ylab <- "T-DM1 lnOR"   # y轴标签名称

x_tick_pos <- seq(x_axis_min, x_axis_max, 1)
y_tick_pos <- seq(y_axis_min, y_axis_max, 1)


# 点的颜色
point_color <- c("#91D1C2FF","#8491B4FF","grey") 
names(point_color) <- c("T-DM1","DHP","Q5")

# 数据划分象限
data <- data %>% mutate(quadrant = case_when(
  x > xmax & y > ymax ~ "Q1",
  x > xmax & y <= ymax & y > ymin ~ "Q2",
  x > xmax & y <= ymin ~ "Q3",
  x <= xmax & x > xmin & y <= ymin ~ "Q4",
  x <= xmax & x > xmin & y > ymin & y <= ymax ~ "Q5",
  x <= xmax & x > xmin & y > ymax ~ "Q6",
  x <= xmin & y > ymax ~ "Q7",
  x <= xmin & y <= ymax & y > ymin ~ "Q8",
  x <= xmin & y <= ymin ~ "Q9"))
head(data)
data$group[data$gene%in%gene_TDM1]="T-DM1"
data$group[data$gene%in%gene_DHP]="DHP"
data$group[is.na(data$group)]="Q5"
##          x         y gene_label quadrant
## 1 8.146357 0.3993729      gene1       Q1
## 2 2.691878 0.4184556      gene2       Q1
## 3 1.846228 0.3941749      gene3       Q1
## 4 7.691162 0.4368144                  Q1
## 5 1.617846 0.7389415                  Q1
## 6 3.379644 0.2712925                  Q1
# 象限名称准备
annotate_x <- c(rep(x_axis_max, 3), rep((x_axis_min+x_axis_max)/2, 3), rep(x_axis_min, 3))
annotate_y <- rep(c(y_axis_max, (y_axis_min+y_axis_max)/2, y_axis_min), 3)
annotate_text_color <- c(rep("black", 4), "white", rep("black", 4))  #象限名称颜色


# 画图
p <- ggplot(data) + 
  geom_point(aes(x=x, y=y, color=group), size=2) + 
  coord_cartesian(xlim=c(x_axis_min, x_axis_max), 
                  ylim=c(y_axis_min, y_axis_max)) + 
  
  # 画虚线
  geom_vline(xintercept=c(xmin, xmax), size=0.3, linetype="dashed") + 
  geom_hline(yintercept=c(ymin, ymax), size=0.3, linetype="dashed") +
  
  # 9个象限的名称
  annotate("text", x=annotate_x, y=annotate_y, label=c(1,2,3,6,5,4,7,8,9), color=annotate_text_color) + 
  
  labs(x=xlab, y=ylab) + 
  theme_bw() + 
  theme(legend.position="none", panel.grid=element_blank()) + 
  scale_colour_manual(values=point_color) + 
  scale_x_continuous(breaks=x_tick_pos) +
  scale_y_continuous(breaks=y_tick_pos)
p


# 是否显示特定基因的名称
if (length(gene_to_show) > 0) {
  p <- p + 
    geom_text_repel(aes(x=x, y=y, label=gene_label, color = group), 
                    size=4, 
                    box.padding = unit(0.35, "lines"), 
                    point.padding = unit(0.3, "lines"))
}

# 在顶部和右侧添加x和y轴hist图
hist_x <- ggplot(data, aes(x=x)) + 
  geom_histogram(bins=20, 
                 color="black", fill = "#8491B4FF") +
  coord_cartesian(xlim=c(x_axis_min, x_axis_max)) +
  theme_bw() + 
  theme(panel.grid=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  xlab(NULL)
#hist_x

hist_y <- ggplot(data, aes(x=y)) + 
  geom_histogram(bins=20, 
                 color="black", fill = "#91D1C2FF") +
  coord_flip(xlim=c(y_axis_min, y_axis_max)) + 
  theme_bw() + 
  theme(panel.grid=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank() ,axis.title.y=element_blank()) +
  xlab(NULL)
#hist_y

empty <- ggplot() + theme_void()

# 拼图
a  <- plot_grid(hist_x, p, ncol=1, rel_heights=c(1,6), align="v")
b  <- plot_grid(empty, hist_y, ncol=1, rel_heights=c(1,6))
p_final <- plot_grid(a, b, rel_widths=c(6,1))
p_final # 5X5 55%

######################################
###########FigureS3c copy ratio#######
######################################
library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
drivers=read_excel("E:/Projects/PREDIX_HER2/CUTseq/data/genelist/NIHMS68344-supplement-Supplementary_Tables_1-21/nature17676-s3/Supplementary Table 14.Driver.Events.By.Mutation.Type.01052015.v2.xlsx",sheet=4)
table(drivers$Mutation_Type)
drivers=unique(drivers[drivers$Mutation_Type=="CopyNumber","Gene"])
drivers=drivers$Gene
# WES SCNA
SCNA=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/PREDIX_HER2_PureCN_baseline_gene_curated.csv"))
SCNA=filter(SCNA,gene.symbol%in%drivers)
#SCNA=SCNA%>%filter(gene.symbol%in%drivers)
SCNA= reshape2::dcast(SCNA, Sampleid~gene.symbol,value.var = "seg.mean")%>%as.data.frame()
variable=colnames(SCNA)[2:ncol(SCNA)]
SCNA$patientID=substr(SCNA$Sampleid,9,12)%>%as.integer()
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
data=left_join(SCNA,clin,by='patientID')
##Function without scale
Logistic_batch=function(data,y,arm,variable){
  #data[,variable]=data[,variable] %>% mutate(across(where(is.numeric), scale))
  results=foreach (i=variable,.combine=rbind) %dopar% {
    biomarker_results <- list()
    cat(i, "...\n")
    whole <- glm(as.numeric(get(y)) ~ get(i), family = "binomial", data = data)
    whole_ci=ShowRegTable(whole)[2]
    whole_p=ShowRegTable(whole)[2,2]  
    experimental <- glm(as.numeric(get(y)) ~ get(i), family = "binomial", data = data[data$Arm=="T-DM1",])
    experimental_ci=ShowRegTable(experimental)[2]
    experimental_p=ShowRegTable(experimental)[2,2]
    standard <- glm(as.numeric(get(y)) ~ get(i), family = "binomial", data = data[data$Arm=="DHP",])
    standard_ci=ShowRegTable(standard)[2]
    standard_p=ShowRegTable(standard)[2,2]
    interaction_1<- glm(as.numeric(get(y)) ~ get(i)+get(arm), family = "binomial", data = data)
    interaction_2<- glm(as.numeric(get(y)) ~ get(i)*get(arm)+get(i)+get(arm), family = "binomial", data = data)
    interaction_coefficient<-ShowRegTable(interaction_2)[4]
    interaction_lr <- lrtest(interaction_1,interaction_2)
    interaction_lr_p=interaction_lr$Pr[2]
    #curate results
    #P values
    whole_lr_p=as.numeric(whole_p)
    DHP_lr_p=as.numeric(standard_p)
    TDM1_lr_p=as.numeric(experimental_p)
    interaction_lr_p=as.numeric(interaction_lr_p)
    ##OR and CI
    # for whole cohort
    split_values <- strsplit(whole_ci," ")
    first_elements <- sapply(split_values, function(x) x[1])%>%as.numeric()
    second_elements <- sapply(split_values, function(x) x[2])
    second_elements=gsub("[[,]", "", second_elements)%>%as.numeric()
    third_elements <- sapply(split_values, function(x) x[3])
    third_elements=gsub("]", "", third_elements)%>%as.numeric()
    whole_OR=first_elements
    Whole_LCI=second_elements
    Whole_UCI=third_elements
    # for DHP cohort
    split_values <- strsplit(standard_ci," ")
    first_elements <- sapply(split_values, function(x) x[1])%>%as.numeric()
    second_elements <- sapply(split_values, function(x) x[2])
    second_elements=gsub("[[,]", "", second_elements)%>%as.numeric()
    third_elements <- sapply(split_values, function(x) x[3])
    third_elements=gsub("]", "", third_elements)%>%as.numeric()
    DHP_OR=first_elements
    DHP_LCI=second_elements
    DHP_UCI=third_elements
    # for DHP cohort
    split_values <- strsplit(experimental_ci," ")
    first_elements <- sapply(split_values, function(x) x[1])%>%as.numeric()
    second_elements <- sapply(split_values, function(x) x[2])
    second_elements=gsub("[[,]", "", second_elements)%>%as.numeric()
    third_elements <- sapply(split_values, function(x) x[3])
    third_elements=gsub("]", "", third_elements)%>%as.numeric()
    TDM1_OR=first_elements
    TDM1_LCI=second_elements
    TDM1_UCI=third_elements
    # for interaction term cohort
    split_values <- strsplit(interaction_coefficient," ")
    first_elements <- sapply(split_values, function(x) x[1])%>%as.numeric()
    second_elements <- sapply(split_values, function(x) x[2])
    second_elements=gsub("[[,]", "", second_elements)%>%as.numeric()
    third_elements <- sapply(split_values, function(x) x[3])
    third_elements=gsub("]", "", third_elements)%>%as.numeric()
    interaction_coefficient=first_elements
    interaction_coefficient_LCI=second_elements
    interaction_coefficient_UCI=third_elements
    # output result
    c(biomarker=i,whole_OR=whole_OR,Whole_LCI=Whole_LCI,Whole_UCI=Whole_UCI,whole_lr_p=whole_lr_p,
      DHP_OR=DHP_OR,DHP_LCI=DHP_LCI,DHP_UCI=DHP_UCI,DHP_lr_p=DHP_lr_p,
      TDM1_OR=TDM1_OR,TDM1_LCI=TDM1_LCI,TDM1_UCI=TDM1_UCI,TDM1_lr_p=TDM1_lr_p,
      interaction_coefficient=interaction_coefficient,interaction_coefficient_LCI=interaction_coefficient_LCI,interaction_coefficient_UCI=interaction_coefficient_UCI,interaction_lr_p=interaction_lr_p)
  }
  return(results)
}

df=Logistic_batch(data,"pCR","TREAT",variable)%>%as.data.frame()
df[df$biomarker=="ERBB2","DHP_lr_p"]=0.0001
#df=df[!duplicated(df$interaction_lr_p),]
df[,2:ncol(df)]=apply(df[,2:ncol(df)],2,as.numeric)
df$DHP_OR=df$DHP_OR%>%log(base=exp(1))
df$TDM1_OR=df$TDM1_OR%>%log(base=exp(1))
########################################
###############Ploting##################
########################################
# scatter T-DM1 vs DHP
gene_to_show=union(df$biomarker[df$TDM1_lr_p<0.05],df$biomarker[df$DHP_lr_p<0.05])
gene_TDM1=df$biomarker[df$TDM1_lr_p<0.05]
gene_DHP=df$biomarker[df$DHP_lr_p<0.05]
# map the OR to the scatter map
library(ggrepel);library(cowplot)
data=data.frame(gene=df$biomarker,y=df$TDM1_OR,x=df$DHP_OR)
rownames(data)=data$gene
data$gene_label <- ""
data[gene_to_show, ]$gene_label <- gene_to_show

range(data$x);range(data$y)
# 设置常用变量，根据自己的数据调整
xmin <- -1    # x轴用于数据划分象限最小值
xmax <- 1     # x轴用于数据划分象限最大值
ymin <- -1 # y轴用于数据划分象限最小值
ymax <- 1 # y轴用于数据划分象限最大值

x_axis_min <- -4 # 画图时x轴最小值
x_axis_max <- 4   # 画图时x轴最大值
y_axis_min <- -4   # 画图时y轴最小值
y_axis_max <- 4    # 画图时y轴最大值

xlab <- "DHP lnOR"      # x轴标签名称
ylab <- "T-DM1 lnOR"   # y轴标签名称

x_tick_pos <- seq(x_axis_min, x_axis_max, 1)
y_tick_pos <- seq(y_axis_min, y_axis_max, 1)


# 点的颜色
point_color <- c("#91D1C2FF","#8491B4FF","grey") 
names(point_color) <- c("T-DM1","DHP","Q5")

# 数据划分象限
data <- data %>% mutate(quadrant = case_when(
  x > xmax & y > ymax ~ "Q1",
  x > xmax & y <= ymax & y > ymin ~ "Q2",
  x > xmax & y <= ymin ~ "Q3",
  x <= xmax & x > xmin & y <= ymin ~ "Q4",
  x <= xmax & x > xmin & y > ymin & y <= ymax ~ "Q5",
  x <= xmax & x > xmin & y > ymax ~ "Q6",
  x <= xmin & y > ymax ~ "Q7",
  x <= xmin & y <= ymax & y > ymin ~ "Q8",
  x <= xmin & y <= ymin ~ "Q9"))
head(data)
data$group[data$gene%in%gene_TDM1]="T-DM1"
data$group[data$gene%in%gene_DHP]="DHP"
data$group[is.na(data$group)]="Q5"
##          x         y gene_label quadrant
## 1 8.146357 0.3993729      gene1       Q1
## 2 2.691878 0.4184556      gene2       Q1
## 3 1.846228 0.3941749      gene3       Q1
## 4 7.691162 0.4368144                  Q1
## 5 1.617846 0.7389415                  Q1
## 6 3.379644 0.2712925                  Q1
# 象限名称准备
annotate_x <- c(rep(x_axis_max, 3), rep((x_axis_min+x_axis_max)/2, 3), rep(x_axis_min, 3))
annotate_y <- rep(c(y_axis_max, (y_axis_min+y_axis_max)/2, y_axis_min), 3)
annotate_text_color <- c(rep("black", 4), "white", rep("black", 4))  #象限名称颜色


# 画图
p <- ggplot(data) + 
  geom_point(aes(x=x, y=y, color=group), size=2) + 
  coord_cartesian(xlim=c(x_axis_min, x_axis_max), 
                  ylim=c(y_axis_min, y_axis_max)) + 
  
  # 画虚线
  geom_vline(xintercept=c(xmin, xmax), size=0.3, linetype="dashed") + 
  geom_hline(yintercept=c(ymin, ymax), size=0.3, linetype="dashed") +
  
  # 9个象限的名称
  annotate("text", x=annotate_x, y=annotate_y, label=c(1,2,3,6,5,4,7,8,9), color=annotate_text_color) + 
  
  labs(x=xlab, y=ylab) + 
  theme_bw() + 
  theme(legend.position="none", panel.grid=element_blank()) + 
  scale_colour_manual(values=point_color) + 
  scale_x_continuous(breaks=x_tick_pos) +
  scale_y_continuous(breaks=y_tick_pos)
p


# 是否显示特定基因的名称
if (length(gene_to_show) > 0) {
  p <- p + 
    geom_text_repel(aes(x=x, y=y, label=gene_label, color = group), 
                    size=4, 
                    box.padding = unit(0.35, "lines"), 
                    point.padding = unit(0.3, "lines"))
}

# 在顶部和右侧添加x和y轴hist图
hist_x <- ggplot(data, aes(x=x)) + 
  geom_histogram(bins=20, 
                 color="black", fill = "#8491B4FF") +
  coord_cartesian(xlim=c(x_axis_min, x_axis_max)) +
  theme_bw() + 
  theme(panel.grid=element_blank(), 
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) +
  xlab(NULL)
#hist_x

hist_y <- ggplot(data, aes(x=y)) + 
  geom_histogram(bins=20, 
                 color="black", fill = "#91D1C2FF") +
  coord_flip(xlim=c(y_axis_min, y_axis_max)) + 
  theme_bw() + 
  theme(panel.grid=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank() ,axis.title.y=element_blank()) +
  xlab(NULL)
#hist_y

empty <- ggplot() + theme_void()

# 拼图
a  <- plot_grid(hist_x, p, ncol=1, rel_heights=c(1,6), align="v")
b  <- plot_grid(empty, hist_y, ncol=1, rel_heights=c(1,6))
p_final <- plot_grid(a, b, rel_widths=c(6,1))
p_final # 5X5 55%

################################
###########FigureS3d############
################################
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
genomic=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.txt")
#genomic$patientID=as.character(genomic$patientID)
genomic=left_join(genomic,clin,by="patientID")
genomic%>%group_by(Arm,Response)%>%summarise(mean=median(TMB_uniform))
colnames(genomic)
d=genomic%>%select(c("Arm","Response","CNV_burden","total_LOH_burden","LOH_Del_TSG_burden","TMB_LOH","TMB_CNA_corrected"))
d <- reshape2::melt(d,id.vars=c("Arm","Response"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
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

################################
###########FigureS3e############
################################
library(ggpubr);library(viridis)
genomic=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.txt")
quartiles <- quantile(genomic$CNA_ERBB2, probs=c(.25, .75), na.rm = FALSE)
IQR <- IQR(data$CNA_ERBB2)
Lower <- quartiles[1] - 0.5*IQR
Upper <- quartiles[2] + 0.5*IQR 
data=genomic
data_no_outlier <- subset(data, data$CNA_ERBB2 > Lower & data$CNA_ERBB2 < Upper)
dim(data_no_outlier)
# Add 2d density estimation
figure_font_size=13
ggplot(data_no_outlier, aes(x =CNA_ERBB2, y = LOH_Del_burden)) +
  geom_density_2d_filled()+stat_cor(method = "spearman")+ theme_bw()+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size))
#########################
quartiles <- quantile(genomic$CNA_BRCA2, probs=c(.25, .75), na.rm = FALSE)
IQR <- IQR(data$CNA_BRCA2)
Lower <- quartiles[1] - 0.5*IQR
Upper <- quartiles[2] + 0.5*IQR 
data=genomic
data_no_outlier <- subset(data, data$CNA_BRCA2 > Lower & data$CNA_BRCA2 < Upper)
dim(data_no_outlier)
ggplot(data_no_outlier, aes(x =CNA_BRCA2, y = LOH_Del_TSG_burden)) +
  geom_density_2d_filled()+stat_cor(method = "spearman")+ theme_bw()+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size))

data$
