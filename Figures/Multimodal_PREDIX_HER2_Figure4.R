#####################################
##########Figure4.b #################
#####################################
# density plot https://r-graph-gallery.com/135-stacked-density-graph.html
library(data.table);library(tidyverse)
source("E:/Projects/PREDIX_HER2/Multimodal/Analyses/NMF/src/main.expression_clustering_simple.WR.R")
tpm=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/TMM-normalized-TPM.rds")
gene=readRDS('E:/Projects/PREDIX_HER2/Multimodal/Analyses/NMF/geneset.rds')
purity=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/PREDIX_HER2_ESTIMATE.txt"))
purity=purity[substr(purity$sample,14,14)==0,]
purity$patientID=substr(purity$sample,9,12)%>%as.integer()
purity=purity[order(purity$purity),];purity=purity[!duplicated(purity$patientID),]
#purity=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/WES/PureCN/PREDIX_HER2_PureCN_baseline_purity_curated.csv"))
#purity$patientID=substr(purity$sampleID,9,12)%>%as.integer()
tpm=tpm[gene,]
abundance=as.matrix(t(tpm))
tpm=scale(abundance,center = TRUE, scale = TRUE)
numeric_matrix <- matrix(as.numeric(unlist(tpm)), nrow = nrow(tpm), ncol = ncol(tpm))
# Set the row and column names
rownames(numeric_matrix) <- rownames(tpm)
colnames(numeric_matrix) <- colnames(tpm)
tpm=as.data.frame(tpm)
rna=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.txt")
dna=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
dna$sampleID=NULL
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")

load(file="E:/Projects/PREDIX_HER2/Multimodal/Analyses/NMF/input_data/BayesNMF_output/PREDIX/g.Bayes.RData")
load(file="E:/Projects/PREDIX_HER2/Multimodal/Analyses/NMF/input_data/BayesNMF_output/PREDIX/H.RData")
HER2ADC <- data.frame(sampleID=names(g.Bayes), HER2ADC=paste0("S",g.Bayes))
saveRDS(HER2ADC,file="E:/Projects/PREDIX_HER2/Multimodal/Analyses/NMF/input_data/BayesNMF_output/PREDIX/predix_her2_Her2ADC.rds")
table(HER2ADC$HER2ADC)
color.axis="grey"
sample.ordering <- colnames(H)[order(g.Bayes,decreasing=F)]
res <- get.sample.association.heatmap(H,g.Bayes,1.0)
p2 <- res[[2]] # normalized
p2
# 10X3
### heatmap ###
# totalTMB, CNA burden, LOH Deletion Burden, ER, HER2 IHC, pCR, event, TP53/PIK3CA/ERBB2 mut, PAM50, HER2DX, Taxane, pik3ca score, ABC transporter, exosome

tpm=tpm[sample.ordering,]
library(data.table);library(tidyverse);library(IOBR);library(readxl);library(tidyverse);library(data.table);library(readr);library(readxl)
library(tableone);library(survival);library(tidyverse);library(survminer);library(data.table)
library(ComplexHeatmap);library(circlize);library(RColorBrewer);library(ggstatsplot)
heatmap_gene=colnames(tpm)
df=tpm
df$sampleID=row.names(tpm)
df$patientID=substr(df$sampleID,9,12)%>%as.integer()
df=left_join(df,HER2ADC,by="sampleID")%>%left_join(clin,by="patientID")%>%left_join(rna,by="patientID")%>%
   left_join(dna,by="patientID")%>%left_join(purity,by="patientID")
row.names(df)=df$sampleID
name=c("HER2DX_pCR_likelihood_score","Taxane_response","pik3ca_sig","ABC_transporter","Exosome")
df[,name]=scale(df[,name])
#### HER2DX group
s1 = as.matrix(t(tpm[df$sampleID[df$HER2ADC=="S1"],]))   # subtype1
S1_meta=df[colnames(s1),]
s2 = as.matrix(t(tpm[df$sampleID[df$HER2ADC=="S2"],]))  # subtype1
S2_meta=df[colnames(s2),]
s3 = as.matrix(t(tpm[df$sampleID[df$HER2ADC=="S3"],]))   # subtype1
S3_meta=df[colnames(s3),]

pheno=S1_meta
ha1=HeatmapAnnotation(foo=anno_block(gp= gpar(fill = "#4DE66699"),labels = c("S1 (n = 69)")),
                      TMB=anno_barplot(pheno$totalTMB,ylim =c(0,60),border =F,show_name = FALSE,bar_width = 0.4),
                      KI67=anno_barplot(pheno$MKI67_baseline,ylim =c(0,15),border =F,show_name = FALSE,bar_width = 0.4),
                      Arm=pheno$Arm,
                      PAM50=pheno$sspbc.subtype,
                      ER=pheno$ER,
                      HER2_IHC=pheno$HER2neu1,
                      Response=pheno$Response,
                      Event=pheno$EFS.status,
                      TP53_mut=pheno$coding_mutation_TP53,
                      PIK3CA_mut=pheno$coding_mutation_PIK3CA,
                      ERBB2_mut=pheno$coding_mutation_ERBB2,
                      Taxane=pheno$Taxane_response,
                      HER2DX=pheno$HER2DX_pCR_likelihood_score,
                      PIK3CA=pheno$pik3ca_sig,
                      ABC_TRANSPORTER=pheno$ABC_transporter,
                      Exosome=pheno$Exosome,
                      col=list(Arm=c("DHP"="#8491B4FF","T-DM1"="#91D1C2FF"),
                               PAM50= c("LumA"="#1f78b4", "LumB"="#a6cee3", "Her2"="#fb9a99","Basal"="#e31a1c"),
                               HER2_IHC=c("2+"="#E8E8E8","3+"="black"),
                               ER=c("negative"="#E8E8E8","positive"="black"),
                               Response=c("pCR"="#E8E8E8","RD"="black"),
                               Event=c("0"="#E8E8E8","1"="black"),
                               TP53_mut=c("0"="#E8E8E8","1"="black"),
                               PIK3CA_mut=c("0"="#E8E8E8","1"="black"),
                               ERBB2_mut=c("0"="#E8E8E8","1"="black"),
                               ABC_TRANSPORTER=colorRamp2(c(-2,-1, 0, 1,2), c("#0A2A7C","#6B7BA4","#EAEAEA","#F98181", "#F41D1D")),
                               Taxane=colorRamp2(c(-2,-1, 0, 1,2), c("#0A2A7C","#6B7BA4","#EAEAEA","#F98181", "#F41D1D")),
                               PIK3CA=colorRamp2(c(-2,-1, 0, 1,2), c("#0A2A7C","#6B7BA4","#EAEAEA","#F98181", "#F41D1D")),
                               HER2DX=colorRamp2(c(-2,-1, 0, 1,2), c("#0A2A7C","#6B7BA4","#EAEAEA","#F98181", "#F41D1D")),
                               Exosome=colorRamp2(c(-2,-1, 0, 1,2), c("#0A2A7C","#6B7BA4","#EAEAEA","#F98181", "#F41D1D"))
                               ),
                      na_col="#808080",show_annotation_name = FALSE,
                      border=F)

pheno=S2_meta
ha2=HeatmapAnnotation(foo=anno_block(gp= gpar(fill = "#F39B7FFF"),labels = c("S2 (n = 74)")),
                      TMB=anno_barplot(pheno$totalTMB,ylim =c(0,60),border =F,show_name = FALSE,bar_width = 0.4),
                      KI67=anno_barplot(pheno$MKI67_baseline,ylim =c(0,15),border =F,show_name = FALSE,bar_width = 0.4),
                      Arm=pheno$Arm,
                      PAM50=pheno$sspbc.subtype,
                      ER=pheno$ER,
                      HER2_IHC=pheno$HER2neu1,
                      Response=pheno$Response,
                      Event=pheno$EFS.status,
                      TP53_mut=pheno$coding_mutation_TP53,
                      PIK3CA_mut=pheno$coding_mutation_PIK3CA,
                      ERBB2_mut=pheno$coding_mutation_ERBB2,
                      Taxane=pheno$Taxane_response,
                      HER2DX=pheno$HER2DX_pCR_likelihood_score,
                      PIK3CA=pheno$pik3ca_sig,
                      ABC_TRANSPORTER=pheno$ABC_transporter,
                      Exosome=pheno$Exosome,
                      col=list(Arm=c("DHP"="#8491B4FF","T-DM1"="#91D1C2FF"),
                               PAM50= c("LumA"="#1f78b4", "LumB"="#a6cee3", "Her2"="#fb9a99","Basal"="#e31a1c"),
                               HER2_IHC=c("2+"="#E8E8E8","3+"="black"),
                               ER=c("negative"="#E8E8E8","positive"="black"),
                               Response=c("pCR"="#E8E8E8","RD"="black"),
                               Event=c("0"="#E8E8E8","1"="black"),
                               TP53_mut=c("0"="#E8E8E8","1"="black"),
                               PIK3CA_mut=c("0"="#E8E8E8","1"="black"),
                               ERBB2_mut=c("0"="#E8E8E8","1"="black"),
                               ABC_TRANSPORTER=colorRamp2(c(-2,-1, 0, 1,2), c("#0A2A7C","#6B7BA4","#EAEAEA","#F98181", "#F41D1D")),
                               Taxane=colorRamp2(c(-2,-1, 0, 1,2), c("#0A2A7C","#6B7BA4","#EAEAEA","#F98181", "#F41D1D")),
                               PIK3CA=colorRamp2(c(-2,-1, 0, 1,2), c("#0A2A7C","#6B7BA4","#EAEAEA","#F98181", "#F41D1D")),
                               HER2DX=colorRamp2(c(-2,-1, 0, 1,2), c("#0A2A7C","#6B7BA4","#EAEAEA","#F98181", "#F41D1D")),
                               Exosome=colorRamp2(c(-2,-1, 0, 1,2), c("#0A2A7C","#6B7BA4","#EAEAEA","#F98181", "#F41D1D"))
                      ),
                      na_col="#808080",show_annotation_name = FALSE,
                      border=F)

pheno=S3_meta
ha3=HeatmapAnnotation(foo=anno_block(gp= gpar(fill = "#4D1A6699"),labels = c("S3 (n = 42)")),
                      TMB=anno_barplot(pheno$totalTMB,ylim =c(0,60),border =F,show_name = FALSE,bar_width = 0.4),
                      KI67=anno_barplot(pheno$MKI67_baseline,ylim =c(0,15),border =F,show_name = FALSE,bar_width = 0.4),
                      Arm=pheno$Arm,
                      PAM50=pheno$sspbc.subtype,
                      ER=pheno$ER,
                      HER2_IHC=pheno$HER2neu1,
                      Response=pheno$Response,
                      Event=pheno$EFS.status,
                      TP53_mut=pheno$coding_mutation_TP53,
                      PIK3CA_mut=pheno$coding_mutation_PIK3CA,
                      ERBB2_mut=pheno$coding_mutation_ERBB2,
                      Taxane=pheno$Taxane_response,
                      HER2DX=pheno$HER2DX_pCR_likelihood_score,
                      PIK3CA=pheno$pik3ca_sig,
                      ABC_TRANSPORTER=pheno$ABC_transporter,
                      Exosome=pheno$Exosome,
                      col=list(Arm=c("DHP"="#8491B4FF","T-DM1"="#91D1C2FF"),
                               PAM50= c("LumA"="#1f78b4", "LumB"="#a6cee3", "Her2"="#fb9a99","Basal"="#e31a1c"),
                               HER2_IHC=c("2+"="#E8E8E8","3+"="black"),
                               ER=c("negative"="#E8E8E8","positive"="black"),
                               Response=c("pCR"="#E8E8E8","RD"="black"),
                               Event=c("0"="#E8E8E8","1"="black"),
                               TP53_mut=c("0"="#E8E8E8","1"="black"),
                               PIK3CA_mut=c("0"="#E8E8E8","1"="black"),
                               ERBB2_mut=c("0"="#E8E8E8","1"="black"),
                               ABC_TRANSPORTER=colorRamp2(c(-2,-1, 0, 1,2), c("#0A2A7C","#6B7BA4","#EAEAEA","#F98181", "#F41D1D")),
                               Taxane=colorRamp2(c(-2,-1, 0, 1,2), c("#0A2A7C","#6B7BA4","#EAEAEA","#F98181", "#F41D1D")),
                               PIK3CA=colorRamp2(c(-2,-1, 0, 1,2), c("#0A2A7C","#6B7BA4","#EAEAEA","#F98181", "#F41D1D")),
                               HER2DX=colorRamp2(c(-2,-1, 0, 1,2), c("#0A2A7C","#6B7BA4","#EAEAEA","#F98181", "#F41D1D")),
                               Exosome=colorRamp2(c(-2,-1, 0, 1,2), c("#0A2A7C","#6B7BA4","#EAEAEA","#F98181", "#F41D1D"))
                      ),
                      na_col="#808080",
                      border=F)
#### Plot heatmap
col_fun = colorRamp2(c(-2.5,-1,0,1,2.5), c("#00FF00", "#008000", "#000000", "#800000", "#FF0000"))
ht_list = Heatmap(as.matrix(s1), col = col_fun, name = "Scaled log2(TPM+1)",
                  cluster_rows =T,show_row_names = FALSE,cluster_columns = F,
                  show_row_dend = T, show_column_dend = F,
                  show_column_names = F,
                  bottom_annotation = ha1,
                  row_title_gp = gpar(col = "#FFFFFF00"), width = unit(6.4, "cm"),height= unit(10, "cm"))+
  Heatmap(as.matrix(s2),col=col_fun, show_column_names = F, cluster_columns = F,
          show_column_dend = F,bottom_annotation=ha2,
          show_heatmap_legend = FALSE, width = unit(6.8, "cm"),height= unit(10, "cm"))+
  Heatmap(as.matrix(s3),col=col_fun, show_column_names = F, cluster_columns = F,
          show_column_dend = F,bottom_annotation=ha3,row_names_gp = gpar(fontsize = 2),
          show_heatmap_legend = FALSE, width = unit(3.9, "cm"),height= unit(10, "cm"))+
  rowAnnotation(link = anno_mark(at = c(1,11,21,26,30,36,46,58,61,79,81,83), 
                                 labels =c("SIX2","S100A9","ABCC3","PIK3C2G",
  "KRT19","CCND1","COX6C","COL17A1","ZNF516","ESR1","IL6ST","RARA")))

p0=draw(ht_list, row_title = paste0("DEGs used for BeyasianNMF"),
        annotation_legend_side = "bottom", heatmap_legend_side = "left") 


# 15 X 10 P

#######Description#######
vars=c("Arm","Response","totalTMB","CNV_burden","LOH_Del_burden","sspbc.subtype","ER","HER2neu1","MKI67_baseline","EFS.status",
       "coding_mutation_TP53","coding_mutation_PIK3CA","coding_mutation_ERBB2","purity",
       "HER2DX_pCR_likelihood_score","Taxane_response","pik3ca_sig","ABC_transporter","Exosome","manuTILS")
catvars=c("Arm","Response","sspbc.subtype","ER","HER2neu1","EFS.status","coding_mutation_TP53","coding_mutation_PIK3CA","coding_mutation_ERBB2")
tableOne <- CreateTableOne(vars=vars,strata = c("HER2ADC"),factorVars = catvars,includeNA=T,data = df)
tab_out<-print(tableOne)
#write.csv(tab_out,file='table_to_figure2.csv')

my_comparisons <- list(c("S1", "S2"), c("S1", "S3"), c("S2", "S3"))
df$manuTILS
# 绘图
p <- ggplot(df, aes(x = HER2ADC, y = purity, fill = HER2ADC)) +
  stat_summary(fun = mean, geom = "bar", color = "black", width = 0.6) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  geom_jitter(width = 0.1, alpha = 0.5) +
  labs(x = "Group", y = "Purity") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "none") +
  scale_fill_manual(values = c("steelblue", "tomato", "darkgreen")) +
  # 显示具体 p 值，改用 Wilcoxon
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.format") +
  stat_compare_means(method = "wilcox.test", label.y = 0.9)  # 总体显著性
print(p)
######## Figure.4c ##############
# compare ESR1 ERBB2 MKI67 TOP2A
#################################
library(tidyverse);library(data.table);library(ggpubr);library(ggplot2);library(ggridges);library(coin)
abundance=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/TMM-normalized-TPM.rds")
OUTPUT_ttype <- "E:/Projects/PREDIX_HER2/Multimodal/Analyses/NMF/input_data/BayesNMF_output/PREDIX/"
load(file=paste(OUTPUT_ttype,"g.Bayes.RData",sep=""))
df2 <- data.frame(sample=names(g.Bayes), g.Bayes=g.Bayes)
res=df2
colnames(res)[1]="sampleID"
res$patientID=substr(res$sampleID,9,12)%>%as.integer()
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
res=left_join(res,clin,by="patientID")%>%cbind(t(abundance[c("ERBB2","ESR1","PGR","MKI67","TOP2A","NR2F1","NFE2L2","AXL","GAS6","WNT5A",
                                                             "CDKN1B", "CDKN1A", "SOX9", "ID1", "ID3", "MAPK14", "TGFBR2", "SMAD3"),]))
#"NFE2L2","NR2F1","CD8A"

res$HER2ADC=paste0("S",res$g.Bayes)
#res=res[,c("sampleID","patientID","HER2ADC")]
#write.table(res,file="E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/HER2ADC_PREDIX.txt",quote = F,row.names =F,sep="\t")

kruskal.test(res$ERBB2,res$HER2ADC,distribution = "exact")
kruskal.test(res$ESR1,res$HER2ADC,distribution = "exact")
kruskal.test(res$PGR,res$HER2ADC,distribution = "exact")
kruskal.test(res$MKI67,res$HER2ADC)
kruskal.test(res$MKI67,res$HER2ADC)

d=res[,c("HER2ADC","ERBB2","ESR1","PGR","MKI67","TOP2A","NR2F1","NFE2L2","AXL","GAS6",
         "CDKN1B", "CDKN1A", "SOX9", "ID1", "ID3", "MAPK14", "TGFBR2", "SMAD3")]
d <- reshape2::melt(d,id.vars=c("HER2ADC"))
colnames(d)=c("Her2ADC","mRNA","log2(TPM+1)")
library(ggplot2)
library(cowplot)
# Identity on x-axis https://github.com/ycl6/StackedVlnPlot
d$Her2ADC=factor(d$Her2ADC,levels = c("S3","S2","S1"))
g <- ggplot(d, aes(`log2(TPM+1)`, factor(Her2ADC), fill = Her2ADC)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_fill_manual(breaks = c("S1", "S2", "S3"), 
                    values=c("#4DE66699", "#F39B7FFF", "#4D1A6699")) +
  facet_grid(cols = vars(mRNA), scales = "free_x") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank()) +
  xlab("log2(TPM+1)") + ylab("")
g+geom_boxplot(width=0.1) # median and quartile

# 8X4 40%

##############################
########### Fig2d  ###########
#########HER2ADC GSEA#########
##############################
library (data.table)
library (edgeR)
library (EnsDb.Hsapiens.v86)
library (fgsea)
library (ggbiplot)
library (ggplot2)
library (ggpmisc)
library (ggridges)
library (MASS)
library (org.Hs.eg.db)
library (Pigengene)
library (ReactomePA)
library (readxl)
library (sm)
library (stringr)
library (viridis)
library(DESeq2)
library(edgeR)
library(tidyverse)
library(dplyr)
figure_font_size=12
ensemblToHugo <- readRDS("E:/Projects/PREDIX_HER2/RNAseq/data/readSalmon/tx2gnee2ENTREZID.rds")
ensemblToHugo=ensemblToHugo[,c("ENTREZID","ENSEMBL","GeneName")]
ensemblToHugo=ensemblToHugo[!duplicated(ensemblToHugo$ENTREZID),]
protein=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Resource/gene_with_protein_product.txt"))
protein=protein$symbol
# read salmon count (with batch correction)
count=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/Salmon_count_withbatchcorrection.rds")
txi=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/PREDIX_HER2_baseline_curated_txi.rds")
gene=intersect(protein,row.names(count))
abundance=txi$abundance[gene,colnames(count)]
count=txi$counts[gene,colnames(count)]
length=txi$length[gene,colnames(count)]
txi$abundance=txi$abundance[gene,colnames(count)]
txi$counts=txi$counts[gene,colnames(count)]
txi$length=txi$length[gene,colnames(count)]
# read meta data
batch=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/PREDIX_HER2_baseline_RNA-seq_curated_meta.rds")
batch$batch=paste0("Run",batch$batch,"_",batch$lane)
# read her2adc subtype
load(file="E:/Projects/PREDIX_HER2/Multimodal/Analyses/NMF/input_data/BayesNMF_output/PREDIX/g.Bayes.RData")
HER2ADC <- data.frame(sampleID=names(g.Bayes), HER2ADC=paste0("S",g.Bayes))
table(HER2ADC$HER2ADC)
HER2ADC$S1="Other"
HER2ADC$S1[HER2ADC$HER2ADC=="S1"]="S1"
HER2ADC$S2="Other"
HER2ADC$S2[HER2ADC$HER2ADC=="S2"]="S2"
HER2ADC$S3="Other"
HER2ADC$S3[HER2ADC$HER2ADC=="S3"]="S3"
#batch$batch=paste0("Run",batch$batch)
batch$patientID=as.integer(batch$patientID)
batch=batch[,c("patientID","batch")]
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
meta=data.frame(sampleID=colnames(count));meta$patientID=substr(meta$sampleID,9,12)%>%as.integer()
meta=left_join(meta,clin,by="patientID")%>%
  left_join(HER2ADC,by="sampleID")%>%
  left_join(batch,by="patientID")
meta$Response=factor(meta$Response,levels = c("RD","pCR"))
row.names(meta)=meta$sampleID
all.equal(meta$sampleID,colnames(count))
all.equal(meta$sampleID,colnames(length))
#######################
## DESeq2
library("DESeq2")
dds <- DESeqDataSetFromTximport(txi,
                                colData = meta,
                                design = ~batch+S1) #batch + Response
dim(dds)

smallestGroupSize <- 10
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
table(keep)
dds <- dds[keep,]
#specifying the reference level
library("BiocParallel")
register(MulticoreParam(4))
dds$S1 <- relevel(dds$S1, ref = "Other")
dds <- DESeq(dds)
resultsNames(dds)
resLFC <- lfcShrink(dds,coef=c("S1_S1_vs_Other"), type="apeglm")
resLFC
summary(resLFC)
# filter DGE
res=as.data.frame(resLFC)
res$gene=row.names(res)
# load MSigDB Hallmarks Gene set from resources folder and run enrichment
Hallmark<- gmtPathways("E:/Projects/PREDIX_HER2/Multimodal/Resource/h.all.v2023.2.Hs.symbols.gmt")
geneList=Hallmark
# fgsea - get Normalised Enrichment Score
df=res[!is.na(res$pvalue),]
#df$enrichment_score <- (2*as.numeric(df$log2FoldChange > 0) - 1) * -log10(df$padj) # signed p-value
df$enrichment_score <-df$log2FoldChange
df <- df[order(df$enrichment_score, decreasing = TRUE),]
ranks <- df$enrichment_score
names(ranks) <- df$gene

g <- geneList
q <- fgsea(pathways = g, stats = ranks,minSize=15, maxSize=500, nperm=100000)
S1_GSEA=as.data.frame(q)
S1_GSEA$HER2ADC='S1'

## DESeq2
library("DESeq2")
dds <- DESeqDataSetFromTximport(txi,
                                colData = meta,
                                design = ~batch+S2) #batch + Response
dim(dds)

smallestGroupSize <- 10
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
table(keep)
dds <- dds[keep,]
#specifying the reference level
library("BiocParallel")
register(MulticoreParam(4))
dds$S2 <- relevel(dds$S2, ref = "Other")
dds <- DESeq(dds)
resultsNames(dds)
resLFC <- lfcShrink(dds,coef=c("S2_S2_vs_Other"), type="apeglm")
resLFC
summary(resLFC)
# filter DGE
res=as.data.frame(resLFC)
res$gene=row.names(res)
# load MSigDB Hallmarks Gene set from resources folder and run enrichment
Hallmark<- gmtPathways("E:/Projects/PREDIX_HER2/Multimodal/Resource/h.all.v2023.2.Hs.symbols.gmt")
geneList=Hallmark
# fgsea - get Normalised Enrichment Score
df=res[!is.na(res$pvalue),]
#df$enrichment_score <- (2*as.numeric(df$log2FoldChange > 0) - 1) * -log10(df$pvalue) # signed p-value
df$enrichment_score <-df$log2FoldChange
df <- df[order(df$enrichment_score, decreasing = TRUE),]
ranks <- df$enrichment_score
names(ranks) <- df$gene

g <- geneList
q <- fgsea(pathways = g, stats = ranks,minSize=15, maxSize=500, nperm=100000)
S2_GSEA=as.data.frame(q)
S2_GSEA$HER2ADC='S2'

## DESeq2
library("DESeq2")
dds <- DESeqDataSetFromTximport(txi,
                                colData = meta,
                                design = ~batch+S3) #batch + Response
dim(dds)

smallestGroupSize <- 10
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
table(keep)
dds <- dds[keep,]
#specifying the reference level
library("BiocParallel")
register(MulticoreParam(4))
dds$S3 <- relevel(dds$S3, ref = "Other")
dds <- DESeq(dds)
resultsNames(dds)
resLFC <- lfcShrink(dds,coef=c("S3_S3_vs_Other"), type="apeglm")
resLFC
summary(resLFC)
# filter DGE
res=as.data.frame(resLFC)
res$gene=row.names(res)
# load MSigDB Hallmarks Gene set from resources folder and run enrichment
Hallmark<- gmtPathways("E:/Projects/PREDIX_HER2/Multimodal/Resource/h.all.v2023.2.Hs.symbols.gmt")
geneList=Hallmark
# fgsea - get Normalised Enrichment Score
df=res[!is.na(res$pvalue),]
#df$enrichment_score <- (2*as.numeric(df$log2FoldChange > 0) - 1) * -log10(df$pvalue) # signed p-value
df$enrichment_score <-df$log2FoldChange
df <- df[order(df$enrichment_score, decreasing = TRUE),]
ranks <- df$enrichment_score
names(ranks) <- df$gene

g <- geneList
q <- fgsea(pathways = g, stats = ranks,minSize=15, maxSize=500, nperm=100000)
S3_GSEA=as.data.frame(q)
S3_GSEA$HER2ADC='S3'

#########save results##########
HER2ADC_GSEA=rbind(S1_GSEA,S2_GSEA,S3_GSEA)%>%as.data.frame()
openxlsx::write.xlsx(HER2ADC_GSEA,file="E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/HER2ADC_GSEA.xlsx")

####Figure3h
df=read_excel('E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/HER2ADC_GSEA.xlsx')
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
h <- ggplot(df, aes(x = HER2ADC, y = pathway, size = PVAL, color = NES)) +
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
h
# 8X6
#ggarrange(g,h,nrow = 1,widths=c(2,3))
# 5X6 P
###############################
##########Figure4.e############
###############################
library(readr);library(ggplot2);library(readxl)
df=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure4/pCR_NMF_by_treatment.xlsx")
df=df[df$N>10,]
df$Name=paste0(df$Trial,df$Arm)
df$pCR_rate=round(df$pCR_rate,1)
#df$Name=factor(df$Name,levels = c("NCT02591966AC-TH","TRIOB07TCL","TRIOB07TCH","TRIOB07TCHL",
#              "I-SPY2T-DM1+P","I-SPY2DHP","I-SPY2Paclitaxel + Neratinib","I-SPY2TH-AKTi",
#              "I-SPY2TH","NCT02326974T-DM1+P(site2)","NCT02326974T-DM1+P(site1)","NCT02326974T-DM1+P",
#              "PREDIXT-DM1","PREDIXDHP"))
df$Name=factor(df$Name,levels = c("NCT02326974T-DM1+P(site2)","NCT02326974T-DM1+P(site1)","NCT02326974T-DM1+P",
                                               "PREDIXT-DM1","PREDIXDHP"))
                                  
#pCR
p1=ggplot(df, aes(`NMF`, Name, fill=pCR_rate)) +
  geom_tile(color="grey2", size=0.7) +
  geom_text(aes(label=pCR_rate), color="gray2", size=4) +
  scale_fill_gradientn(
    colours=c('#2D6DB1', 'white', '#DC1623')
  ) +
  labs(x=NULL, y=NULL, fill="") +
  theme_minimal() +
  theme(axis.text=element_text(colour='black', size=9)) +
  theme(axis.text.x=element_text(angle=60, hjust=0.9, vjust=0.9, size=9))
#size
library(readxl);library(ggplot2);library(RColorBrewer);library(ggpubr)
df=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure4/pCR_NMF_by_treatment.xlsx")
df$Name=paste0(df$Trial,df$Arm)
#df=df[df$Name%in%c("NCT02591966AC-TH","TRIOB07TCL","TRIOB07TCH","TRIOB07TCHL",
#                   "I-SPY2T-DM1+P","I-SPY2DHP","I-SPY2Paclitaxel + Neratinib","I-SPY2TH-AKTi",
#                   "I-SPY2TH","NCT02326974T-DM1+P(site2)","NCT02326974T-DM1+P(site1)","NCT02326974T-DM1+P",
#                   "PREDIXT-DM1","PREDIXDHP"),]
#df$Name=factor(df$Name,levels = c("NCT02591966AC-TH","TRIOB07TCL","TRIOB07TCH","TRIOB07TCHL",
#                                  "I-SPY2T-DM1+P","I-SPY2DHP","I-SPY2Paclitaxel + Neratinib","I-SPY2TH-AKTi",
#                                  "I-SPY2TH","NCT02326974T-DM1+P(site2)","NCT02326974T-DM1+P(site1)","NCT02326974T-DM1+P",
#                                  "PREDIXT-DM1","PREDIXDHP"))
df$Name=factor(df$Name,levels = c("NCT02326974T-DM1+P(site2)","NCT02326974T-DM1+P(site1)","NCT02326974T-DM1+P",
                                  "PREDIXT-DM1","PREDIXDHP"))


discrete_colors <- brewer.pal(3, 'Greens')
continuous_palette <- colorRampPalette(discrete_colors)

p2=ggplot(df, aes(`NMF`, Name, fill=N)) +
  geom_tile(color="grey2", size=0.7) +
  geom_text(aes(label=N), color="gray2", size=4) +
  scale_fill_gradientn(
    colours=c("#E5F5E0","#A1D99B","#31A354")
  ) +
  labs(x=NULL, y=NULL, fill="") +
  theme_minimal() +
  theme(axis.text=element_text(colour='black', size=9)) +
  theme(axis.text.x=element_text(angle=60, hjust=0.9, vjust=0.9, size=9))

ggarrange(p2,p1,nrow = 1)
# 8X4 L
# 8X2 L 70%
###############################
##########Figure3.############
###############################
# https://r-graph-gallery.com/37-barplot-with-number-of-observation.html
library(openxlsx);library(ggplot2)
res=openxlsx::read.xlsx("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/pCR_NMF_by_treatment.xlsx")
res$pCR_rate=round(res$pCR_rate)
res=res[res$N>9,]
PREDIX=res[res$Trial=='PREDIX',]
NCT02326974=res[res$Trial=='NCT02326974',]
ISPY2=res[res$Trial=='I-SPY2',]
TRIOB07=res[res$Trial=='TRIOB07',]
NCT02591966=res[res$Trial=='NCT02591966',]
# Increase bottom margin
par(mar=c(6,4,4,4))
# Basic Barplot
my_bar <- barplot(PREDIX$pCR_rate, border=F , names.arg=PREDIX$Arm, 
                  las=2 , 
                  col=PREDIX$color , 
                  ylim=c(0,80) , ylab = "pCR rate (%)",
                  main="PREDIX HER2")
# Add abline
abline(v=c(3.7) , col="grey")

# Add the text 
text(my_bar, PREDIX$pCR_rate+4, paste0(PREDIX$pCR_rate,"%") ,cex=1) 


# Increase bottom margin
par(mar=c(6,4,4,4))
# Basic Barplot
my_bar <- barplot(NCT02326974$pCR_rate, border=F , names.arg=NCT02326974$Arm, 
                  las=2 , 
                  col=NCT02326974$color , 
                  ylim=c(0,80) , ylab = "pCR rate (%)",
                  main="NCT02326974")
# Add the text 
text(my_bar, NCT02326974$pCR_rate+4, paste0(NCT02326974$pCR_rate,"%") ,cex=1) 
abline(v=c(3.7,7.3) , col="grey")
#Legende
legend("topright", legend = c("S1","S2","S3" ) , 
       col = c("#4D1A6699","grey","#4DE66699") , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.8, horiz = FALSE, inset = c(0.05, 0.05))


#########
# ISPY2 #
#########
par(mar=c(6,4,4,4))
# Basic Barplot
my_bar <- barplot(ISPY2$pCR_rate, border=F , names.arg=ISPY2$Arm, 
                  las=2 , 
                  col=ISPY2$color, 
                  ylim=c(0,100) , ylab = "pCR rate (%)",
                  main="I-SPY2")
# Add abline
abline(v=c(1.3,3.7,7.3,9.7) , col="grey")

# Add the text 
text(my_bar, ISPY2$pCR_rate+4, paste0(ISPY2$pCR_rate,"%") ,cex=1) 

#Legende
legend("topright", legend = c("S2","S3" ) , 
       col = c("grey","#4DE66699") , 
       bty = "n", pch=20 , pt.cex = 2, cex = 0.8, horiz = FALSE, inset = c(0.05, 0.05))

#########
# TRIOB07 #
#########
par(mar=c(6,4,4,4))
# Basic Barplot
my_bar <- barplot(TRIOB07$pCR_rate, border=F , names.arg=TRIOB07$Arm, 
                  las=2 , 
                  col=TRIOB07$color, 
                  ylim=c(0,100) , ylab = "pCR rate (%)",
                  main="TRIOB07")
# Add abline
abline(v=c(2.5,4.9) , col="grey")

# Add the text 
text(my_bar, TRIOB07$pCR_rate+4, paste0(TRIOB07$pCR_rate,"%") ,cex=1) 


#########
# NCT02591966 #
#########

par(mar=c(6,4,4,4))
# Basic Barplot
my_bar <- barplot(NCT02591966$pCR_rate, border=F , names.arg=NCT02591966$Arm, 
                  las=2, 
                  col=NCT02591966$color, 
                  ylim=c(0,100) , ylab = "pCR rate (%)",
                  main="NCT02591966")
# Add the text 
text(my_bar, NCT02591966$pCR_rate+4, paste0(NCT02591966$pCR_rate,"%") ,cex=1)

###############################
# Logistic regression P value #
###############################
library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
library(forestplot);library(ggpubr)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
#genomic=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.txt")%>%as.data.frame()
rna=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.txt")%>%as.data.frame()
variable=c("KEGG_ABC_TRANSPORTERS","ERBB2_fusion","chr17q12_fusion","sspbc.subtype","Taxane_response",
           "HER2DX_pCR_likelihood_score","pik3ca_sig","Cytolysis_score")
genomic=left_join(rna,clin,by="patientID")%>%as.data.frame()
str(genomic)
colnames(genomic)
results=Logistic_batch_adjER(genomic,"pCR","Arm",variable,"ER")%>%as.data.frame()

colnames(results)
biomarker=union(results$biomarker[results$TDM1_lr_p<0.105],results$biomarker[results$DHP_lr_p<0.105])
TDM1=results[results$biomarker%in%biomarker,c("biomarker","TDM1_OR","TDM1_lr_p")]
TDM1$group="TDM1"
DHP=results[results$biomarker%in%biomarker,c("biomarker","DHP_OR","DHP_lr_p")]
DHP$group="DHP"
colnames(TDM1)=c("Transcriptomic","OR","Pvalue","group")
colnames(DHP)=c("Transcriptomic","OR","Pvalue","group")
df=rbind(TDM1,DHP)
df$OR=as.numeric(df$OR)
df[1,2]=1/df[1,2]
df$Pvalue=as.numeric(df$Pvalue)
df$log10P=-log10(df$Pvalue)
df$logOR=log(as.numeric(df$OR))
df$group[df$Pvalue>0.105]="NS"
df$group=factor(df$group,levels = c("TDM1","DHP","NS"))
unique(df$Transcriptomic)
df$Transcriptomic=factor(df$Transcriptomic,levels = c("sspbc.subtype","HER2DX_pCR_likelihood_score","ERBB2_fusion","chr17q12_fusion",
                                                      "pik3ca_sig","Taxane_response","KEGG_ABC_TRANSPORTERS"))
df=df[order(df$Transcriptomic),]

baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
ggplot(data=df,aes(Transcriptomic,logOR,fill=group))+
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.7,size=0.25)+
  scale_fill_modelbarplot(name="Treatment Arm")+
  labs(x="Genomic Metrics",y="lnOR")+
  theme_manuscript(base_size = figure_font_size)+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))+
  coord_flip()