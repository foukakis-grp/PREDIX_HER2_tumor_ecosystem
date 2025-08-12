#######################################
##############Total patients###########
#########Supplemental figure5##########
#######################################
rm (list=ls())

#load packages
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
# load list of breast cancer driver genes from resources directory
driverGenes<- scan("E:/Projects/PREDIX_HER2/Multimodal/Resource/breast-cancer-driver-genes.txt", what=character()) #breast cancer
cancergenes<-fread('E:/Projects/PREDIX_HER2/Multimodal/Resource/Census_allMon Nov 28 13_48_06 2022.csv')
cancergenes=cancergenes$`Gene Symbol`
driverGenes=c(driverGenes,cancergenes)
# load Gene Ensembl ID to Hugo ID dictionary
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
#batch$batch=paste0("Run",batch$batch)
batch$patientID=as.integer(batch$patientID)
batch=batch[,c("patientID","batch")]
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
meta=data.frame(sampleID=colnames(count));meta$patientID=substr(meta$sampleID,9,12)%>%as.integer()
meta=left_join(meta,clin,by="patientID")%>%
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
                                design = ~batch+Response) #batch + Response
dim(dds)
#dds=dds[,dds$Arm=="T-DM1"]

smallestGroupSize <- 10
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
table(keep)
dds <- dds[keep,]
#specifying the reference level
library("BiocParallel")
register(MulticoreParam(4))
dds$Response <- relevel(dds$Response , ref = "RD")
dds <- DESeq(dds)
resultsNames(dds)
resLFC <- lfcShrink(dds,coef=c("Response_pCR_vs_RD"), type="apeglm")
resLFC
summary(resLFC)
# filter DGE
res=as.data.frame(resLFC)
res$gene=row.names(res)
#DGE=filter(res,abs(log2FoldChange)>0.5,abs(padj)<0.05)
#DGE[intersect(driverGenes,row.names(DGE)),]

# Annotate full gene list with their Hugo Gene name
annotatedResults <- res
annotatedResults[annotatedResults$padj>=0.05|is.na(annotatedResults$padj),"expression"] <- "notDE"
annotatedResults$expression[is.na(annotatedResults$expression)]="DE"
annotatedResults[annotatedResults$log2FoldChange>0&annotatedResults$expression=="DE","expression"] <- "overexpressed"
annotatedResults[annotatedResults$log2FoldChange<0&annotatedResults$expression=="DE","expression"] <- "underexpressed"
annotatedResults$GeneName=row.names(annotatedResults)
annotatedResults <- left_join(annotatedResults,ensemblToHugo, by="GeneName")
table(is.na(annotatedResults$ENTREZID),annotatedResults$expression)
#=========================================================================
# Differential gene expression and enrichment
# Figure 3a
#=========================================================================
# Characterise expression landscape of driver genes and plot Figure 3a
driverExpression <- annotatedResults[annotatedResults$expression!="notDE",]
driverExpression <- driverExpression[driverExpression$GeneName %in% c(driverGenes),]
driverExpression <- driverExpression[abs(driverExpression$log2FoldChange)>0.5,]
driverExpression[order(driverExpression$log2FoldChange,decreasing = T),]
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
fig3a <- 
  ggplot(driverExpression,aes(x=reorder(GeneName,log2FoldChange),y=(log2FoldChange),color=expression))+
  geom_point(aes(size= -log10(padj)))+
  geom_hline(yintercept=0,color="black", linetype="dotted")+
  labs(x="",y="log FC pCR vs RD")+
  coord_flip()+
  scale_colour_manual(values = c("#FB6542","#375E97"))+
  scale_size_continuous(name=expression(italic(-log[10]~padj)),breaks = c(2:4))+
  guides(color="none")+
  theme_manuscript(base_size = figure_font_size)+
  theme(panel.grid.major.y = element_blank(),
        axis.text.y = element_text(face="italic"),
        plot.margin = unit(c(1.1,1.5,0,0.5), "lines"))
fig3a

#5X4, 55%
#=========================================================================
# Differentially expressed pathways/gene sets in pCR vs RD
# Figure 3b, Extended Figure 5a
#=========================================================================

# a. MSigDB enrichment

# load MSigDB Hallmarks Gene set from resources folder and run enrichment
Hallmark<- gmtPathways("E:/Projects/PREDIX_HER2/Multimodal/Resource/h.all.v2023.2.Hs.symbols.gmt")
geneList=Hallmark
# fgsea - get Normalised Enrichment Score
df=annotatedResults[!is.na(annotatedResults$pvalue),]
#df$enrichment_score <- (2*as.numeric(df$log2FoldChange > 0) - 1) * -log10(df$pvalue) # signed p-value
df$enrichment_score <-df$log2FoldChange
df <- df[order(df$enrichment_score, decreasing = TRUE),]
ranks <- df$enrichment_score
names(ranks) <- df$GeneName

g <- geneList
q <- fgsea(pathways = g, stats = ranks,minSize=15, maxSize=500, nperm=100000)
#topPathwaysUp <- q[ES > 0][head(order(pval), n=15), pathway]
#topPathwaysDown <- q[ES < 0][head(order(pval), n=15), pathway]
q=as.data.frame(q)

# create Figure 3b - MSigDB gene set enrichment
q$pathway <- gsub("HALLMARK_"," ",q$pathway)
q$pathway <- gsub("_"," ",q$pathway)
q$pathway <- str_to_sentence(q$pathway)
q$pathway <- gsub("E2f","E2F",q$pathway)
q$pathway <- gsub("Ifn","IFN",q$pathway)
q$pathway <- gsub("G2m","G2M",q$pathway)
q$pathway <- gsub("Il6 jak stat3","IL6 JAK STAT3",q$pathway)
q$pathway <- gsub("Myc","MYC",q$pathway)
q$pathway <- gsub("Mtorc1","MTORC1",q$pathway)
q$pathway <- gsub("Tnfa signaling via nfkb","TNFA signaling via NFKB",q$pathway)
q$pathway <- gsub("Il2 stat","IL2 STAT",q$pathway)
q <- q[order(q$NES),]
q$yax <- ifelse(q$NES > 0, -0.02, 0.02)
q$col <- ifelse(q$NES > 0, "blue","red")

q=q[q$pval<0.2,]
fig3b <- 
  ggplot(q,aes(y=NES,x=reorder(pathway,NES),label =pathway))+
  geom_text(aes(y = yax,hjust = NES > 0),size=(figure_font_size)/(14/5))+ 
  geom_bar(stat="identity",aes(fill=col))+
  geom_hline(yintercept = 0)+
  coord_flip(ylim = c(-2,2.5))+
  labs(y = "Normalised enrichment score", x = "",title="")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="pCR vs RD",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_manuscript(base_size = figure_font_size)+
  theme(panel.grid.major = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(0.1,0,0.05,0), "lines"))
fig3b
#4X5,portrait,60%

# b. Reactome enrichment
res <- annotatedResults[!is.na(annotatedResults$pvalue),]
#res$enrichment_score <- (2*as.numeric(df$log2FoldChange > 0) - 1) * -log10(df$pvalue) # signed p-value
df$enrichment_score <-df$log2FoldChange
res<- res[order(res$enrichment_score, decreasing = TRUE),]
res=res[!is.na(res$ENTREZID),]
res=res[!duplicated(res$ENTREZID),]
ranks=res$enrichment_score
names(ranks)=res$ENTREZID
gse <- gsePathway(ranks , nPerm=1000, minGSSize=120, pvalueCutoff=0.05, pAdjustMethod="BH", verbose=FALSE)
res <- as.data.frame(gse)

# sort reactome enrichment by NES
e     <- res[order(res$NES),]
e$col <- ifelse(e$NES > 0, "blue","red")
e$yax <- ifelse(e$NES > 0, -0.02, 0.02)
e     <- e[abs(e$NES)>1.6,]
e$Description <- str_to_sentence(e$Description)
e$Description <- gsub("Dna","DNA",e$Description)
e$Description <- gsub("Esr","ESR",e$Description)
e$Description <- gsub("hiv","HIV",e$Description)
e$Description <- gsub("g1","G1",e$Description)
e$Description <- gsub("s phase","S phase",e$Description)
e$Description <- gsub("s trans","S trans",e$Description)
e$Description<-gsub("and","&",e$Description)
e$Description<-gsub("Rrna","rRNA",e$Description)
e$Description<-gsub("Neurotransmitter receptors &","",e$Description)

# Plot Extended Figure 5a: Reactome pathways associated with pCR vs RD
eFig5a <- 
  ggplot(e,aes(y=NES,x=reorder(Description,NES),label = Description))+
  geom_text(aes(y = yax,hjust = NES > 0),size=(figure_font_size-1)/(14/5))+ 
  geom_bar(stat="identity",aes(fill=col),width = 0.9)+
  geom_hline(yintercept = 0)+
  coord_flip()+
  labs(y = "Normalised enrichment score", x = "",title="pCR vs residual disease")+
  scale_fill_manual(values = c("#FB6542","#375E97"),name="Expression: pCR",labels=c("Up","Down"))+
  scale_x_discrete(breaks = NULL)+
  theme_manuscript(base_size = figure_font_size)+
  theme(panel.grid.major = element_blank(),
        axis.line.y = element_blank(),
        plot.margin = unit(c(0.1,0,0.2,0), "lines"))
eFig5a

##################################################
####################DEGs##########################
##################################################
#Figure3.A
#collect DEG from two arm, respectively
#load packages
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
# load list of breast cancer driver genes from resources directory
driverGenes<- scan("E:/Projects/PREDIX_HER2/Multimodal/Resource/breast-cancer-driver-genes.txt", what=character()) #breast cancer
cancergenes<-fread('E:/Projects/PREDIX_HER2/Multimodal/Resource/Census_allMon Nov 28 13_48_06 2022.csv')
cancergenes=cancergenes$`Gene Symbol`
driverGenes=c(driverGenes,cancergenes)
# load Gene Ensembl ID to Hugo ID dictionary
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
batch$batch=paste0("Run",batch$batch,"_",batch$lane);
#batch$batch=paste0("Run",batch$batch);
batch$patientID=as.integer(batch$patientID)
batch=batch[,c("patientID","batch")]
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
meta=data.frame(sampleID=colnames(count));meta$patientID=substr(meta$sampleID,9,12)%>%as.integer()
meta=left_join(meta,clin,by="patientID")%>%
  left_join(batch,by="patientID")
meta$Response=factor(meta$Response,levels = c("RD","pCR"))
row.names(meta)=meta$sampleID
all.equal(meta$sampleID,colnames(count))
all.equal(meta$sampleID,colnames(length))
library("DESeq2")
dds <- DESeqDataSetFromTximport(txi,
                                colData = meta,
                                design = ~batch+Response) #batch + Response
dds$Response <- relevel(dds$Response , ref = "RD")
smallestGroupSize <- 10
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
table(keep)
dds <- dds[keep,]
TDM1=dds[,dds$Arm=="T-DM1"]
DHP=dds[,dds$Arm=="DHP"]
#specifying the reference level
library("BiocParallel")
register(MulticoreParam(4))
# for T-DM1
res_TDM1 <- DESeq(TDM1)
res_TDM1 <- lfcShrink(res_TDM1,coef=c("Response_pCR_vs_RD"), type="apeglm")
plotMA(res_TDM1, ylim=c(-2,2))
summary(res_TDM1)
res_TDM1=as.data.frame(res_TDM1)
res_TDM1$gene=row.names(res_TDM1)
# for T-DM1 ER-adjusted
TDM1$ER=as.factor(TDM1$ER)
design(TDM1) <- ~ batch+ER+Response
TDM1_ER=DESeq(TDM1)
TDM1_ER<- lfcShrink(TDM1_ER,coef=c("Response_pCR_vs_RD"), type="apeglm")
summary(TDM1_ER)
res_TDM1_ER=as.data.frame(TDM1_ER)
res_TDM1_ER$gene=row.names(res_TDM1_ER)
# for DHP
res_DHP <- DESeq(DHP)
res_DHP <- lfcShrink(res_DHP ,coef=c("Response_pCR_vs_RD"), type="apeglm")
summary(res_DHP)
res_DHP=as.data.frame(res_DHP)
res_DHP$gene=row.names(res_DHP)
# for DHP ER-adjusted
DHP$ER=as.factor(DHP$ER)
design(DHP) <- ~ batch+ER+Response
DHP_ER=DESeq(DHP)
DHP_ER <- lfcShrink(DHP_ER,coef=c("Response_pCR_vs_RD"), type="apeglm")
summary(DHP_ER)
res_DHP_ER=as.data.frame(DHP_ER)
res_DHP_ER$gene=row.names(res_DHP_ER)
# gene interest 
genes_to_showname =intersect(union(res_TDM1$gene[abs(res_TDM1$log2FoldChange)>0.5&res_TDM1$padj<0.05],
                                   res_DHP$gene[abs(res_DHP$log2FoldChange)>0.5&res_DHP$padj<0.05]),driverGenes) 
genes_to_showname_ER=intersect(union(res_TDM1_ER$gene[abs(res_TDM1_ER$log2FoldChange)>0.5&res_TDM1_ER$padj<0.05],
                                     res_DHP_ER$gene[abs(res_DHP_ER$log2FoldChange)>0.5&res_DHP_ER$padj<0.05])) 

setdiff(union(res_TDM1$gene[abs(res_TDM1$log2FoldChange)>0.5&res_TDM1$padj<0.05],
              res_DHP$gene[abs(res_DHP$log2FoldChange)>0.5&res_DHP$padj<0.05]))


all.equal(res_TDM1$gene,res_TDM1_ER$gene,res_DHP$gene,res_DHP_ER$gene)
require(openxlsx)
list_of_datasets <- list("DEG_total"=res,"DEG_TDM1" = res_TDM1, "DEG_DHP" = res_DHP,
                         "DEG_TDM1_ERadjusted"=res_TDM1_ER,"DEG_DHP_ERadjusted"=res_DHP_ER)
write.xlsx(list_of_datasets, file ="E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEG_All_RNAseq.xlsx")
###########################################
#   prepare for scatter plot Figure 3A    #
###########################################
require(openxlsx);library(ggplot2);library(data.table);library(tidyverse);library(glmmSeq)
driverGenes<- scan("E:/Projects/PREDIX_HER2/Multimodal/Resource/breast-cancer-driver-genes.txt", what=character()) #breast cancer
cancergenes<-fread('E:/Projects/PREDIX_HER2/Multimodal/Resource/Census_allMon Nov 28 13_48_06 2022.csv')
cancergenes=cancergenes$`Gene Symbol`
driverGenes=c(driverGenes,cancergenes)
res_TDM1=openxlsx::read.xlsx("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEG_All_RNAseq.xlsx",sheet=4)
res_TDM1$driver="No"
res_TDM1$driver[res_TDM1$gene%in%driverGenes]="Yes"
res_DHP=openxlsx::read.xlsx("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEG_All_RNAseq.xlsx",sheet=6)
nrow(res_TDM1[!is.na(res_TDM1$padj)&res_TDM1$padj<0.05&res_TDM1$log2FoldChange>0.5,])
nrow(res_TDM1[!is.na(res_TDM1$padj)&res_TDM1$padj<0.05&res_TDM1$log2FoldChange<(-0.5),])
nrow(res_DHP[!is.na(res_DHP$padj)&res_DHP$padj<0.05&res_DHP$log2FoldChange>0.5,])
nrow(res_DHP[!is.na(res_DHP$padj)&res_DHP$padj<0.05&res_DHP$log2FoldChange<(-0.5),])

genes_to_showname =intersect(union(res_TDM1$gene[abs(res_TDM1$log2FoldChange)>0.5&res_TDM1$padj<0.05],
                                   res_DHP$gene[abs(res_DHP$log2FoldChange)>0.5&res_DHP$padj<0.05]),driverGenes) 
genes_to_showname=union(genes_to_showname,union(res_TDM1$gene[order(res_TDM1$log2FoldChange)][1:5],
                        res_TDM1$gene[order(-res_TDM1$log2FoldChange)][1:5]))
genes_to_showname=union(genes_to_showname,union(res_DHP$gene[order(res_DHP$log2FoldChange)][1:5],
                                                res_DHP$gene[order(-res_DHP$log2FoldChange)][1:5]))
genes_to_showname<- intersect(union(res_TDM1$gene[abs(res_TDM1$log2FoldChange)>0.5&res_TDM1$padj<0.05],
                                                     res_DHP$gene[abs(res_DHP$log2FoldChange)>0.5&res_DHP$padj<0.05]),genes_to_showname) 
genes_to_showname=union(genes_to_showname,c("ABCC12","ABCC3"))
all.equal(res_TDM1$gene,res_DHP$gene)
data=data.frame(ID=res_TDM1$gene,x=res_DHP$log2FoldChange,x_q=res_DHP$padj,y=res_TDM1$log2FoldChange,y_q=res_TDM1$padj)
data[is.na(data)] <- 0.9
genes_to_showname
intersect(genes_to_showname,data$ID)

genes_to_showname=intersect(genes_to_showname,data$ID)
data$combined_sig=0
data$combined_sig[data$x_q<0.05&data$y_q>0.05&abs(data$x)>0.5]=1 # sig only in one arm
data$combined_sig[data$x_q>0.05&data$y_q<0.05&abs(data$y)>0.5]=1 # sig only in one arm
data$combined_sig[data$x_q<0.05&data$y_q<0.05&abs(data$x)>0.5&abs(data$y)>0.5]=2 # sig for both arm
table(data$combined_sig)
data$group[data$combined_sig==0]="Notsig"
data$group[data$x_q<0.05&data$combined_sig==1]="unique DEG in DHP arm"
data$group[data$y_q<0.05&data$combined_sig==1]="unique DEG in T-DM1 arm"
data$group[data$combined_sig==2&data$x*data$y<0]="opposite DEG"
data$group[data$combined_sig==2&data$x*data$y>0]="shared DEG"
table(data$group)
range(data$x);range(data$y)
data$gene_label <- ""
row.names(data)=data$ID
data[genes_to_showname, ]$gene_label <- genes_to_showname
colours = c('grey', 'green3', 'gold3', 'blue')
fig3a=data

annot <- lapply(genes_to_showname, function(i) {
  row <- data[i, ]
  x <- row$x
  y <- row$y
  z <- sqrt(x^2 + y^2)
  list(x = x, y = y,
       text = i, textangle = 0, ax = x/z*75, ay = -y/z*75,
       font = list(color = "black", size =11),
       arrowcolor = "black", arrowwidth = 0.5, arrowhead = 0, arrowsize = 1.5,
       xanchor = "auto", yanchor = "auto")
})
library(plotly)
p <- plot_ly(data = data, x = ~x, y = ~y, type = 'scatter', 
             mode = 'markers',
             color = ~group, colors = colours,
             marker = list(size = 7, 
                           line = list(width = 0.25, color = 'white')),
             text = data$gene_label, hoverinfo = 'text') %>%
  layout(annotations = annot,
         xaxis = list(title = "DHP Arm logFC(pCR vs RD)",
                      color = 'black'),
         yaxis = list(title = "T-DM1 Arm logFC(pCR vs RD)",
                      color = 'black'),
         font = list(size = 12),
         legend = list(x = 0, y = 1, font = list(color = 'black'))) %>%
  config(edits = list(annotationPosition = FALSE,
                      annotationTail = TRUE,
                      annotationText = TRUE),
         toImageButtonOptions = list(format = "svg"))
p



###########################
#############GSEA##########
#########Figure3.B#########
###########################
library(fgsea)
Hallmark<- gmtPathways("E:/Projects/PREDIX_HER2/Multimodal/Resource/h.all.v2023.2.Hs.symbols.gmt")
geneList=Hallmark
# fgsea - get Normalised Enrichment Score
# For DHP
df=res_DHP[!is.na(res_DHP$padj),]
#df$enrichment_score <- (2*as.numeric(df$log2FoldChange > 0) - 1) * -log10(df$padj) # signed p-value
df$enrichment_score <-df$log2FoldChange
df <- df[order(df$enrichment_score, decreasing = TRUE),]
ranks <- df$enrichment_score
names(ranks) <- df$gene
q <- fgsea(pathways = geneList, stats = ranks,minSize=15, maxSize=500, nperm=100000)
gse_dhp=as.data.frame(q)
gse_dhp$pathway <- gsub("HALLMARK_","",gse_dhp$pathway)
row.names(gse_dhp)=gse_dhp$pathway
# For T-DM1
df=res_TDM1[!is.na(res_TDM1$padj),]
#df$enrichment_score <- (2*as.numeric(df$log2FoldChange > 0) - 1) * -log10(df$padj) # signed p-value
df$enrichment_score <-df$log2FoldChange
df <- df[order(df$enrichment_score, decreasing = TRUE),]
ranks <- df$enrichment_score
names(ranks) <- df$gene
q <- fgsea(pathways = geneList, stats = ranks,minSize=15, maxSize=500, nperm=100000)
gse_tdm1=as.data.frame(q)
gse_tdm1$pathway <- gsub("HALLMARK_","",gse_tdm1$pathway)
row.names(gse_tdm1)=gse_tdm1$pathway
# integrate
term=union(gse_dhp$pathway[gse_dhp$padj<0.05&abs(gse_dhp$NES)>1],gse_tdm1$pathway[gse_tdm1$padj<0.05&abs(gse_tdm1$NES)>1])
gse_dhp=gse_dhp[term,]
gse_tdm1=gse_tdm1[term,]
fig3b=cbind(gse_dhp,gse_tdm1)

# Define the desired order for 'pathway'
gse_dhp=gse_dhp[order(gse_dhp$NES),]
desired_order <-gse_dhp$pathway  # Define your pathway names in the desired order
# Reorder 'pathway' based on the desired order
gse_dhp$pathway <- factor(gse_dhp$pathway, levels = desired_order)
gse_tdm1$pathway <- factor(gse_tdm1$pathway, levels = desired_order)
# ploting
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
dhp <- 
  ggplot(gse_dhp,aes(x=pathway,y=(NES),color=NES))+
  geom_point(aes(size= -log10(padj)))+
  geom_hline(yintercept=0,color="black", linetype="dotted")+
  labs(x="",y="NES")+
  coord_flip()+
  scale_colour_gradient2(low = "#375E97",mid = "white", high = "#FB6542")+
  scale_size_continuous(name=expression(italic(-log[10]~padj)),breaks = c(1:3))+
  theme_manuscript(base_size = figure_font_size)+
  theme(axis.text.y = element_text(face="italic"),
        plot.margin = unit(c(1.1,1.5,0,0.5), "lines"))
tdm1 <- 
  ggplot(gse_tdm1,aes(x=pathway,y=(NES),color=NES))+
  geom_point(aes(size= -log10(padj)))+
  geom_hline(yintercept=0,color="black", linetype="dotted")+
  labs(x="",y="NES")+
  coord_flip()+
  scale_colour_gradient2(low = "#375E97",mid = "white", high = "#FB6542")+
  scale_size_continuous(name=expression(italic(-log[10]~padj)),breaks = c(1:3))+
  theme_manuscript(base_size = figure_font_size)+
  theme(axis.text.y = element_text(face="italic"),
        plot.margin = unit(c(1.1,1.5,0,0.5), "lines"))

ggarrange(dhp,tdm1,nrow = 1,
          font.label = list(size = figure_font_size, family="Helvetica"),
          common.legend =T)

# landscape 11X6
###############################
###########KEGG pathway########
###############################
library(IOBR);library(data.table)
tpm=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/TMM-normalized-TPM.rds")
txi=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/PREDIX_HER2_baseline_curated_txi.rds")
tpm=txi$abundance
res=calculate_sig_score(eset = tpm, signature = kegg, method = "ssgsea",parallel.size=16)
res[,2:ncol(res)]=scale(res[,2:ncol(res)])
kegg=colnames(res)[2:ncol(res)]
res$patientID=substr(res$ID,9,12)%>%as.integer()
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
res=left_join(res,clin,by="patientID")
res_TDM1=batch_wilcoxon(
  res[res$Arm=="T-DM1",],
  target = "Response",
  feature =kegg,
  feature_manipulation = FALSE
)

res_DHP=batch_wilcoxon(
  res[res$Arm=="DHP",],
  target = "Response",
  feature =kegg,
  feature_manipulation = FALSE
)
require(openxlsx)
list_of_datasets <- list("KEGG_TDM1"=res_TDM1,
                         "KEGG_DHP"=res_DHP)
write.xlsx(list_of_datasets, file ="E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/KEGG_ssGSEA.xlsx")


#install KEGG_db https://github.com/YuLab-SMU/clusterProfiler/issues/561
library(clusterProfiler);library(KEGG.db);library(ReactomePA)
#KEGG pathway over-representation analysis
protein=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Resource/gene_with_protein_product.txt"))
protein$gene=protein$symbol
res_DHP=left_join(res_DHP,protein[,c("gene",'entrez_id')],by="gene")
res_TDM1=left_join(res_TDM1,protein[,c("gene",'entrez_id')],by="gene")
genelist=res_DHP$entrez_id[res_DHP$log2FoldChange>(0.5)&res_DHP$padj<0.05]%>%as.character()
genelist=res_TDM1$entrez_id[res_TDM1$log2FoldChange>0&res_TDM1$padj<0.05]%>%na.omit()%>%as.character()
kegg <- enrichKEGG(gene = genelist,use_internal_data=T,pvalueCutoff = 0.05,pAdjustMethod="BH")
kegg=kegg@result%>%as.data.frame()%>%dplyr::filter(p.adjust<0.1)

x <- enrichPathway(gene=genelist, pvalueCutoff = 0.05, readable=TRUE)
x=x@result
##########Figure3.c##########
library(networkD3);library(tidyr);library(tibble);library(dplyr)    # data manipulation
library(UpSetR);library(ggplot2);library(plyr);library(gridExtra);library(ggpubr);library(data.table)
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.txt")
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
df=left_join(df,clin,by="patientID")
table(df$sspbc.subtype,df$Response)
fig3c=df[,c("patientID","Arm","sspbc.subtype","Response")]

df$sspbc.bin="Non-HER2"
df$sspbc.bin[df$sspbc.subtype=="Her2"]="HER2-enriched"
df$sspbc.bin=factor(df$sspbc.bin,levels = c("Non-HER2","HER2-enriched"))
res<- glm(as.numeric(pCR) ~ sspbc.bin+ER, family = "binomial", data = df[df$Arm=="T-DM1",])
ShowRegTable(res)

table(df$sspbc.subtype,df$Response)
df1=df[,c("Arm","sspbc.subtype")];colnames(df1)=c('source','target')
df2=df[,c("sspbc.subtype","Response")];colnames(df2)=c('source','target')
links=rbind(df1,df2)
links$color=links$source

my_color <- 'd3.scaleOrdinal() .domain(["DHP","T-DM1","LumA","LumB","Her2","Basal","RD","pCR"]) .range(["#8491B4FF","#91D1C2FF","#1f78b4","#a6cee3","#fb9a99","#e31a1c","#fdb462","#00A087FF"])'

links$source <- as.character(links$source)
links$target<- as.character(links$target)
table(links$target)
nodes <- data.frame(name = unique(c(links$source,links$target)))
nodes$color=nodes$name

links$source <- match(links$source, nodes$name) - 1
links$target <- match(links$target, nodes$name) - 1
links$value <- 1 # add also a value

sankeyNetwork(Links = links, Nodes = nodes, Source = 'source',colourScale=my_color,LinkGroup = NULL,
              Target = 'target', Value = 'value', NodeID = 'name',NodeGroup="color",
              height ="600",width ="800",nodeWidth ="50")


###############################
##########Figure3.d.e##########
###############################
library(data.table);library(tidyverse);library(ggpubr)
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.txt")
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
df=left_join(df,clin,by="patientID")
#d=df%>%select(c("Arm","Response","HER2DX_IGG","HER2DX_prolif","HER2DX_luminal","HER2DX_HER2_amplicon"))
d=df[,c("Arm","Response","HER2DX_IGG","HER2DX_prolif","HER2DX_luminal","HER2DX_HER2_amplicon")]
colnames(d)=c("Arm","Response","IGG","Proliferation","Luminal","HER2_amplicon")
d <- reshape2::melt(d,id.vars=c("Arm","Response"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
Fig3c <-
  ggplot(d,aes(x=Arm,y=value,fill=Response))+
  geom_boxplot(outlier.size = 0.2, width=0.8)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Treatment Response")+
  stat_compare_means(aes(group=Response),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.98)+
  labs(y="HER2DX signature score",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.2,0.2,0.5), "lines"))
Fig3c
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.txt")
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
df=left_join(df,clin,by="patientID")
Fig3d=ggdensity(df[df$Arm=="DHP",], title="DHP Arm",x=c("HER2DX_pCR_likelihood_score"), 
                add = "mean", rug = TRUE,
                combine = T,
                color = "Response", fill = "Response",
                palette = c("#fdb462","#00A087FF"))
wilcox.test(HER2DX_pCR_likelihood_score~Response,df[df$Arm=="DHP",])

Fig3e=ggdensity(df[df$Arm=="T-DM1",], title="T-DM1 Arm",x=c("HER2DX_pCR_likelihood_score"), 
                add = "mean", rug = TRUE,
                combine = T,
                color = "Response", fill = "Response",
                palette = c("#fdb462","#00A087FF"))
wilcox.test(HER2DX_pCR_likelihood_score~Response,df[df$Arm=="T-DM1",])

ggarrange(
  Fig3c,                # First row with line plot
  # Second row with box and dot plots
  ggarrange(Fig3d, Fig3e, ncol = 2, labels = c("D", "E")), 
  nrow = 2, 
  labels = "C"       # Label of the line plot
) 

res<- glm(as.numeric(pCR) ~ HER2DX_pCR_likelihood_score, family = "binomial", data = df[df$Arm=='DHP',])
ShowRegTable(res)

res<- glm(as.numeric(pCR) ~ HER2DX_pCR_likelihood_score, family = "binomial", data = df[df$Arm=='T-DM1',])
ShowRegTable(res)
#8X8 30%
#####################################
##########Figure3.f #################
#####################################
library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
library(forestplot);library(ggpubr)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
rna=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.txt')
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
data=left_join(rna,clin,by="patientID")%>%as.data.frame()
data$ERBB2_fusion=as.factor(data$ERBB2_fusion)
data$chr17q12_fusion=as.factor(data$chr17q12_fusion)
variable=c("Taxane_response","HER2DX_pCR_likelihood_score","pik3ca_sig","ABC_transporter","Apoptosis","Lysosome","Endocytosis","EMT","Exosome",
           "Oxidative_phosphorylation","Purine_metabolism","Citrate_cycles","Glutathione_metabolism","Fatty_acid_metabolism",
           "Glycolysis","Hypoxia")
norm_variable=c("pik3ca_sig","ABC_transporter","Apoptosis","Lysosome","Endocytosis","EMT","Exosome",
                "Oxidative_phosphorylation","Purine_metabolism","Citrate_cycles","Glutathione_metabolism","Fatty_acid_metabolism",
                "Glycolysis","Hypoxia")
data[,norm_variable]=scale(data[,norm_variable]) 


results=Logistic_batch_adjER(data,"pCR","Arm",variable,"ER")%>%as.data.frame()
fig3f=results
TDM1=results[,c("biomarker","TDM1_OR","TDM1_lr_p")]
TDM1$group="TDM1"
DHP=results[,c("biomarker","DHP_OR","DHP_lr_p")]
DHP$group="DHP"
colnames(TDM1)=c("Signature","OR","Pvalue","group")
colnames(DHP)=c("Signature","OR","Pvalue","group")
df=rbind(TDM1,DHP)
df$Pvalue=as.numeric(df$Pvalue)
df$log10P=-log10(df$Pvalue)
df$logOR=log(as.numeric(df$OR))
df$group=factor(df$group,levels = c("TDM1","DHP"))
unique(df$Signature)
df$Signature=factor(df$Signature,levels =c("Hypoxia","Glycolysis","Fatty_acid_metabolism","Glutathione_metabolism",
                                           "Citrate_cycles","Purine_metabolism","Oxidative_phosphorylation","Apoptosis","EMT","Exosome",
                                           "Endocytosis","Lysosome","ABC_transporter","pik3ca_sig",
                                           "HER2DX_pCR_likelihood_score","Taxane_response"))
df=df[order(df$Signature),]
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=13
ggplot(data=df,aes(Signature,logOR,fill=group))+
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.7,size=0.25)+
  scale_fill_modelbarplot(name="Treatment Arm")+
  labs(x="Gene Signature",y="lnOR")+
  theme_manuscript(base_size = figure_font_size)+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.5,0.5,1), "lines"))+
  coord_flip()+scale_y_continuous(breaks = c(-0.5, 0, 0.5))


#####################################
##########Figure3.g #################
#####################################
library(data.table)
library(tidyverse)
library(tableone)
library(forestploter)
require(openxlsx)
## continuous variable ##
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
rna=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.txt')
genomic=inner_join(rna,clin,by="patientID")%>%as.data.frame()
variable=c("Taxane_response","HER2DX_pCR_likelihood_score","pik3ca_sig","ABC_transporter","Exosome","Hypoxia","EMT")
genomic=slice_metric(genomic,variable,20,80,5)
## re-run the batch logistic
selected_columns <- grep("_per_", colnames(genomic), value = TRUE)
Interact_result=Logistic_batch_adjER(genomic,"pCR","Arm",selected_columns ,"ER")%>%as.data.frame()
genomic$Arm=factor(genomic$Arm,levels = c("DHP","T-DM1"))
genomic=as.data.frame(genomic)
res=Logistic_batch_continuous_subgroup(genomic,selected_columns)%>%as.data.frame()

res_continuous=res[res$biomarker%in%c("pik3ca_sig_per_80","ABC_transporter_per_65","Exosome_per_45","Hypoxia_per_80"),]
Interact_continuous=Interact_result[Interact_result$biomarker%in%c("pik3ca_sig_per_80","ABC_transporter_per_65","Exosome_per_45","Hypoxia_per_80"),]
continuous=cbind(res_continuous,Interact_continuous)
# merge
library(openxlsx)
list_of_datasets <- list("rna" = continuous)
write.xlsx(list_of_datasets, file = "E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure3/Predictive_biomarker.xlsx")
table(genomic$Hypoxia_per_80,genomic$Arm)
#whole <- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = genomic[genomic$pik3ca_sig_per_80=="Low",])
#ShowRegTable(whole)
#whole <- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = genomic[genomic$pik3ca_sig_per_80=="High",])
#ShowRegTable(whole)
# forest plot #
library(data.table)
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure3/Subgroup_Forestplot.csv")
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



library(openxlsx)
data_list=list("fig3a"=fig3a%>%as.data.frame(),
               "fig3b"=fig3b%>%as.data.frame(),
               "fig3c"=fig3c%>%as.data.frame(),
               "fig3f"=fig3f%>%as.data.frame())
openxlsx::write.xlsx(data_list,file='E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/Figure3.xlsx')








