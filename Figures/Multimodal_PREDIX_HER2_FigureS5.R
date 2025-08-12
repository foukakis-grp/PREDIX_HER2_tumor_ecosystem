########HER2-enriched########
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
df=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.txt")
#batch$batch=paste0("Run",batch$batch)
batch$patientID=as.integer(batch$patientID)
batch=batch[,c("patientID","batch")]
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
meta=data.frame(sampleID=colnames(count));meta$patientID=substr(meta$sampleID,9,12)%>%as.integer()
meta=left_join(meta,clin,by="patientID")%>%
  left_join(batch,by="patientID")%>%
  left_join(df,by="patientID")
meta$Response=factor(meta$Response,levels = c("RD","pCR"))
row.names(meta)=meta$sampleID
all.equal(meta$sampleID,colnames(count))
all.equal(meta$sampleID,colnames(length))
########for DHP HER2-enriched########
## DESeq2
library("DESeq2")
meta$Response=as.factor(meta$Response)
#meta$batch=as.factor(meta$batch)
dds <- DESeqDataSetFromTximport(txi,
                                colData = meta,
                                design = ~batch+Response) #batch + Response
dim(dds)
smallestGroupSize <- 10
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
table(keep)
dds <- dds[keep,]
TDM1_Her2=dds[,dds$Arm=="T-DM1"&dds$sspbc.subtype=="Her2"]
TDM1_Her2$batch=as.factor(TDM1_Her2$batch)
design(TDM1_Her2) <- ~ Response

DHP_Her2=dds[,dds$Arm=="DHP"&dds$sspbc.subtype=="Her2"]
DHP_Her2$batch=as.factor(DHP_Her2$batch)
design(DHP_Her2) <- ~ Response
########for T-DM1 HER2-enriched########
res_TDM1_Her2 <- DESeq(TDM1_Her2)
res_TDM1_Her2 <- lfcShrink(res_TDM1_Her2,coef=c("Response_pCR_vs_RD"), type="apeglm")
summary(res_TDM1_Her2)
res_TDM1=as.data.frame(res_TDM1_Her2)
res_TDM1$gene=row.names(res_TDM1)

########for DHP HER2-enriched########
res_DHP_Her2 <- DESeq(DHP_Her2)
res_DHP_Her2 <- lfcShrink(res_DHP_Her2,coef=c("Response_pCR_vs_RD"), type="apeglm")
summary(res_DHP_Her2)
res_DHP=as.data.frame(res_DHP_Her2)
res_DHP$gene=row.names(res_DHP)

driverGenes<- scan("E:/Projects/PREDIX_HER2/Multimodal/Resource/breast-cancer-driver-genes.txt", what=character()) #breast cancer
#cancergenes<-fread('E:/Projects/PREDIX_HER2/Multimodal/Resource/Census_allMon Nov 28 13_48_06 2022.csv')
#cancergenes=cancergenes$`Gene Symbol`
#driverGenes=c(driverGenes,cancergenes)

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

# figS5h i 
library(gprofiler2);library(readr);library(ggplot2)
res_TDM1=openxlsx::read.xlsx("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEG_All_RNAseq.xlsx",sheet=4)
res_TDM1=res_TDM1[!is.na(res_TDM1$padj),]
res_TDM1=res_TDM1[order(-res_TDM1$log2FoldChange),]
res_DHP=openxlsx::read.xlsx("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEG_All_RNAseq.xlsx",sheet=6)
res_DHP=res_DHP[!is.na(res_DHP$padj),]
res_DHP=res_DHP[order(-res_DHP$log2FoldChange),]

DHP_up=res_DHP$gene[res_DHP$log2FoldChange>0.5&res_TDM1$padj<0.05]
TDM1_up=res_TDM1$gene[res_TDM1$log2FoldChange>0.5&res_TDM1$padj<0.05]
overexpression<- gost(query = list("overexpression (DHP)" = DHP_up,"overexpression (T-DM1)" =TDM1_up), 
                      sources = c("GO", "KEGG", "REAC"),
                      multi_query = TRUE)
head(overexpression$result,20)
p=gostplot(overexpression, capped = F, interactive = F)
pp <- publish_gostplot(p, highlight_terms = c("GO:0016020", "GO:0007267","GO:0012505"), 
                       width = NA, height = NA, filename = NULL)

# 10X5
DHP_down=res_DHP$gene[res_DHP$log2FoldChange<(-0.5)&res_TDM1$padj<0.05]
TDM1_down=res_TDM1$gene[res_TDM1$log2FoldChange<(-0.5)&res_TDM1$padj<0.05]
underexpression<- gost(query = list("underexpression (DHP)" = DHP_down,"underexpression (T-DM1)" =TDM1_down), 
                      sources = c("GO", "KEGG", "REAC"),
                      multi_query = TRUE)
head(underexpression$result,20)
p=gostplot(underexpression, capped = F, interactive = F)
pp <- publish_gostplot(p, highlight_terms = c("GO:0071944", "GO:0005886","GO:0023052","GO:0007154"), 
                       width = NA, height = NA, filename = NULL)


#############################$
# validation NCT02326974 data
#############################$
library(tidyverse);library(data.table);library(DESeq2);library(edgeR)
rm(list=ls())
#exp=read.csv("E:/Projects/PREDIX_HER2/Multimodal/Validation/NCT02326974/GSE243375_All_Samples_TMM_log2_CPM.csv",row.names = 1)%>%as.data.frame()
cts=read.csv("E:/Projects/PREDIX_HER2/Multimodal/Validation/NCT02326974/GSE243375_All_Sample_Raw_Counts.csv",row.names = 1)%>%na.omit()%>%as.data.frame()
txi=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/RNAseq/PREDIX_HER2_tximport_symbol.rds")
length=txi$length%>%as.data.frame()
cts=cts[intersect(row.names(cts),row.names(length)),]
#order=data.frame(sampleID=colnames(cts),N_count=rowSums(t(cts)))
#order=order[order(order$N_count,decreasing =F),]
#order$sampleID<- gsub("site_(\\d+)", "site\\1", order$sampleID)
#order$sampleID<- gsub("core_(\\d+)", "core\\1", order$sampleID)
length=length[intersect(row.names(cts),row.names(length)),1]
length=as.data.frame(matrix(rep(length,ncol(cts)), ncol = ncol(cts)))
row.names(length)=row.names(cts);colnames(length)=colnames(cts)
normMat <- length[row.names(cts),colnames(cts)]
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat
library(edgeR)
eff.lib <- calcNormFactors(normCts,method="TMM") * colSums(normCts)
head(eff.lib )
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)
y <- DGEList(cts)
y <- scaleOffset(y, normMat)
rpkms <- edgeR::rpkm(y$counts, gene.length=length[row.names(y$counts),c(1)], log = FALSE)
tpm<- t( t(rpkms) / colSums(rpkms) ) * 1e6
tpm=log2(tpm+1)
colnames(tpm)<- gsub("site_(\\d+)", "site\\1", colnames(tpm))
colnames(tpm) <- gsub("core_(\\d+)", "core\\1", colnames(tpm))
meta=openxlsx::read.xlsx("E:/Projects/PREDIX_HER2/Multimodal/Validation/NCT02326974/176454-JCI-CMED-RV-3_sd_749148.xlsx")
meta=meta%>%filter(tpt=="Pre")
meta$sampleID=paste0("X",meta$sampleID)
#meta=left_join(meta,order,by="sampleID")
#meta=meta[order(meta$N_count,decreasing =F),]
meta=meta[!duplicated(meta$patientID),]
tpm=tpm[,intersect(colnames(tpm),meta$sampleID)]
meta=meta[meta$sampleID%in%intersect(meta$sampleID,colnames(tpm)),]
all.equal(meta$sampleID,colnames(tpm))
table(meta$HR)
############################
#   GGI PIK3CA GENIUS       #
# https://www.bioconductor.org/packages/release/bioc/vignettes/genefu/inst/doc/genefu.html
#############################
library(genefu);library(org.Hs.eg.db)
genes <- rownames(tpm) %>% as.data.frame()
colnames(genes) <- "probe"
genes$EntrezGene.ID <-   mapIds(org.Hs.eg.db,
                                keys=genes$probe,
                                column="ENTREZID",
                                keytype="SYMBOL",
                                multiVals="first")
row.names(genes)=genes$probe
data(sig.ggi);data(sig.pik3cags);data(sig.genius)
#ggi
ggi_sig<- ggi(data=t(tpm), annot=genes, do.mapping=TRUE)
#PI3KCA
pik3ca_sig <- pik3cags(data=t(tpm), annot=genes, do.mapping=TRUE)
genefu=data.frame(ggi_sig=ggi_sig$score,pik3ca_sig=pik3ca_sig)
genefu$sampleID=row.names(genefu)
####ABC_TRANSPORTERS#######
library(IOBR)
input=as.matrix(tpm)%>%round(2)
ABC=calculate_sig_score(eset=input, signature = kegg, method = "ssgsea",parallel.size=16)
ABC=ABC[,c("ID","KEGG_ABC_TRANSPORTERS","KEGG_APOPTOSIS","KEGG_LYSOSOME","KEGG_ENDOCYTOSIS","KEGG_OXIDATIVE_PHOSPHORYLATION",
           "KEGG_PURINE_METABOLISM","KEGG_CITRATE_CYCLE_TCA_CYCLE","KEGG_GLUTATHIONE_METABOLISM","KEGG_FATTY_ACID_METABOLISM")]
colnames(ABC)=c('sampleID',"ABC_transporter","Apoptosis","Lysosome","Endocytosis","Oxidative_phosphorylation",
                "Purine_metabolism","Citrate_cycles","Glutathione_metabolism","Fatty_acid_metabolism")

sig=calculate_sig_score(eset=input, signature = signature_collection, method = "ssgsea",parallel.size=16)
sig=sig[,c("ID","Glycolysis","Nature_metabolism_Hypoxia","EMT2","CMLS_Review_Exosome")]
colnames(sig)=c('sampleID',"Glycolysis","Hypoxia","EMT","Exosome")

Sig= left_join(ABC,genefu,by="sampleID")%>%left_join(sig,by="sampleID")

##################
# HER2 signaling #
##################
#library(tidyverse);library(data.table);library(GSVA);library(GSEABase);library(dplyr)
#geneSet=getGmt("E:/Projects/PREDIX_HER2/Multimodal/Resource/REACTOME_SIGNALING_BY_ERBB2.v2023.2.Hs.gmt")
#HER2=gsva(as.matrix(tpm),geneSet,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)%>%t()%>%scale()%>%as.data.frame() #ssGSEA计算
#HER2$sampleID=row.names(HER2)

HER2DX=data.frame(sampleID=colnames(tpm),
                  HER2DX_IGG=as.numeric(tpm["CD27",]+tpm["CD79A",]+tpm["HLA-C",]+tpm["IGKC",]+
                                          tpm["IGLC6",]+tpm["IGLV3-25",]+tpm["IL2RG",]+tpm["CXCL8",]+tpm["LAX1",]+tpm["NTN3",]+
                                          tpm["PIM2",]+tpm["POU2AF1",]+tpm["TNFRSF17",]),
                  HER2DX_prolif=as.numeric(tpm["EXO1",]+tpm["ASPM",]+tpm["NEK2",]+tpm["KIF23",]),
                  HER2DX_luminal=as.numeric(tpm["BCL2",]+tpm["DNAJC12",]+tpm["AGR3",]+tpm["AFF3",]+tpm["ESR1",]),
                  HER2DX_HER2_amplicon=as.numeric(tpm["ERBB2",]+tpm["GRB7",]+tpm["STARD3",]+tpm["TCAP",]) 
)%>%mutate(HER2DX_pCR_likelihood_score=HER2DX_IGG/14+HER2DX_prolif/4+HER2DX_HER2_amplicon/5-HER2DX_luminal/4,
           HER2DX_risk_score=HER2DX_IGG/14+HER2DX_luminal/4-HER2DX_prolif/4)

Mitotic=as.data.frame(t(tpm[c("BUB1B","CDK1","AURKB","TTK"),]))
Ceramide=as.data.frame(t(tpm[c("UGCG","COL4A3BP"),]))
Taxane=data.frame(sampleID=colnames(tpm),Taxane_response=rowMeans(Mitotic)-rowMeans(Ceramide))

rna=left_join(Sig,HER2DX,by="sampleID")%>%left_join(Taxane,by="sampleID")

library(tableone);library(tidyverse);library(caret);library(foreach);library(stats);library(lmtest);library(data.table);library(readxl)
library(forestplot);library(ggpubr)
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
data=left_join(rna,meta,by="sampleID")%>%as.data.frame()

data$Arm=NA
data$Arm <-regmatches(data$sampleID, regexpr("site\\d+", data$sampleID))
data$Arm[data$Arm=="site1"]="DHP"
data$Arm[data$Arm=="site2"]="T-DM1"
table(data$Arm)
data$pCR[data$pCR=="Yes"]=1
data$pCR[data$pCR=="No"]=0

variable=c("Taxane_response","HER2DX_pCR_likelihood_score","pik3ca_sig","ABC_transporter","Apoptosis","Lysosome","Endocytosis","EMT","Exosome",
           "Oxidative_phosphorylation","Purine_metabolism","Citrate_cycles","Glutathione_metabolism","Fatty_acid_metabolism",
           "Glycolysis","Hypoxia")
norm_variable=c("Taxane_response","HER2DX_pCR_likelihood_score","pik3ca_sig","ABC_transporter","Apoptosis","Lysosome","Endocytosis","EMT","Exosome",
                "Oxidative_phosphorylation","Purine_metabolism","Citrate_cycles","Glutathione_metabolism","Fatty_acid_metabolism",
                "Glycolysis","Hypoxia")
data[,norm_variable]=scale(data[,norm_variable]) 
data[,variable]

results=Logistic_batch_adjER(as.data.frame(data),"pCR","Arm",variable,"HR")%>%as.data.frame()
TDM1=results[,c("biomarker","whole_OR","whole_lr_p")]
TDM1$group="TDM1"
colnames(TDM1)=c("Signature","OR","Pvalue","group")
df=TDM1
df$Pvalue=as.numeric(df$Pvalue)
df$log10P=-log10(df$Pvalue)
df$logOR=log(as.numeric(df$OR))
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
















