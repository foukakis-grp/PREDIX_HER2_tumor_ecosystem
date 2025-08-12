##################################
###############Fig4a##############
##################################
library(data.table);library(tidyverse);library(ggpubr)
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
clin$patientID=as.character(clin$patientID)
meta=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Meta_adjust_V2.tsv')%>%as.data.frame()
meta$sampleID=meta$Sample_ID
meta$patientID=as.character(meta$patientID)
MS=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Quant_adjust_V2.tsv')%>%as.data.frame()
row.names(MS)=MS$Gene_name
MS$Gene_name=NULL
all.equal(colnames(MS),meta$Sample_ID)
colnames(MS)=meta$patientID
transport=c("SLC12A2","SLC25A10","RAB11B","EEA1","ARL1","FLOT1","VAMP3")
MS=MS[transport,]%>%t()%>%as.data.frame()
MS$patientID=row.names(MS)
df=left_join(MS,clin,by="patientID")
df=df[,c(transport,"Response","Arm")]
d <- reshape2::melt(df,id.vars=c("Response","Arm"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
eFig7a <-
  ggplot(d,aes(x=Arm,y=value,fill=Response))+
  geom_boxplot(outlier.size = 0.5, width=0.7)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Response")+
  stat_compare_means(aes(group=Response),label = "p.format", hide.ns = F,size=5,
                     color="black", label.y.npc = 0.97)+
  labs(y="Normalized protein abundance",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))
eFig7a 

##16X5
###############################
############Fig4b##############
###############################
# MS Description
library(data.table);library(tidyverse);library(readxl)
MS=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Quant_adjust_HarmonizR.tsv")
N_NA=colSums(is.na(MS[,2:ncol(MS)]))
quantile(N_NA)
#load mRNA-protein correlation results
x=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/mRNA_protein_cor.xlsx")
res_TDM1=openxlsx::read.xlsx("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEG_All_RNAseq.xlsx",sheet=4)
gene1=res_TDM1$gene[!is.na(res_TDM1$padj)&res_TDM1$padj<0.05&abs(res_TDM1$log2FoldChange)>0.5]
res_DHP=openxlsx::read.xlsx("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEG_All_RNAseq.xlsx",sheet=6)
gene2=res_DHP$gene[!is.na(res_DHP$padj)&res_DHP$padj<0.05&abs(res_DHP$log2FoldChange)>0.5]
row.names(x)=x$gene
length(intersect(x$gene,c(gene1,gene2)))
x=x[ intersect(x$gene,c(gene1,gene2)),]
fig4b=x
write.table(fig4b,file="E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure4/fig4b.txt",quote = F,row.names =F,sep="\t")
library(ggplot2)
br <- c(seq(min(x$r,na.rm = TRUE)-0.1,-0.000000001,by=0.05),
        seq(0,max(x$r,na.rm = TRUE),by=0.05))
max_count <- max(table(cut(x$r,breaks = br)))
positive_cor_ratio <- sum(x$r>0,na.rm = TRUE)/nrow(x)

sig_positive_cor_ratio <- sum(x$FDR <= 0.05,na.rm = TRUE)/nrow(x)

mean_cor <- mean(x$r,na.rm = TRUE)
median_cor <- median(x$r,na.rm = TRUE)
n_pairs <- nrow(x)

positive_cor_ratio <- sprintf("%.2f%%",100*positive_cor_ratio)
sig_positive_cor_ratio <- sprintf("%.2f%%",100*sig_positive_cor_ratio) #%>% paste("*\\'%\\'~",sep = "")

x %>% mutate(col=ifelse(r>0,"#BC3C29FF","#2D6DB1")) %>%
  ggplot(aes(x=r,fill=col)) +
  geom_histogram(breaks=br,color="black",size=0.1)+
  xlab("Spearman Correlation Coefficient")+
  ylab("Frequency")+
  ggpubr::theme_pubr()+
  theme(legend.position = "none")+
  scale_fill_manual(breaks = c("#BC3C29FF", "#2D6DB1"),
                    values=c("#BC3C29FF", "#2D6DB1"))+
  geom_vline(xintercept = mean_cor,linetype=2)+
  annotate("text",x=mean_cor,y=1.05*max_count,label=paste("Median = ",sprintf("%.2f",median_cor),sep = ""),
           hjust=-0.05,size=3.5)+
  annotate("text",x=min(br)+0.01,y=1.05*max_count,
           label=paste0(n_pairs," pairs\n",
                        positive_cor_ratio," positive correlation\n",
                        sig_positive_cor_ratio," significant\npositive correlation\n(adjusted P <= 0.05)"),
           vjust=1,hjust=0,size=3.5)
#4X4.5 50% P

##################################
###############Fig5a##############
##################################
require(openxlsx);library(ggplot2);library(data.table);library(tidyverse);library(glmmSeq)
driverGenes<- scan("E:/Projects/PREDIX_HER2/Multimodal/Resource/breast-cancer-driver-genes.txt", what=character()) #breast cancer
res_TDM1=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/DEqMS_output_TDM1_pCR_RD_Covar_ER.tsv")%>%as.data.frame()
res_DHP=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/DEqMS_output_DHP_pCR_RD_Covar_ER.tsv")%>%as.data.frame()
row.names(res_TDM1)=res_TDM1$gene
row.names(res_DHP)=res_DHP$gene
res_DHP=res_DHP[row.names(res_TDM1),]

genes_to_showname =intersect(union(res_TDM1$gene[abs(res_TDM1$logFC)>0.5&res_TDM1$P.Value<0.05],
                                   res_DHP$gene[abs(res_DHP$logFC)>0.5&res_DHP$P.Value<0.05]),driverGenes) 
a=res_TDM1[res_TDM1$P.Value<0.05,]
genes_to_showname=union(genes_to_showname,union(a$gene[order(a$logFC)][1:2],
                                                a$gene[order(-a$logFC)][1:2]))
b=res_DHP[res_DHP$P.Value<0.05,]
genes_to_showname=union(genes_to_showname,union(b$gene[order(b$logFC)][1:2],
                                                b$gene[order(-b$logFC)][1:2]))
all.equal(res_TDM1$gene,res_DHP$gene)
data=data.frame(ID=res_TDM1$gene,x=res_DHP$logFC,x_q=res_DHP$P.Value,y=res_TDM1$logFC,y_q=res_TDM1$P.Value)
genes_to_showname
intersect(genes_to_showname,data$ID)

genes_to_showname=union(genes_to_showname,c("GRB7","HLA-B","IL18","MAST4","FBN1","COX11",
                                            "ZNF652","RAB5C","RAB11B","ATP2C1"))
genes_to_showname=intersect(genes_to_showname,data$ID)
data$combined_sig=0
data$combined_sig[data$x_q<0.05&data$y_q>0.05&abs(data$x)>0.5]=1 # sig only in one arm
data$combined_sig[data$x_q>0.05&data$y_q<0.05&abs(data$y)>0.5]=1 # sig only in one arm
data$combined_sig[data$x_q<0.05&data$y_q<0.05&abs(data$x)>0.5&abs(data$y)>0.5]=2 # sig for both arm
table(data$combined_sig)
data$group[data$combined_sig==0]="Notsig"
data$group[data$x_q<0.05&data$combined_sig==1]="unique DEG in DHP arm"
data$group[data$y_q<0.05&data$combined_sig==1]="unique DEG in T-DM1 arm"
data$group[data$combined_sig==2]="shared DEG"
fig5a=data
genes_to_showname=union(genes_to_showname,data[data$group=="shared DEG","ID"])

table(data$group)
range(data$x);range(data$y)
data$gene_label <- ""
row.names(data)=data$ID
data[genes_to_showname, ]$gene_label <- genes_to_showname
colours = c('grey', 'green3', 'gold3', 'blue')

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
# 6.5X5 50%

##################################
############Figure4A.1##############
##################################
# MS Description
library(data.table);library(tidyverse);library(readxl)
MS=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Quant_adjust_HarmonizR.tsv")
N_NA=colSums(is.na(MS[,2:ncol(MS)]))
quantile(N_NA)
#load mRNA-protein correlation results
x=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/mRNA_protein_cor.xlsx")
fig5b=x
library(ggplot2)
br <- c(seq(min(x$r,na.rm = TRUE)-0.1,-0.000000001,by=0.01),
        seq(0,max(x$r,na.rm = TRUE),by=0.01))
max_count <- max(table(cut(x$r,breaks = br)))
positive_cor_ratio <- sum(x$r>0,na.rm = TRUE)/nrow(x)

sig_positive_cor_ratio <- sum(x$FDR <= 0.05,na.rm = TRUE)/nrow(x)

mean_cor <- mean(x$r,na.rm = TRUE)
median_cor <- median(x$r,na.rm = TRUE)
n_pairs <- nrow(x)

positive_cor_ratio <- sprintf("%.2f%%",100*positive_cor_ratio)
sig_positive_cor_ratio <- sprintf("%.2f%%",100*sig_positive_cor_ratio) #%>% paste("*\\'%\\'~",sep = "")

x %>% mutate(col=ifelse(r>0,"#BC3C29FF","#2D6DB1")) %>%
  ggplot(aes(x=r,fill=col)) +
  geom_histogram(breaks=br,color="black",size=0.1)+
  xlab("Spearman Correlation Coefficient")+
  ylab("Frequency")+
  ggpubr::theme_pubr()+
  theme(legend.position = "none")+
  scale_fill_manual(breaks = c("#BC3C29FF", "#2D6DB1"),
                    values=c("#BC3C29FF", "#2D6DB1"))+
  geom_vline(xintercept = mean_cor,linetype=2)+
  annotate("text",x=mean_cor,y=1.05*max_count,label=paste("Median = ",sprintf("%.2f",median_cor),sep = ""),
           hjust=-0.05,size=3.5)+
  annotate("text",x=min(br)+0.01,y=1.05*max_count,
           label=paste0(n_pairs," pairs\n",
                        positive_cor_ratio," positive correlation\n",
                        sig_positive_cor_ratio," significant\npositive correlation\n(adjusted P <= 0.05)"),
           vjust=1,hjust=0,size=3.5)
#6X4 50%
##################################
############Figure4A.2##############
##################################
library(clusterProfiler);library(readxl);library(tidyverse);library(limma)
#Hallmark<- gmtPathways("E:/Projects/PREDIX_HER2/Multimodal/Resource/h.all.v2023.2.Hs.symbols.gmt")
x=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/mRNA_protein_cor.xlsx")
eg = bitr(x$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
colnames(eg)[1]="gene"
x=left_join(x,eg,by="gene")
x=x[!is.na(x$ENTREZID),]
# fgsea - get Normalised Enrichment Score
#x=x%>%mutate(Signed_Log10_P_Value = ifelse(r > 0,
#                                     -log10(p),
#                                     log10(p)))
#x=x%>%mutate(Signed_Log10_P_Value =-log10(p))
#x$Signed_Log10_P_Value[x$Signed_Log10_P_Value=="Inf"]=16
median_cor=median(x$r)
df <- x[order(x$r, decreasing = T),]
df$group[df$r<median_cor]="Low"
df$group[df$r>=median_cor]="High" 
df_low_cor=df[df$group=="Low",]
df_high_cor=df[df$group=="High",]
### low correlation
ranks <-df_low_cor$r
names(ranks) <- df_low_cor$ENTREZID

kk2 <- gseKEGG(geneList     = ranks,
               organism     = 'hsa',
               minGSSize    = 20,
               pvalueCutoff = 1,
               pAdjustMethod = "BH",
               verbose      = FALSE)
kegg_low=kk2@result
#kegg_low=kegg_low[kegg_low$Description%in%c("NF-kappa B signaling pathway","Biosynthesis of cofactors","Basal transcription factors","Lysosome","TNF signaling pathway","ErbB signaling pathway","Ribosome"),]
### high correlation
ranks <- df_high_cor$r
names(ranks) <- df_high_cor$ENTREZID

kk2 <- gseKEGG(geneList     = ranks,
               organism     = 'hsa',
               minGSSize    = 20,
               pvalueCutoff = 1,
               pAdjustMethod = "BH",
               verbose      = FALSE)
kegg_high=kk2@result

#kegg_high=kegg_high[kegg_high$Description%in%c("Arginine and proline metabolism","PPAR signaling pathway","Th1 and Th2 cell differentiation",
#                                               "Natural killer cell mediated cytotoxicity","Th17 cell differentiation","Tryptophan metabolism"),]
#######visualize the GSEA results###########
library(ComplexHeatmap)
library(circlize)  # 用于colorRamp2函数
df <- df[order(df$r, decreasing = F),]
#tab <- getGeneKEGGLinks(species="hsa")
#tab$kegg=gsub("path:","",tab$PathwayID)
## enriched in high cor group
df$Arginine_and_proline_metabolism=0;df$Arginine_and_proline_metabolism[df$r>0&df$group=="High"&df$ENTREZID%in%kk2@geneSets$hsa00330]=1;df$Arginine_and_proline_metabolism[df$r<0&df$group=="High"&df$ENTREZID%in%kk2@geneSets$hsa00330]=-1
df$PPAR_signaling_pathway=0;df$PPAR_signaling_pathway[df$r>0&df$group=="High"&df$ENTREZID%in%kk2@geneSets$hsa04650]=1;df$PPAR_signaling_pathway[df$r<0&df$group=="High"&df$ENTREZID%in%kk2@geneSets$hsa04650]=-1
df$Natural_killer_cell_mediated_cytotoxicity=0;df$Natural_killer_cell_mediated_cytotoxicity[df$r>0&df$group=="High"&df$ENTREZID%in%kk2@geneSets$hsa03320]=1;df$Natural_killer_cell_mediated_cytotoxicity[df$r<0&df$group=="High"&df$ENTREZID%in%kk2@geneSets$hsa03320]=-1

#hsa00330                            Arginine and proline metabolism      23       0.5769132  2.311003 0.0000169738
#hsa04650                  Natural killer cell mediated cytotoxicity      33       0.4636195  2.060593 0.0003356842
#hsa03320                                     PPAR signaling pathway      34       0.4447419  1.998715 0.0014225730
#df$Ribosome=0;df$Ribosome[df$r>0&df$ENTREZID%in%unlist(strsplit(kegg_high$core_enrichment[kegg_high$Description=="Ribosome"], "/"))]=1;df$Ribosome[df$r<0&df$ENTREZID%in%unlist(strsplit(kegg_high$core_enrichment[kegg_high$Description=="Ribosome"], "/"))]=-1

## enriched in low cor group
df$Lysosome=0;df$Lysosome[df$r>0&df$group=="Low"&df$ENTREZID%in%kk2@geneSets$hsa04142]=1;df$Lysosome[df$r<0&df$group=="Low"&df$ENTREZID%in%kk2@geneSets$hsa04142]=-1
df$Ribosome=0;df$Ribosome[df$r>0&df$group=="Low"&df$ENTREZID%in%kk2@geneSets$hsa03010]=1;df$Ribosome[df$r<0&df$group=="Low"&df$ENTREZID%in%kk2@geneSets$hsa03010]=-1
df$Spliceosome=0;df$Spliceosome[df$r>0&df$group=="Low"&df$ENTREZID%in%kk2@geneSets$hsa03040]=1;df$Spliceosome[df$r<0&df$group=="Low"&df$ENTREZID%in%kk2@geneSets$hsa03040]=-1

# hsa03010 hsa03010                                                   Ribosome      21      -0.3380910 -2.083179
# hsa03040 hsa03040                                                Spliceosome
# hsa04142 Lysosome

## high cor group
data <- df[,c("Arginine_and_proline_metabolism","PPAR_signaling_pathway","Natural_killer_cell_mediated_cytotoxicity")]%>%t()%>%as.matrix()

# 定义基因名
gene_names <-gsub("_"," ",c("Arginine_and_proline_metabolism","PPAR_signaling_pathway","Natural_killer_cell_mediated_cytotoxicity")) 

# 将0和1映射到颜色
col_fun <- colorRamp2(c(0, 1, -1), c("white","#BC3C29FF","#2D6DB1"))

# 创建热图对象
ht <- Heatmap(data, 
              name = "ht",  # 为热图对象命名
              col = col_fun,
              cluster_rows = FALSE, 
              cluster_columns = FALSE,
              show_heatmap_legend = FALSE,
              row_labels = gene_names,  # 添加行名称
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (data[i, j] == 1) {
                  grid.rect(x, y, width, height, 
                            gp = gpar(col = "#BC3C29FF", fill = "#BC3C29FF", lwd = 0.1))
                } else if (data[i, j] == -1) {
                  grid.rect(x, y, width, height, 
                            gp = gpar(col = "#2D6DB1", fill = "#2D6DB1", lwd = 0.1))
                }
              })

# 绘制热图并添加黑色框
draw(ht)
ht1=decorate_heatmap_body("ht", {
  grid.rect(gp = gpar(col = "black", fill = NA, lwd = 1))
})
# 10X2

## low cor group
#data <- df[,c("NF_kappa_B_signaling_pathway","Metabolic_pathways","Biosynthesis_of_cofactors","Basal_transcription_factors","Lysosome","TNF_signaling_pathway","ErbB_signaling_pathway")]%>%t()%>%as.matrix()
data <- df[,c("Spliceosome","Ribosome","Lysosome")]%>%t()%>%as.matrix()

# 定义基因名
gene_names <- gsub("_"," ",c("Spliceosome","Ribosome","Lysosome"))

# 将0和1映射到颜色
col_fun <- colorRamp2(c(0, 1, -1), c("white","#BC3C29FF","#2D6DB1"))

# 创建热图对象
ht <- Heatmap(data, 
              name = "ht",  # 为热图对象命名
              col = col_fun,
              cluster_rows = FALSE, 
              cluster_columns = FALSE,
              show_heatmap_legend = FALSE,
              row_labels = gene_names,  # 添加行名称
              cell_fun = function(j, i, x, y, width, height, fill) {
                if (data[i, j] == 1) {
                  grid.rect(x, y, width, height, 
                            gp = gpar(col = "#BC3C29FF", fill = "#BC3C29FF",lwd = 0.1))
                } else if (data[i, j] == -1) {
                  grid.rect(x, y, width, height, 
                            gp = gpar(col = "#2D6DB1", fill = "#2D6DB1",lwd = 0.1))
                }
              })

# 绘制热图并添加黑色框
draw(ht)
ht2=decorate_heatmap_body("ht", {
  grid.rect(gp = gpar(col = "black", fill = NA, lwd = 1))
})

# 10X2  35%

###############################################
#Veen correlation between cna-rna cna-protein#
###############################################
library(readxl);library(data.table);library(tidyverse)
seg_count=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/sample_seg_counts.txt")
clinical_anno=read_excel("E:/Projects/PREDIX_HER2/Clin/data_tidy/PREDIX_HER2_23Jan2023.xlsx")%>%dplyr::filter(TREAT%in%c("Experimentell","Standard"))%>%mutate(patientID=as.character(patientID))
seg_count=seg_count[order(seg_count$segment_count),]
seg_count=seg_count[!duplicated(substr(seg_count$sample,1,4)),]
seg_count=seg_count[substr(seg_count$sample,1,4)%in%clinical_anno$patientID,]
sampleid=seg_count$sample

AMP=c("1p31.3","1q32.2","3p11.1","5p11","6q21","7q11.21","8p11.23","8q23.3","10q22.3","11q13.3",
      "12p11.1","12q15","17q12","19q11","20q13.2")
#DEL=c("1p36.32","1p21.1","1p12","1q21.1","1q21.1","2p11.2","2q13","3p21.2","4p16.1","4q12","4q35.2",  
#      "5q13.2","5q35.2","6q27","7p11.2","7q11.21","7q35","8p23.1","8p21.2","8q24.3","9p24.3","9q11",    
#      "9q13","9q21.11","10q11.22","10q11.22","10q23.2","11p15.5","11p15.4","11q14.3","11q25","12q24.33",
#      "13q11","13q14.13","13q34","14q11.2","14q32.11","15q25.2","16p12.3","16p11.2","16q22.2","17p13.1",
#      "17p11.2","17q21.31","18q23","19p13.3","20q11.1","21p11.2","21q22.12","22p11.2","22q13.33"
#)

gene=as.data.frame(fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/all_thresholded.by_genes.txt"))
amp=gene[gene$Cytoband%in%AMP,c("Gene Symbol","Locus ID","Cytoband",sampleid)]
amp_gene=amp$`Gene Symbol`
CNA_mRNA=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/CNA_mRNA_protein_cor.xlsx",sheet=1)
CNA_protein=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/CNA_mRNA_protein_cor.xlsx",sheet=2)
#Located in amplification focal peaks 
amp_gene=intersect(CNA_mRNA$gene,amp_gene)%>%unique()
CNA_mRNA=CNA_mRNA$gene[CNA_mRNA$gene%in%amp_gene&CNA_mRNA$cis_group=="cis_cor"]
CNA_protein=CNA_protein$gene[CNA_protein$gene%in%amp_gene&CNA_protein$cis_group=="cis_cor"]

x=list(CNA_mRNA_correlation=CNA_mRNA,
       CNA_protein_correlation=CNA_protein) 

library(ggplot2)
library(ggVennDiagram)
ggVennDiagram(x) + scale_fill_gradient(low="#0073C2FF",high = "#CD534CFF")

###############################################
#               CNA drivers (AMP)             #
###############################################
amp=amp[amp$`Gene Symbol`%in%intersect(CNA_mRNA,CNA_protein),]
# select gain or amp
index <- apply(amp[,4:ncol(amp)], 1, function(row) sum(row > 1))>10
amp=amp[index,] 
amp$AMP_rate=apply(amp[,4:ncol(amp)], 1, function(row) sum(row > 1))/179
CNA_mRNA=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/CNA_mRNA_protein_cor.xlsx",sheet=1)
row.names(CNA_mRNA)=CNA_mRNA$gene
CNA_mRNA=CNA_mRNA[amp$`Gene Symbol`,]
CNA_protein=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/CNA_mRNA_protein_cor.xlsx",sheet=2)
row.names(CNA_protein)=CNA_protein$gene
CNA_protein=CNA_protein[amp$`Gene Symbol`,]
#### add correlation ####
amp$CNA_RNA_cor=CNA_mRNA$r
amp$CNA_protein_cor=CNA_protein$r
#### add lnOR #####
seg_count=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/sample_seg_counts.txt")
seg_count=seg_count[order(seg_count$segment_count),]
seg_count=seg_count[!duplicated(substr(seg_count$sample,1,4)),]
seg_count=seg_count[substr(seg_count$sample,1,4)%in%clinical_anno$patientID,]
sampleid=seg_count$sample
#gene-level
CNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/all_data_by_genes.txt")%>%as.data.frame()
CNA=CNA[CNA$`Gene Symbol`%in%amp$`Gene Symbol`,c("Gene Symbol",sampleid)]
row.names(CNA)=CNA$`Gene Symbol`
CNA$`Gene Symbol`=NULL
CNA=as.data.frame(t(CNA))
CNA$patientID=substr(row.names(CNA),1,4)%>%as.integer() 
clin=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.txt")
data=left_join(CNA,clin,by='patientID')
# input for 
variable=colnames(CNA)[1:(ncol(CNA)-1)]
##Function without scale
source("E:/Projects/PREDIX_HER2/Multimodal/Code/Logistic_batch.R")
CNA_logi=Logistic_batch_adjER(data,"pCR","Arm",variable,"ER")%>%as.data.frame()
CNA_logi[,2:ncol(CNA_logi)]=apply(CNA_logi[,2:ncol(CNA_logi)],2,as.numeric)
CNA_logi$DHP_OR=CNA_logi$DHP_OR%>%log(base=exp(1))%>%round(digits = 2)
CNA_logi$TDM1_OR=CNA_logi$TDM1_OR%>%log(base=exp(1))%>%round(digits = 2)
#### add logFC RNA and protein #####
library(readxl)
RNA_FC_TDM1=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEG_All_RNAseq.xlsx",sheet = 4)%>%as.data.frame()
RNA_FC_DHP=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEG_All_RNAseq.xlsx",sheet = 6)%>%as.data.frame()
Protein_FC_TDM1=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/DEqMS_output_TDM1_pCR_RD_Covar_ER.tsv")%>%as.data.frame()
Protein_FC_DHP=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/DEqMS_output_DHP_pCR_RD_Covar_ER.tsv")%>%as.data.frame()
row.names(RNA_FC_TDM1)=RNA_FC_TDM1$gene;row.names(RNA_FC_DHP)=RNA_FC_DHP$gene
row.names(Protein_FC_TDM1)=Protein_FC_TDM1$gene;row.names(Protein_FC_DHP)=Protein_FC_DHP$gene
RNA_FC_TDM1=RNA_FC_TDM1[CNA_logi$biomarker,];RNA_FC_DHP=RNA_FC_DHP[CNA_logi$biomarker,];Protein_FC_TDM1=Protein_FC_TDM1[CNA_logi$biomarker,];Protein_FC_DHP=Protein_FC_DHP[CNA_logi$biomarker,]
## visualization ##
library(tidyverse)
library(ggh4x)
library(ggnewscale)
library(ggfun)
df=amp[,c("Gene Symbol","Cytoband","AMP_rate","CNA_RNA_cor","CNA_protein_cor")]
df$AMP_rate=round(100*df$AMP_rate,digits =0)
df$CNA_RNA_cor=round(df$CNA_RNA_cor,digits =2)
df$CNA_protein_cor=round(df$CNA_protein_cor,digits =2)
all.equal(df$`Gene Symbol`,CNA_logi$biomarker,RNA_FC_TDM1$gene,RNA_FC_DHP$gene,Protein_FC_TDM1$gene,Protein_FC_DHP$gene)
fig5d=df
df$CNA_DHP_lnOR=CNA_logi$DHP_OR%>%round(digits = 2);df$CNA_DHP_p=CNA_logi$DHP_lr_p%>%as.numeric()%>%round(digits = 2);
df$CNA_TDM1_lnOR=CNA_logi$TDM1_OR%>%round(digits = 2);df$CNA_TDM1_p=CNA_logi$TDM1_lr_p%>%as.numeric()%>%round(digits = 2);df$CNA_interact_p=CNA_logi$interaction_lr_p%>%as.numeric()%>%round(digits = 2);
df$RNA_DHP_logFC=RNA_FC_DHP$log2FoldChange%>%round(digits = 2);df$RNA_DHP_padj=RNA_FC_DHP$padj%>%round(digits = 2)
df$RNA_TDM1_logFC=RNA_FC_TDM1$log2FoldChange%>%round(digits = 2);df$RNA_TDM1_padj=RNA_FC_TDM1$padj%>%round(digits = 2)
df$Prot_DHP_logFC=Protein_FC_DHP$logFC%>%round(digits = 2);df$Prot_DHP_p=Protein_FC_DHP$P.Value%>%round(digits = 2)
df$Prot_TDM1_logFC=Protein_FC_TDM1$logFC%>%round(digits = 2);df$Prot_TDM1_p=Protein_FC_TDM1$P.Value%>%round(digits = 2)
#### transfer data ####
df2 <- df %>% tidyr::pivot_longer(cols = -c("Gene Symbol","Cytoband","CNA_DHP_p","CNA_TDM1_p","CNA_interact_p",
                                            "RNA_DHP_padj","RNA_TDM1_padj","Prot_DHP_p","Prot_TDM1_p"), 
                                  names_to = "Group", values_to = "Value")
df2$`Gene Symbol`=factor(df2$`Gene Symbol`,levels = unique(df2$`Gene Symbol`))
library(ggplot2)
library(cowplot)

# 数据子集
df_amp_rate <- subset(df2, Group == "AMP_rate")
df_amp_rate$`Gene Symbol`=factor(df_amp_rate$`Gene Symbol`,levels = unique(df2$`Gene Symbol`))
df_cna_rna_cor <- subset(df2, Group == "CNA_RNA_cor")
df_cna_rna_cor$`Gene Symbol`=factor(df_cna_rna_cor$`Gene Symbol`,levels = unique(df2$`Gene Symbol`))
df_cna_protein_cor <- subset(df2, Group == "CNA_protein_cor")
df_cna_protein_cor$`Gene Symbol`=factor(df_cna_protein_cor$`Gene Symbol`,levels = unique(df2$`Gene Symbol`))
df_FC=subset(df2, Group %in% c("CNA_DHP_lnOR","CNA_TDM1_lnOR","RNA_DHP_logFC","RNA_TDM1_logFC","Prot_DHP_logFC","Prot_TDM1_logFC"))
df_FC$Group=factor(df_FC$Group,levels = c("Prot_TDM1_logFC","Prot_DHP_logFC","RNA_TDM1_logFC","RNA_DHP_logFC","CNA_TDM1_lnOR","CNA_DHP_lnOR"))
df_FC$`Gene Symbol`=factor(df_FC$`Gene Symbol`,levels = unique(df2$`Gene Symbol`))
df_FC$Label=NA
df_FC$Label[df_FC$Group=="Prot_TDM1_logFC"] <- with(df_FC[df_FC$Group=="Prot_TDM1_logFC",], ifelse(Prot_TDM1_p < 0.01, paste(Value, "\n***", sep=""),
                                                                                                   ifelse(Prot_TDM1_p < 0.05, paste(Value, "\n**", sep=""),
                                                                                                          ifelse(Prot_TDM1_p < 0.1, paste(Value, "\n*", sep=""), as.character(Value)))))
df_FC$Label[df_FC$Group=="Prot_DHP_logFC"] <- with(df_FC[df_FC$Group=="Prot_DHP_logFC",], ifelse(Prot_DHP_p < 0.01, paste(Value, "\n***", sep=""),
                                                                                                 ifelse(Prot_DHP_p < 0.05, paste(Value, "\n**", sep=""),
                                                                                                        ifelse(Prot_DHP_p < 0.1, paste(Value, "\n*", sep=""), as.character(Value)))))
df_FC$Label[df_FC$Group=="RNA_TDM1_logFC"] <- with(df_FC[df_FC$Group=="RNA_TDM1_logFC",], ifelse(RNA_TDM1_padj < 0.01, paste(Value, "\n***", sep=""),
                                                                                                 ifelse(RNA_TDM1_padj < 0.05, paste(Value, "\n**", sep=""),
                                                                                                        ifelse(RNA_TDM1_padj < 0.1, paste(Value, "\n*", sep=""), as.character(Value)))))
df_FC$Label[df_FC$Group=="RNA_DHP_logFC"] <- with(df_FC[df_FC$Group=="RNA_DHP_logFC",], ifelse(RNA_DHP_padj < 0.01, paste(Value, "\n***", sep=""),
                                                                                               ifelse(RNA_DHP_padj < 0.05, paste(Value, "\n**", sep=""),
                                                                                                      ifelse(RNA_DHP_padj < 0.1, paste(Value, "\n*", sep=""), as.character(Value)))))
df_FC$Label[df_FC$Group=="CNA_TDM1_lnOR"] <- with(df_FC[df_FC$Group=="CNA_TDM1_lnOR",], ifelse(CNA_TDM1_p < 0.01, paste(Value, "\n***", sep=""),
                                                                                               ifelse(CNA_TDM1_p < 0.05, paste(Value, "\n**", sep=""),
                                                                                                      ifelse(CNA_TDM1_p <= 0.11, paste(Value, "\n*", sep=""), as.character(Value)))))
df_FC$Label[df_FC$Group=="CNA_DHP_lnOR"] <- with(df_FC[df_FC$Group=="CNA_DHP_lnOR",], ifelse(CNA_DHP_p < 0.01, paste(Value, "\n***", sep=""),
                                                                                             ifelse(CNA_DHP_p < 0.05, paste(Value, "\n**", sep=""),
                                                                                                    ifelse(CNA_DHP_p <= 0.11, paste(Value, "\n*", sep=""), as.character(Value)))))
# 绘制每个子集
p1 <- ggplot(df_amp_rate, aes(`Gene Symbol`, Group, fill=Value)) +
  geom_tile(color="grey2", size=0.7) +
  geom_text(aes(label=paste0(Value,"%")), color="black", size=4) +
  scale_fill_gradientn(colours=colorRampPalette(c('#e9e9e9', '#d4d4d4', '#707070'))(111)) +
  labs(x=NULL, y=NULL, fill="AMP rate (%)") +
  theme_minimal() +
  theme(axis.text=element_text(colour='black', size=9)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

p2 <- ggplot(df_cna_rna_cor, aes(`Gene Symbol`, Group, fill=Value)) +
  geom_tile(color="grey2", size=0.7) +
  geom_text(aes(label=Value), color="gray2", size=4) +
  scale_fill_gradientn(colours=colorRampPalette(c('#d4c0d9', '#a982b4', '#643d6e'))(111)) +
  labs(x=NULL, y=NULL, fill="CNA-RNA correlation") +
  theme_minimal() +
  theme(axis.text=element_text(colour='black', size=9)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

p3 <- ggplot(df_cna_protein_cor, aes(`Gene Symbol`, Group, fill=Value)) +
  geom_tile(color="grey2", size=0.7) +
  geom_text(aes(label=Value), color="gray2", size=4) +
  scale_fill_gradientn(colours=colorRampPalette(c('#FFCDB2', '#D06814'))(111)) +
  labs(x=NULL, y=NULL, fill="CNA-protein correlation") +
  theme_minimal() +
  theme(axis.text=element_text(colour='black', size=9)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
library(scales)
p4 <- ggplot(df_FC, aes(`Gene Symbol`, Group, fill=Value)) +
  geom_tile(color="grey2", size=0.7) +
  geom_text(aes(label=Label), color="gray2", size=4) +
  scale_fill_gradientn(
    colours=c('#2D6DB1', 'white', '#DC1623'),
    values=rescale(c(-0.5,  0, 2))
  ) +
  labs(x=NULL, y=NULL, fill="lnOR/logFC") +
  theme_minimal() +
  theme(axis.text=element_text(colour='black', size=9)) +
  theme(axis.text.x=element_text(angle=60, hjust=0.9, vjust=0.9, size=9))

# 提取图例，并设置图例文本和标题的大小
legend1 <- get_legend(p1) 
legend2 <- get_legend(p2)
legend3 <- get_legend(p3)
legend4 <- get_legend(p4)
# 创建一个包含所有图例的行
combined_legend <- plot_grid(legend1, legend2, legend3, legend4, ncol = 4, align = 'h')

# 去掉各个子图的图例
p1 <- p1 + theme(legend.position="none")
p2 <- p2 + theme(legend.position="none")
p3 <- p3 + theme(legend.position="none")
p4 <- p4 + theme(legend.position="none")
# 使用 cowplot 组合图
combined_plot <- plot_grid(combined_legend, 
                           plot_grid(p1, p2, p3, p4, ncol = 1, align = "v", rel_heights = c(1, 1, 1, 4)),
                           ncol = 1, rel_heights = c(0.2, 0.8))

# 绘制图例
combined_plot

# 12X7.5



###############################################
#               Figure5e            #
###############################################
#Pseudo her2+
library(data.table);library(tidyverse);library(ggpubr)
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Meta_adjust_V2.tsv')%>%as.data.frame()
meta$sampleID=meta$Sample_ID
MS=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Quant_adjust_V2.tsv')%>%as.data.frame()
row.names(MS)=MS$Gene_name
MS$Gene_name=NULL
all.equal(colnames(MS),meta$Sample_ID)
colnames(MS)=meta$patientID
MS=MS[c("ERBB2","GRB7","MIEN1"),]%>%t()%>%as.data.frame()
MS$patientID=row.names(MS)
MS$HER2_amplicon=rowSums(MS[,1:3])
pid=intersect(MS$patientID,clin$patientID[clin$ISH_HER2_copy<6&clin$HER2neu1=="2+"])
mean=mean(MS[pid,"HER2_amplicon"])
sd=sd(MS[pid,"HER2_amplicon"])
MS$HER2_low="No"
MS$HER2_low[MS$patientID%in%clin$patientID[clin$ISH_HER2_copy<6&clin$HER2neu1=="2+"]]="Yes"
table(MS$HER2_low)
MS$Zscore=(MS$HER2_amplicon-mean)/sd
MS$HER2_prot=NA
MS$HER2_prot[MS$HER2_low=="Yes"]="HER2_low"
MS$HER2_prot[is.na(MS$HER2_prot)&MS$Zscore<2]="Pseudo_HER2"
MS$HER2_prot[is.na(MS$HER2_prot)]="HER2_pos"
table(MS$HER2_prot)
meta=meta[,c("sampleID","patientID")]
MS$patientID=as.integer(MS$patientID)
MS=left_join(meta,MS,by="patientID")
MS$sampleID=NULL
fig5e=MS
saveRDS(MS,file="E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/ERBB2_MS.rds")
gghistogram(MS, x = "Zscore", bins = 41, rug = TRUE,
            color = "HER2_low", fill = "HER2_low",
            palette = c("darkgreen","darkorchid4")) +
  geom_vline(xintercept = c(-2,0,2,4,6,8), linetype = "dashed", color = "black", size = 0.2)+
  scale_x_continuous(breaks = seq(-2, 10, by = 2)) +
  scale_y_continuous(breaks = seq(2, 10, by = 2)) +
  # Add density plot for HER2_low == "Yes"
  geom_density(data = subset(MS, HER2_low == "Yes"), aes(x = Zscore, y = ..density.. *15), 
               fill = "darkorchid4", alpha = 0.3, size = 0.5) 
# 5X3


blue2red <- colorRampPalette(c("midnightblue","#003BC0","#6689d9","white","#eb8694","#D80D29","firebrick4"))
green2purple <- colorRampPalette(c("darkgreen","white","darkorchid4"))

grey2purple <- colorRampPalette(c("grey80","#8ecee2","#A074B6","darkorchid4"))
####Prot and pCR ####
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
data=left_join(clin,MS,by="patientID")
table(data$Response,data$HER2_prot,data$Arm)
chisq.test(data$Response[data$HER2neu1=="2+"],data$Arm[data$HER2neu1=="2+"])
chisq.test(data$Response[data$HER2neu1=="3+"],data$Arm[data$HER2neu1=="3+"])

###############################################
#               Figure5f            #
###############################################
library(tableone);library(data.table);library(tidyverse)
MS=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/ERBB2_MS.rds")
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
clin$HER2neu1[clin$ISH_HER2_copy<6]
rna=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.rds")
rna$patientID=as.double(rna$patientID)
mut=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.rds")
mut=mut[,c("patientID","coding_mutation_ERBB2")]
mut$patientID=as.double(mut$patientID)
mut$ERBB2_mut[mut$coding_mutation_ERBB2==1]="Mut"
mut$ERBB2_mut[mut$coding_mutation_ERBB2==0]="Wild-type"
mut$coding_mutation_ERBB2=NULL
# CNA CUTseq
seg_count=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/sample_seg_counts.txt")
seg_count=seg_count[order(seg_count$segment_count),]
seg_count=seg_count[!duplicated(substr(seg_count$sample,1,4)),]
seg_count=seg_count[substr(seg_count$sample,1,4)%in%clin$patientID,]
sampleid=seg_count$sample
CNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/CUTseq/GISTIC2/PREDIX_HER2_baseline/all_data_by_genes.txt")%>%as.data.frame()
amp=c("ERBB2")
CNA=CNA[CNA$`Gene Symbol`%in%amp,c("Gene Symbol",sampleid)]
row.names(CNA)=CNA$`Gene Symbol`
CNA$`Gene Symbol`=NULL
CNA=as.data.frame(t(CNA))
CNA$patientID=substr(row.names(CNA),1,4)%>%as.integer() 
colnames(CNA)[1]="ERBB2_CN"
# CNA WES
WES_CNA=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/GISTIC2_WES/all_data_by_genes.txt")%>%as.data.frame()
amp=c("ERBB2")
WES_CNA=WES_CNA[WES_CNA$`Gene Symbol`%in%amp,]
row.names(WES_CNA)=WES_CNA$`Gene Symbol`
WES_CNA$`Gene Symbol`=NULL;WES_CNA$`Gene ID`=NULL;WES_CNA$Cytoband=NULL
WES_CNA=as.data.frame(t(WES_CNA))
WES_CNA$patientID=substr(row.names(WES_CNA),9,12)%>%as.integer() 
colnames(WES_CNA)[1]="ERBB2_CN"
pid=setdiff(WES_CNA$patientID,CNA$patientID)
WES_CNA=WES_CNA[WES_CNA$patientID%in%pid,]
# merge CNA
CNA=rbind(CNA,WES_CNA)
data=left_join(MS,clin,by="patientID")%>%left_join(rna,by="patientID")%>%left_join(CNA,by="patientID")%>%left_join(mut,by="patientID")
data$ERBB2_amp=NA
data$ERBB2_amp[data$ERBB2_CN<2]="No"
data$ERBB2_amp[data$ERBB2_CN>2]="Yes"
table(data$ERBB2_amp)
table(data$ERBB2_amp,data$HER2_prot)
table(data$Response,data$HER2_prot)
table(data$sspbc.subtype,data$HER2_prot)
data$ERBB2_PG="Negative"
data$ERBB2_PG[data$HER2_prot=='HER2_pos'&data$ERBB2_amp=="Yes"]="Positive"
table(data$ERBB2_PG)
data$ERBB2_class=NA
data$ERBB2_class[(data$HER2_prot!='HER2_pos'|data$ERBB2_amp!="Yes")&data$sspbc.subtype=="Her2"]="Her2E ERBB2 PG-"
data$ERBB2_class[(data$HER2_prot!='HER2_pos'|data$ERBB2_amp!="Yes")&data$sspbc.subtype!="Her2"]="Other ERBB2 PG-"
data$ERBB2_class[data$HER2_prot=='HER2_pos'&data$ERBB2_amp=="Yes"&data$sspbc.subtype=="Her2"]="Her2E ERBB2 PG+"
data$ERBB2_class[data$HER2_prot=='HER2_pos'&data$ERBB2_amp=="Yes"&data$sspbc.subtype!="Her2"]="Other ERBB2 PG+"
data$ERBB2_class=factor(data$ERBB2_class,levels = c("Her2E ERBB2 PG-","Other ERBB2 PG-","Her2E ERBB2 PG+","Other ERBB2 PG+"))
table(data$ERBB2_class,data$Response)
table(is.na(data$ERBB2_class))
saveRDS(data,file="E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/ERBB2_PG.rds")
library(data.table);library(tidyverse);library(readxl);library(tidyverse);library(data.table);library(readr);library(readxl)
library(tableone);library(survival);library(tidyverse);library(survminer);library(data.table)
library(ComplexHeatmap);library(circlize);library(RColorBrewer);library(ggstatsplot)
df=data[!is.na(data$ERBB2_class),]
df=df[order(df$ERBB2_class),]
row.names(df)=df$patientID
df=df[,c("ERBB2_class","Arm","Response","EFS.status","sspbc.subtype","Her2ADC","ERBB2_PG","ER","HER2neu1","ERBB2_amp","ERBB2_mut","ERBB2","GRB7","MIEN1")]
#### HER2DX group
protein=c("ERBB2","GRB7","MIEN1")
s1 = as.matrix(t(df[df$ERBB2_class=="Her2E ERBB2 PG-",protein]))   # subtype1
S1_meta=df[colnames(s1),]
s2 = as.matrix(t(df[df$ERBB2_class=="Other ERBB2 PG-",protein]))   # subtype1
S2_meta=df[colnames(s2),]
s3 = as.matrix(t(df[df$ERBB2_class=="Her2E ERBB2 PG+",protein]))   # subtype1
S3_meta=df[colnames(s3),]
s4 = as.matrix(t(df[df$ERBB2_class=="Other ERBB2 PG+",protein]))   # subtype1
S4_meta=df[colnames(s4),]

pheno=S1_meta
ha1=HeatmapAnnotation(ERBB2_class=pheno$ERBB2_class,
                      Arm=pheno$Arm,
                      Response=pheno$Response,
                      Event=pheno$EFS.status,
                      Intrinsic_subtype=pheno$sspbc.subtype,
                      Her2ADC=pheno$Her2ADC,
                      ERBB2_PG=pheno$ERBB2_PG,
                      ER=pheno$ER,
                      HER2_IHC=pheno$HER2neu1,
                      ERBB2_amp=pheno$ERBB2_amp,
                      ERBB2_mut=pheno$ERBB2_mut,
                      col=list(ERBB2_class=c("Her2E ERBB2 PG-"="#F6BD60","Other ERBB2 PG-"="#F7EDE2","Her2E ERBB2 PG+"="#F5CAC3","Other ERBB2 PG+"="#84A59D"),
                               Arm=c("DHP"="#8491B4FF","T-DM1"="#91D1C2FF"),
                               Response=c("pCR"="#E8E8E8","RD"="black"),
                               Event=c("0"="#E8E8E8","1"="black"),
                               Intrinsic_subtype= c("LumA"="#1f78b4", "LumB"="#a6cee3", "Her2"="#fb9a99","Basal"="#e31a1c"),
                               Her2ADC=c("S1"="#4DE66699","S2"="#F39B7FFF","S3"="#4D1A6699"),
                               ERBB2_PG=c("Negative"="#E8E8E8","Positive"="black"),
                               ER=c("negative"="#E8E8E8","positive"="black"),
                               HER2_IHC=c("2+"="#E8E8E8","3+"="black"),
                               ERBB2_amp=c("No"="#E8E8E8","Yes"="black"),
                               ERBB2_mut=c("Wild-type"="#E8E8E8","Mut"="black")),
                      height = unit(3, "cm"),na_col="#808080",show_annotation_name = FALSE,
                      border=F)

pheno=S2_meta
ha2=HeatmapAnnotation(ERBB2_class=pheno$ERBB2_class,
                      Arm=pheno$Arm,
                      Response=pheno$Response,
                      Event=pheno$EFS.status,
                      Intrinsic_subtype=pheno$sspbc.subtype,
                      Her2ADC=pheno$Her2ADC,
                      ERBB2_PG=pheno$ERBB2_PG,
                      ER=pheno$ER,
                      HER2_IHC=pheno$HER2neu1,
                      ERBB2_amp=pheno$ERBB2_amp,
                      ERBB2_mut=pheno$ERBB2_mut,
                      col=list(ERBB2_class=c("Her2E ERBB2 PG-"="#F6BD60","Other ERBB2 PG-"="#F7EDE2","Her2E ERBB2 PG+"="#F5CAC3","Other ERBB2 PG+"="#84A59D"),                               Arm=c("DHP"="#8491B4FF","T-DM1"="#91D1C2FF"),
                               Response=c("pCR"="#E8E8E8","RD"="black"),
                               Event=c("0"="#E8E8E8","1"="black"),
                               Intrinsic_subtype= c("LumA"="#1f78b4", "LumB"="#a6cee3", "Her2"="#fb9a99","Basal"="#e31a1c"),
                               Her2ADC=c("S1"="#4DE66699","S2"="#F39B7FFF","S3"="#4D1A6699"),
                               ERBB2_PG=c("Negative"="#E8E8E8","Positive"="black"),
                               ER=c("negative"="#E8E8E8","positive"="black"),
                               HER2_IHC=c("2+"="#E8E8E8","3+"="black"),
                               ERBB2_amp=c("No"="#E8E8E8","Yes"="black"),
                               ERBB2_mut=c("Wild-type"="#E8E8E8","Mut"="black")),
                      height = unit(3, "cm"),na_col="#808080",show_annotation_name = FALSE,
                      border=F)

pheno=S3_meta
ha3=HeatmapAnnotation(ERBB2_class=pheno$ERBB2_class,
                      Arm=pheno$Arm,
                      Response=pheno$Response,
                      Event=pheno$EFS.status,
                      Intrinsic_subtype=pheno$sspbc.subtype,
                      Her2ADC=pheno$Her2ADC,
                      ERBB2_PG=pheno$ERBB2_PG,
                      ER=pheno$ER,
                      HER2_IHC=pheno$HER2neu1,
                      ERBB2_amp=pheno$ERBB2_amp,
                      ERBB2_mut=pheno$ERBB2_mut,
                      col=list(ERBB2_class=c("Her2E ERBB2 PG-"="#F6BD60","Other ERBB2 PG-"="#F7EDE2","Her2E ERBB2 PG+"="#F5CAC3","Other ERBB2 PG+"="#84A59D"),                               Arm=c("DHP"="#8491B4FF","T-DM1"="#91D1C2FF"),
                               Response=c("pCR"="#E8E8E8","RD"="black"),
                               Event=c("0"="#E8E8E8","1"="black"),
                               Intrinsic_subtype= c("LumA"="#1f78b4", "LumB"="#a6cee3", "Her2"="#fb9a99","Basal"="#e31a1c"),
                               Her2ADC=c("S1"="#4DE66699","S2"="#F39B7FFF","S3"="#4D1A6699"),
                               ERBB2_PG=c("Negative"="#E8E8E8","Positive"="black"),
                               ER=c("negative"="#E8E8E8","positive"="black"),
                               HER2_IHC=c("2+"="#E8E8E8","3+"="black"),
                               ERBB2_amp=c("No"="#E8E8E8","Yes"="black"),
                               ERBB2_mut=c("Wild-type"="#E8E8E8","Mut"="black")),
                      height = unit(3, "cm"),na_col="#808080",show_annotation_name = FALSE,
                      border=F)

pheno=S4_meta
ha4=HeatmapAnnotation(ERBB2_class=pheno$ERBB2_class,
                      Arm=pheno$Arm,
                      Response=pheno$Response,
                      Event=pheno$EFS.status,
                      Intrinsic_subtype=pheno$sspbc.subtype,
                      Her2ADC=pheno$Her2ADC,
                      ERBB2_PG=pheno$ERBB2_PG,
                      ER=pheno$ER,
                      HER2_IHC=pheno$HER2neu1,
                      ERBB2_amp=pheno$ERBB2_amp,
                      ERBB2_mut=pheno$ERBB2_mut,
                      col=list(ERBB2_class=c("Her2E ERBB2 PG-"="#F6BD60","Other ERBB2 PG-"="#F7EDE2","Her2E ERBB2 PG+"="#F5CAC3","Other ERBB2 PG+"="#84A59D"),                               Arm=c("DHP"="#8491B4FF","T-DM1"="#91D1C2FF"),
                               Response=c("pCR"="#E8E8E8","RD"="black"),
                               Event=c("0"="#E8E8E8","1"="black"),
                               Intrinsic_subtype= c("LumA"="#1f78b4", "LumB"="#a6cee3", "Her2"="#fb9a99","Basal"="#e31a1c"),
                               Her2ADC=c("S1"="#4DE66699","S2"="#F39B7FFF","S3"="#4D1A6699"),
                               ERBB2_PG=c("Negative"="#E8E8E8","Positive"="black"),
                               ER=c("negative"="#E8E8E8","positive"="black"),
                               HER2_IHC=c("2+"="#E8E8E8","3+"="black"),
                               ERBB2_amp=c("No"="#E8E8E8","Yes"="black"),
                               ERBB2_mut=c("Wild-type"="#E8E8E8","Mut"="black")),
                      height = unit(3, "cm"),na_col="#808080",show_annotation_name = T,
                      border=F)
#### Plot heatmap
col_fun = colorRamp2(c(-3,-2,-1,0,1,2,3), c("midnightblue","#003BC0","#6689d9","white","#eb8694","#D80D29","firebrick4"))

ht_list = Heatmap(as.matrix(s1), col = col_fun, name = "Protein abundance",
                  cluster_rows =T,show_row_names = FALSE,cluster_columns = F,
                  show_row_dend = T, show_column_dend = F,
                  show_column_names = F,
                  top_annotation = ha1,
                  row_title_gp = gpar(col = "#FFFFFF00"), width = unit(2.1, "cm"),height= unit(3, "cm"))+
  Heatmap(as.matrix(s2),col=col_fun, show_column_names = F, cluster_columns = F,
          show_column_dend = F,top_annotation=ha2,
          show_heatmap_legend = FALSE, width = unit(3.2, "cm"),height= unit(3, "cm"))+
  Heatmap(as.matrix(s3),col=col_fun, show_column_names = F, cluster_columns = F,
          show_column_dend = F,top_annotation=ha3,row_names_gp = gpar(fontsize = 2),
          show_heatmap_legend = FALSE, width = unit(5.9, "cm"),height= unit(3, "cm"))+
  Heatmap(as.matrix(s4),col=col_fun, show_column_names = F, cluster_columns = F,
          show_column_dend = F,top_annotation=ha4,row_names_gp = gpar(fontsize = 8),
          show_heatmap_legend = FALSE, width = unit(2, "cm"),height= unit(3, "cm"))
                                 
p0=draw(ht_list, row_title = paste0("protein abundance"),
        annotation_legend_side = "bottom", heatmap_legend_side = "bottom") 
p0



# 8X8
###############################################
#               Figure5g            #
###############################################
library(tidyverse);library(tableone);library(data.table);library(lmtest);library(forestploter)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/ERBB2_PG.rds")
data$ERBB2_class=factor(data$ERBB2_class,levels = c("Other ERBB2 PG-","Her2E ERBB2 PG-","Her2E ERBB2 PG+","Other ERBB2 PG+"))
table(data$ERBB2_PG,data$Arm)
whole<- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = data[data$ERBB2_PG=="Negative",])
ShowRegTable(whole)
whole<- glm(as.numeric(pCR) ~ Arm+ER, family = "binomial", data = data[data$ERBB2_PG=="Positive",])
ShowRegTable(whole)
interaction_1<- glm(as.numeric(pCR) ~ ERBB2_PG+ER+Arm, family = "binomial", data = data)
interaction_2<- glm(as.numeric(pCR) ~ ERBB2_PG+ER+Arm+Arm*ERBB2_PG, family = "binomial", data = data)
lrtest(interaction_1,interaction_2)

df=fread("E:/Projects/PREDIX_HER2/Multimodal/Figures/Figure5/Subgroup_Forestplot.csv")
df$DHP <- ifelse(is.na(df$DHP), "", df$DHP)
df$`T-DM1` <- ifelse(is.na(df$`T-DM1`), "", df$`T-DM1`)
df$`P for interaction` <- ifelse(is.na(df$`P for interaction`), "", df$`P for interaction`)
df$` ` <- paste(rep(" ", 10), collapse = " ")
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

# 6X3

###############################################
#               Figure5h           #
###############################################
library(data.table);library(tidyverse);library(ggpubr)
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
clin$patientID=as.character(clin$patientID)
meta=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Meta_adjust_V2.tsv')%>%as.data.frame()
meta$sampleID=meta$Sample_ID
meta$patientID=as.character(meta$patientID)
MS=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Quant_adjust_V2.tsv')%>%as.data.frame()
row.names(MS)=MS$Gene_name
MS$Gene_name=NULL
all.equal(colnames(MS),meta$Sample_ID)
colnames(MS)=meta$patientID
library(IOBR)
library(GSVA)
input=as.matrix(MS)%>%round(2)
genelist <- list(c("RAB11B","RAB11FIP1","RAB5A","ERLIN2","SLC12A2","SLC25A10","ABCC12"))
names(genelist)="ADC"
ADC=calculate_sig_score_ssgsea(eset = input, signature = genelist)
colnames(ADC)[1]="patientID"
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/ERBB2_PG.rds")
data$patientID=as.character(data$patientID)
data=data[,c("patientID","ERBB2_PG")]

Sig= left_join(meta,ADC,by="patientID")%>%left_join(data,by="patientID")
#Sig=Sig[Sig$Arm=="DHP",]
df=Sig[,c("ADC","ERBB2_PG")]
d <- reshape2::melt(df,id.vars=c("ERBB2_PG"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
eFig7a <-
  ggplot(d,aes(x=ERBB2_PG,y=value,fill=ERBB2_PG))+
  geom_boxplot(outlier.size = 0.5, width=0.7)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Response")+
  stat_compare_means(aes(group=ERBB2_PG),label = "p.format", hide.ns = F,size=5,
                     color="black", label.y.npc = 0.97)+
  labs(y="Normalized protein abundance",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))
eFig7a 





library(data.table);library(tidyverse);library(ggpubr)
clin=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
clin$patientID=as.character(clin$patientID)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/ERBB2_PG.rds")
data$patientID=as.character(data$patientID)
data=data[,c("patientID","ERBB2_PG")]
genomic=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/genomic_metrics_PREDIX_HER2.txt")
genomic$patientID=as.character(genomic$patientID)
meta=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Meta_adjust_V2.tsv')%>%as.data.frame()
meta$sampleID=meta$Sample_ID
meta$patientID=as.character(meta$patientID)
MS=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Quant_adjust_V2.tsv')%>%as.data.frame()
row.names(MS)=MS$Gene_name
MS$Gene_name=NULL
all.equal(colnames(MS),meta$Sample_ID)
colnames(MS)=meta$patientID
transport=c("SLC12A2","SLC25A10","RAB11B","ABCC3","EEA1","ARL1","FLOT1","VAMP3")
MS=MS[transport,]%>%t()%>%as.data.frame()
MS$patientID=row.names(MS)
df=left_join(MS,clin,by="patientID")%>%left_join(data,by="patientID")%>%left_join(genomic,by="patientID")
df=df[,c(transport,"ERBB2_PG")]
d <- reshape2::melt(df,id.vars=c("ERBB2_PG"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
eFig7a <-
  ggplot(d,aes(x=ERBB2_PG,y=value,fill=ERBB2_PG))+
  geom_boxplot(outlier.size = 0.5, width=0.7)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Response")+
  stat_compare_means(aes(group=ERBB2_PG),label = "p.format", hide.ns = F,size=5,
                     color="black", label.y.npc = 0.97)+
  labs(y="Normalized protein abundance",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0,0.2,1), "lines"))
eFig7a 








library(survival);library(survminer)
data$group="Others"
data$group[data$ERBB2_class=="Her2E ERBB2 PG+"]="Her2E ERBB2 PG+"
table(data$group)
fit <- survfit(Surv(EFS.time,EFS.status) ~ group, data = data)
ggsurvplot(fit,data = data, palette = "jama",
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))



data$group="Others"
data$group[data$Arm=="T-DM1"&data$ERBB2_PG=="Negative"]="DHP,PG-negative"
data$HER2_prot
fit <- survfit(Surv(EFS.time,EFS.status) ~ HER2_prot, data = data)
ggsurvplot(fit,data = data, palette = "jama",
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))
table(data$EFS.status)

fit <- survfit(Surv(EFS.time,EFS.status) ~ HER2_prot, data = data)
ggsurvplot(fit,data = data, palette = "jama",
           pval = TRUE,
           break.time.by = 10,
           xlim = c(0, 60))

whole<- glm(as.numeric(pCR) ~ HER2_prot+ER, family = "binomial", data = data[data$HER2_prot%in%c("HER2_pos","Pseudo_HER2"),])
ShowRegTable(whole)




##################################
############Figure4C##############
##################################
library(readxl);library(data.table);library(tidyverse)
gene=c("EGFR","ERBB2","IGF1R","MUC4","PIK3R1","PTEN","PAK1","AKT1","AKT3","MAP3K1","MAP2K4","JUN",
       "IKBKB","CCND1","TP53BP2","ATM","MDM2","MDM4","CDKN1A","CCNB1","CDK4","CDK6","RB1","BRCA1","BRCA2")
RNA_FC_TDM1=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEG_All_RNAseq.xlsx",sheet = 4)
row.names(RNA_FC_TDM1)=RNA_FC_TDM1$gene;RNA_FC_TDM1=RNA_FC_TDM1[gene,c("log2FoldChange",'padj')];colnames(RNA_FC_TDM1)=paste0("RNA_TDM1",colnames(RNA_FC_TDM1))
row.names(RNA_FC_TDM1)=gene
RNA_FC_DHP=read_excel("E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/DEG_All_RNAseq.xlsx",sheet = 6)
row.names(RNA_FC_DHP)=RNA_FC_DHP$gene;RNA_FC_DHP=RNA_FC_DHP[gene,c("log2FoldChange",'padj')];colnames(RNA_FC_DHP)=paste0("RNA_DHP",colnames(RNA_FC_DHP))
row.names(RNA_FC_DHP)=gene
Protein_FC_TDM1=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/DEqMS_output_TDM1_pCR_RD.tsv")%>%as.data.frame()
row.names(Protein_FC_TDM1)=Protein_FC_TDM1$gene;Protein_FC_TDM1=Protein_FC_TDM1[gene,c("logFC","P.Value")];colnames(Protein_FC_TDM1)=paste0("Protein_TDM1",colnames(Protein_FC_TDM1))
row.names(Protein_FC_TDM1)=gene
Protein_FC_DHP=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/DEqMS_output_DHP_pCR_RD.tsv")%>%as.data.frame()
row.names(Protein_FC_DHP)=Protein_FC_DHP$gene;Protein_FC_DHP=Protein_FC_DHP[gene,c("logFC","P.Value")];colnames(Protein_FC_DHP)=paste0("Protein_DHP",colnames(Protein_FC_DHP))
row.names(Protein_FC_DHP)=gene
df=cbind(RNA_FC_DHP,Protein_FC_DHP,RNA_FC_TDM1,Protein_FC_TDM1)%>%select(c("RNA_DHPlog2FoldChange","Protein_DHPlogFC","RNA_TDM1log2FoldChange","Protein_TDM1logFC"))

library(circlize);library(ComplexHeatmap)
col_fun = colorRamp2(c(-1, 0, 1), c("#377EB8", "white", "#E41A1C"))
Heatmap(as.matrix(df), name = "mat", col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.2f", df[i, j]), x, y, gp = gpar(fontsize = 10))
        })
# 4X10 45%


### Pseudo 
library(data.table);library(tidyverse)
meta=fread('E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/PREDIX_Prot_Meta_adjust_V2.tsv')%>%as.data.frame()
meta$Sample=meta$Sample_ID
pseudo=fread("E:/Projects/PREDIX_HER2/Multimodal/Data/Proteomic/Pseudo_cluster_membership.tsv")
data=left_join(meta,pseudo,by="Sample")
table(data$Pseudo_Cluster,data$HER2_low)
table(data$Response[data$Pseudo_Cluster==1&data$HER2_low==0])
library(survival)
cox.test <- coxph(Surv(EFS.time,EFS.status)~ as.numeric(`CXCL12`)+Arm+as.factor(ER)+as.factor(TUMSIZE)+as.factor(ANYNODES), data=meta) ##DII_density_with_supp
ShowRegTable(cox.test)
library(tableone)
data$pCR=0
data$pCR[data$Response=="pCR"]=1
log_model <- glm(pCR~ Pseudo_Cluster, data = data[data$HER2_low=="0"&data$Arm=="DHP",], family = "binomial")
ShowRegTable(log_model)

log_model <- glm(pCR~ Pseudo_Cluster, data = data[data$HER2_low=="0"&data$Arm=="T-DM1",], family = "binomial")
ShowRegTable(log_model)


library(openxlsx)
data_list=list("fig5a"=fig5a%>%as.data.frame(),
               "fig5b"=fig5b%>%as.data.frame(),
               "fig5d"=fig5d%>%as.data.frame(),
               "fig5e"=fig5e%>%as.data.frame())
openxlsx::write.xlsx(data_list,file='E:/Projects/PREDIX_HER2/Multimodal/Figures/SourceData/Figure5.xlsx')


