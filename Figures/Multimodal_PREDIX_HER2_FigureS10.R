##############
# FigureS10A #
##############
library(data.table);library(ggplot2);library(scales);library(RColorBrewer);library(tidyverse)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/Xenium_baselineMeta_cell_state_curated.rds")
pid=fread("E:/Projects/PREDIX_HER2/Multimodal/Data_repository/PREDIX_HER2_multiomics_meta.txt")
pid$patientID=NA
pid$patientID=substr(pid$sampleID.wes,9,12)%>%as.character()
pid$patientID[is.na(pid$patientID)]=substr(pid$sampleID.rna[is.na(pid$patientID)],9,12)%>%as.character()
meta=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta$patientID=as.character(meta$patientID)
data=left_join(data,meta,by="patientID")%>%left_join(pid,by="patientID")
data=data[order(data$Arm),]
table(is.na(data$study_id))
table(data$study_id[data$Arm=="T-DM1"])
table(data$study_id[data$Arm=="DHP"])

meta_summary <- data %>%
  group_by(study_id,patientID) %>%
  summarise(num_cells_detected = n(), .groups = "drop")%>%left_join(meta,by="patientID")
meta_summary$study_id=factor(meta_summary$study_id,levels = c("P103","P108","P115","P119","P133","P135",
                                              "P137","P150","P153","P186","P188","P189","P3","P33","P48","P5","P69","P75","P98",
                                              "P114","P132","P182","P184","P41","P43","P63","P84"))

p1 <- ggplot(meta_summary, aes(x = study_id, y = num_cells_detected)) +
  geom_bar(stat = "identity", position = "dodge", fill = "#3C77AF") +
  labs(title = "Number of cells detected", x = "sample ID", y = "cells detected") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  annotate("text",x="P133", y = 120000, 
           label = paste("Median=",median(meta_summary$num_cells_detected)), 
           size = 4, color = "black")
p1

data=fread("E:/Projects/PREDIX_HER2_longitudinal/data/Xenium/merged_metrics_summary.csv")
data$patientID=substr(data$region_name,1,4)
data$tpt=substr(data$region_name,6,7)
data$tpt[data$tpt%in%c("BL","BO")]="pre"
data=data[data$tpt=="pre",]
pid=fread("E:/Projects/PREDIX_HER2/Multimodal/Data_repository/PREDIX_HER2_multiomics_meta.txt")
pid$patientID=NA
pid$patientID=substr(pid$sampleID.wes,9,12)%>%as.character()
pid$patientID[is.na(pid$patientID)]=substr(pid$sampleID.rna[is.na(pid$patientID)],9,12)%>%as.character()
meta=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta$patientID=as.character(meta$patientID)
data=left_join(data,meta,by="patientID")%>%left_join(pid,by="patientID")
data$study_id=factor(data$study_id,levels = c("P103","P108","P115","P119","P133","P135",
                                                              "P137","P150","P153","P186","P188","P189","P3","P33","P48","P5","P69","P75","P98",
                                                              "P114","P132","P182","P184","P41","P43","P63","P84"))
p2 <- ggplot(data, aes(x = study_id, y = median_transcripts_per_cell)) +
  geom_bar(stat = "identity", position = "dodge", fill = "#3C77AF") +
  labs(title = "Median number of transcripts detected per cell", x = "sample ID", y = "Median number of transcripts per Cell") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  annotate("text",x="P133", y = 500, 
           label = paste("Median (total)=",median(data$median_transcripts_per_cell)), 
           size = 4, color = "black")
p2

p=p1/p2

ggsave(p, file="E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS10/FigS10a_b.pdf", width=4.5, height=5)

################
# FigureS10B-C #
################
library(Seurat);library(scCustomize)
library(data.table);library(ggplot2);library(scales);library(RColorBrewer);library(tidyverse)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/Xenium_baselineMeta_cell_state_curated.rds")
data=data[data$cell_state%in%c("Basal_SC","Her2E_SC","LumA_SC","LumB_SC"),]
data <- data %>%
  group_by(patientID, cell_state) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patientID) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
data <- data %>%
  select(patientID, cell_state, proportion) %>%
  pivot_wider(names_from = cell_state, values_from = proportion, values_fill = 0)
pid=fread("E:/Projects/PREDIX_HER2/Multimodal/Data_repository/PREDIX_HER2_multiomics_meta.txt")
pid$patientID=NA
pid$patientID=substr(pid$sampleID.wes,9,12)%>%as.character()
pid$patientID[is.na(pid$patientID)]=substr(pid$sampleID.rna[is.na(pid$patientID)],9,12)%>%as.character()
meta=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta$patientID=as.character(meta$patientID)
rna=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Curated_metrics/transcriptomic_metrics_PREDIX_HER2.rds")
data=left_join(data,meta,by="patientID")%>%left_join(pid,by="patientID")%>%left_join(rna,by="patientID")
table(data$sspbc.subtype)
data$sspbc=NA
#data$sspbc[data$sspbc.subtype%in%c("LumA","LumB")]="LumAorB"
#data$sspbc[data$sspbc.subtype%in%c("Her2")]="Her2"

d=data%>%select(c("sspbc.subtype","LumA_SC","LumB_SC","Her2E_SC","Basal_SC")) #"Basal_SC",,"LumA_SC","LumB_SC"
d <- reshape2::melt(d,id.vars=c("sspbc.subtype"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
Figs10b <-
  ggplot(d,aes(x=variable,y=value,fill=variable))+
  geom_boxplot(outlier.size = 0.2, width=0.8)+
  facet_wrap(sspbc.subtype~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_manual(values = c(
    "LumAorB" = "#1f78b4",
    "Her2" = "#fb9a99"
  ))+
  stat_compare_means(aes(group=variable),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.9)+
  labs(y="tumor cell proportion",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.2,0.2,0.5), "lines"))
Figs10b


library(data.table);library(ggplot2);library(scales);library(RColorBrewer);library(tidyverse)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/Xenium_baselineMeta_cell_state_curated.rds")
table(data$cell_type_major_curated,data$sampleID)
tumor=data[data$cell_state%in%c("LumA_SC","LumB_SC","Her2E_SC","Basal_SC"),]
data$ADC_trafficking1
table(tumor$cell_state)
table(tumor$cell_state,tumor$sampleID)
tumor_cell_proportion <- tumor %>%
  group_by(patientID, cell_state) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patientID) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
tumor_cell_proportion_wide <- tumor_cell_proportion %>%
  select(patientID, cell_state, proportion) %>%
  pivot_wider(names_from = cell_state, values_from = proportion, values_fill = 0)
# read meta
meta=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta$patientID=as.character(meta$patientID)
tumor_cell_proportion_wide=left_join(tumor_cell_proportion_wide,meta,by="patientID")
#colnames(tumor_cell_proportion_wide)=gsub(" cancer (Cell|cell|cells)", "", colnames(tumor_cell_proportion_wide))
d=tumor_cell_proportion_wide%>%select(c("Arm","Response","LumA_SC","LumB_SC","Her2E_SC","Basal_SC")) #
d <- reshape2::melt(d,id.vars=c("Arm","Response"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
Figs10c <-
  ggplot(d,aes(x=Arm,y=value,fill=Response))+
  geom_boxplot(outlier.size = 0.2, width=0.8)+
  facet_wrap(variable~.,scales="free_y",nrow=1)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Treatment Response")+
  stat_compare_means(aes(group=Response),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.9)+
  labs(y="tumor cell proportion",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.2,0.2,0.5), "lines"))
Figs10c

p=Figs10b/Figs10c
ggsave(p,file="E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS10/FigS10bc.pdf", width=7, height=8)

################
# FigureS10d #
################
library(Seurat);library(scCustomize)
options(future.globals.maxSize = 50 * 1024^3)  
object=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/Seurat_cell_state_curated_ondisk.rds")
object$cell_state[object$cell_state%in%c("LumA_SC","LumB_SC","Her2E_SC","Basal_SC")]="Tumor epithelial"
Idents(object)=object$cell_state
object=subset(object,idents = c("Adipocytes","Mast cells","Plasma"), invert = TRUE)
table(Idents(object))
Idents(object)=factor(Idents(object),levels =c("Tumor epithelial","Luminal progenitor","Myoepithelial",
          "Cycling T","Exhausted T","CTLs","CD8 Tem","Th1","Th2","Th17","Treg","NK",
          "Memory B","Naive B","Cycling myeloid","Monocyte","cDC","pDC","LAM","M1 macrophages","M2 macrophages","MDSC","Neutrophils",
          "Cycling Endothelial","Lymphatic Endothelial","Vascular Endothelial","EMT-like CAF","iCAF","myCAF","PVLs") )
Epi=subset(object,idents = c("Tumor epithelial","Luminal progenitor","Myoepithelial"))
Tcell=subset(object,idents = c("Cycling T","Exhausted T","CTLs","CD8 Tem","Th1","Th2","Th17","Treg","NK"))
Bcell=subset(object,idents = c("Memory B","Naive B"))
Myeloid=subset(object,idents = c("Cycling myeloid","Monocyte","cDC","pDC","LAM","M1 macrophages","M2 macrophages","MDSC","Neutrophils"))
Endo=subset(object,idents = c("Cycling Endothelial","Lymphatic Endothelial","Vascular Endothelial"))
Mesenchymal=subset(object,idents = c("EMT-like CAF","iCAF","myCAF","PVLs"))
table(object$cell_state[object$cell_type_major_curated=="T cells"])
#ACTA2 TP63 Myoepithelial
#KIT KRT19 Luminal progenitor 
# T cell res0.5
# Exhausted T cell PDCD1 
# Cycling T cell MKI67 

# CD8T CD8A 
#### Naive CD8 T cell, CCR7 SELL XXX
# CTLs GZMK 
# CD8 Tem CCR7,IL7R 

# CD4T 
# CD4 naive SELL CCR7 IL7R  
#Th1 IFNG CXCR3 CXCL9 CCL20
#Th2 GATA3,IL4R,STAT6,CCL17
#Th17 RORC 
#Treg FOXP3 

# NK cell KLRD1,NCAM1

# Cycling myeloid MKI67
# Monocyte CD14 CD16 (FCGR3A)   
# M1 macrophages CD68 CD86, CD80   
# M2 macrophages CD163,MRC1
# LAM TREM2,FABP4 
# cDC CD1C, ITGAX, CLEC9A 
# pDC CLEC4C
# MDSC CD16b CD33 

# PVLs MCAM 
# EMT-like CAF 
# iCAF CXCL12/14  
# epi

#####FigS10e-h
library(patchwork)
load("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/subcluster/Epithelial_subcluster.Rdata") 
merged.object@meta.data$cell_state=NA
merged.object$cell_state[merged.object$res_0.3%in%c(10)]="Myoepithelial"
merged.object$cell_state[merged.object$res_0.3%in%c(11,12)]="Luminal progenitor"
merged.object$cell_state[merged.object$res_0.3%in%c(0:9)]="Tumor epithelial"
merged.object$cell_state=factor(merged.object$cell_state,levels = c("Tumor epithelial",
             "Luminal progenitor","Myoepithelial"))
p1=DimPlot(merged.object,reduction = "umap_subcluster",group.by = "cell_state",cols="alphabet",
        alph=0.7,label = T)+NoLegend()
features = c("ACTA2","ERBB2","ESR1","KIT","KRT19","TP63")
Idents(Epi)=factor(Idents(Epi),levels = c("Tumor epithelial",
         "Luminal progenitor","Myoepithelial"))
p2=Clustered_DotPlot(seurat_object =Epi, features = features,
                     show_ident_legend = FALSE)
p2=p2[[2]]
ggsave(p1,file="E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS10/subcluster_umap/Epithelial_subcluster_umap.pdf",width=4,height=5)



# Cycling myeloid MKI67 3,12
# Monocyte CD14 CD16 (FCGR3A)   6
# M1 macrophages CD68 CD86, CD80   4,5,9,13,
# M2 macrophages CD163,MRC1   2,11,15,18
# LAM TREM2,FABP4  10,17
# cDC CD1C, ITGAX, CLEC9A  1,16
# pDC CLEC4C  14
# MDSC CD16b CD33 0,7
# Neutrophil CXCR1 CXCR2 FCGR3B FPR2 8
load("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/subcluster/Myeloid_subcluster.Rdata") 
features=c("MKI67","FCGR3A","CD14","CD68","CD86","CD80","CD163","MRC1",
           "TREM2","FABP4","CD1C","ITGAX","CLEC9A","CLEC4C","CD33",
           "CXCR1","CXCR2","FCGR3B","FPR2")
p2=Clustered_DotPlot(seurat_object = Myeloid, features = features,
                     show_ident_legend = FALSE)
p2=p2[[2]]
p2

# mesenchymal
features=c("PDGFRA","PDGFRB","MCAM","ACTA2","FAP","CXCL12","CDH2","SNAI2",
           "POSTN","FN1")
p2=Clustered_DotPlot(seurat_object = Mesenchymal, features = features,
                     show_ident_legend = FALSE)
p2=p2[[2]]
p2

# Endo
# Endothelial res0.3
# Cycling endo 6,12
# Vascular  PECAM1 VWF 0:5,7:8,10:11
# Lymphatic PROX1 9
features=c("MKI67","VWF","PECAM1","PROX1")
p2=Clustered_DotPlot(seurat_object = Endo, features = features,
                     show_ident_legend = FALSE)
p2=p2[[2]]
p2

features=c("CD79A","MS4A1","CD27")
p2=Clustered_DotPlot(seurat_object = Bcell, features = features,
                     show_ident_legend = FALSE)
p2=p2[[2]]
p2

################
# FigureS10k #
################
library(data.table);library(ggplot2);library(scales);library(RColorBrewer);library(tidyverse)
data=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/Xenium_baselineMeta_cell_state_curated.rds")
table(data$cell_state,data$sampleID)
unique(data$cell_type_major_curated)
table(data$cell_type_major_curated,data$sampleID)
table(data$cell_state)
immune=c("T cells","Myeloid cells","Plasma","B cells","Mast cells","Endothelial cells","Mesenchymal cells")
Immune=data[data$cell_type_major_curated%in%immune,]
table(Immune$cell_type_major_curated,Immune$sampleID)
table(Immune$sampleID)
#Immune <- Immune[!Immune$sampleID %in% c("1228_BO"), ] #"1404_BL","1512_BL"
#Immune <- Immune[!Immune$sampleID %in% c("1404_BL"), ] # can include
Immune_proportion <- Immune %>%
  group_by(patientID, cell_state) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(patientID) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup()
Immune_wide <- Immune_proportion %>%
  select(patientID, cell_state, proportion) %>%
  pivot_wider(names_from = cell_state, values_from = proportion, values_fill = 0)
# read meta
meta=readRDS("E:/Projects/PREDIX_HER2/Multimodal/Data/Clin/PREDIX_HER2_clin_curated.rds")
meta$patientID=as.character(meta$patientID)
Immune_wide =left_join(Immune_wide ,meta,by="patientID")
#colnames(tumor_cell_proportion_wide)=gsub(" cancer (Cell|cell|cells)", "", colnames(tumor_cell_proportion_wide))
d=Immune_wide%>%select(c("Arm","Response", setdiff(unique(data$cell_state[data$cell_type_major_curated%in%immune]),"Mast cells"))) # 
d <- reshape2::melt(d,id.vars=c("Arm","Response"))
d$variable=factor(d$variable,levels = c("Plasma","Naive B","Memory B","Cycling T","Exhausted T","Th1","Th2","Th17","CD8 Tem","CTLs",
      "Treg","NK","Cycling myeloid","M1 macrophages","M2 macrophages","Monocyte","LAM","cDC","pDC","Neutrophils","MDSC",
      "myCAF","iCAF","EMT-like CAF","PVLs","Cycling Endothelial","Lymphatic Endothelial","Vascular Endothelial"))
# change the contents of variable baseDir to the root analysis folder 
baseDir <- "E:/Projects/PREDIX_HER2/Multimodal/"
# load directory structure downloaded from github site
source (paste0(baseDir,"/Code/theme.R"))
figure_font_size=12
str(d)
Fig3c <-
  ggplot(d,aes(x=Arm,y=value,fill=Response))+
  geom_boxplot(outlier.size = 0.2, width=0.8)+
  facet_wrap(variable~.,scales="free_y",nrow=2)+
  theme_manuscript(base_size = figure_font_size)+
  scale_fill_pCR_RD(name="Treatment Response")+
  stat_compare_means(aes(group=Response),label = "p.format", hide.ns = F,size=4,
                     color="#000000", label.y.npc = 0.95)+
  labs(y="Immune/stomal cell proportion",x="")+
  theme(strip.background = element_blank(),
        panel.grid.major.x  = element_blank(),
        axis.line = element_blank(),
        strip.text = element_text(size = figure_font_size),
        panel.border = element_rect(colour = "black", fill=NA, size=0.8),
        plot.margin = unit(c(0.5,0.2,0.2,0.5), "lines"))
Fig3c
ggsave(Fig3c, file="E:/Projects/PREDIX_HER2/Multimodal/Figures/FigureS10/FigS10k.pdf", width=28, height=8)


table(d$variable)%>%unique()
