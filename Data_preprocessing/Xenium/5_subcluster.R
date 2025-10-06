library(Seurat);library(SeuratObject);library(SeuratDisk);library(BPCells)
library(dplyr);library(ggplot2);library(patchwork)
options(future.globals.maxSize = 200 * 1024^3)  
options(max.print = 10000)
 
#Epi res0.4
#ACTA2 TP63 Myoepithelial 10
#KIT KRT19 Luminal progenitor 11,12
#####ESR1 PGR FOXA1 Mature Luminal 2 XXXX
load("/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Epithelial_subcluster.Rdata")
DimPlot(merged.object,reduction = "umap.subcluster.harmony",group.by = "harm.graph_res.0.4",
        alph=0.7,label = T)+NoLegend()
features=c("ACTA2","MKI67","EPCAM","TP63","KIT","KRT19","ERBB2","ESR1","PGR","FOXA1")
FeaturePlot(merged.object,reduction="umap.subcluster.harmony", features = features)
DotPlot(merged.object, features = features,group.by ="harm.graph_res.0.4") + RotatedAxis()
#merged.object <- JoinLayers(merged.object)
#markers <- FindAllMarkers(
#  merged.object,
#  group.by = "harm.graph_res.0.4",
#  test.use = 'wilcox',
#  logfc.threshold = 0.25,
#  only.pos = TRUE
#)
merged.object@meta.data$cell_type=NA
merged.object$cell_type="Tumor epithelial"
merged.object$cell_type[merged.object$harm.graph_res.0.4%in%c(11)]="Normal epithelial"
pal=c("#f3c300","#875692","#f38400","#a1caf1","#be0032","#c2b280","#008856","grey","#e68fac","#0067a5","#f99379","#604e97","#f6a600")
Epi=DimPlot(merged.object,reduction = "umap.subcluster.harmony",group.by = "cell_state",cols=pal,alph=0.7,label = T)+ ggplot2::theme(legend.position = "none")
Epi
ggsave(Epi,file="/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Epi_subcluster.pdf",width=5,height=5)

Epi_meta=merged.object@meta.data[,c("orig.ident","cell_type")]
Epi_meta$cell.name=row.names(Epi_meta)
table(Epi_meta$cell_type)
# add scSubtype
scSubtype=readRDS("/proj/sens2022005/Xenium/PREDIX_HER2/result/scSubtype/scSubtype_result_sep25.rds")
scSubtype$cell.name=row.names(scSubtype)
scSubtype=scSubtype[,c("cell.name","SCSubtypeCall")]
Epi_meta=left_join(Epi_meta,scSubtype,by="cell.name")
table(Epi_meta$SCSubtypeCall,Epi_meta$cell_type)
Epi_meta$cell_state=Epi_meta$SCSubtypeCall
Epi_meta$SCSubtypeCall=NULL
Epi_meta$cell_state[Epi_meta$cell_type=="Normal epithelial"]="Normal epithelial"
table(Epi_meta$cell_state)
colnames(Epi_meta)
Epi_meta=Epi_meta[,c("cell.name","cell_type","cell_state")]
saveRDS(Epi_meta,file="/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Epi_meta.rds")

# B cell res0.3
# Naive B cell  TCL1A 0,4,5
# Memory B cell  CD27,AICDA  1,2,3,6,7,8
load("/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Bcell_subcluster.Rdata")
gc()
DimPlot(merged.object,reduction = "umap.subcluster.harmony",group.by = "harm.graph_res.0.2",
        alph=0.7,label = T)+NoLegend()
features=c("TCL1A", "CD27")
FeaturePlot(merged.object,reduction="umap.subcluster.harmony", features = features,
            min.cutoff = 0, max.cutoff = 10,)
DotPlot(merged.object, features = features,group.by ="harm.graph_res.0.2") + RotatedAxis()

merged.object <- JoinLayers(merged.object)
markers <- FindAllMarkers(
  merged.object,
  group.by = "harm.graph_res.0.2",
  test.use = 'wilcox',
  logfc.threshold = 0.25,
  only.pos = TRUE
)
merged.object@meta.data$cell_type="Bcell"
merged.object@meta.data$cell_state=NA
merged.object$cell_state="Memory B"
merged.object$cell_state[merged.object$harm.graph_res.0.2%in%c(2)]="Naive B"
pal=c("#f3c300","#875692","#f38400","#a1caf1","#be0032","#c2b280","#008856","#e68fac","#0067a5","#f99379","#604e97","#f6a600")
Bcell=DimPlot(merged.object,reduction = "umap.subcluster.harmony",group.by = "cell_state",cols=pal,alph=0.7,label = T)+ ggplot2::theme(legend.position = "none")
Bcell
ggsave(Bcell,file="/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Bcell_subcluster.pdf",width=5,height=5)

Bcell_meta=merged.object@meta.data[,c("cell_type","cell_state")]
Bcell_meta$cell.name=row.names(Bcell_meta)
Bcell_meta=Bcell_meta[,c("cell.name","cell_type","cell_state")]
saveRDS(Bcell_meta,file="/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Bcell_meta.rds")


# Endothelial res0.3
# Cycling endo 6,12
# Vascular  PECAM1 VWF 0:5,7:8,10:11
# Lymphatic PROX1 9
load("/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Endothelial_subcluster.Rdata")
gc()
DimPlot(merged.object,reduction = "umap.subcluster.harmony",group.by = "harm.graph_res.0.2",
        alph=0.7,label = T)+NoLegend()
features=c("PECAM1", "VWF", "PROX1","MKI67")
FeaturePlot(merged.object,reduction="umap.subcluster.harmony", features = features,
            min.cutoff = 0, max.cutoff = 10,)
DotPlot(merged.object, features = features,group.by ="harm.graph_res.0.2") + RotatedAxis()
merged.object@meta.data$cell_type=NA
merged.object$cell_type="Endothelial"
merged.object@meta.data$cell_state=NA
merged.object$cell_state[merged.object$harm.graph_res.0.2%in%c(4)]="Cycling"
merged.object$cell_state[merged.object$harm.graph_res.0.2%in%c(2)]="Inflammatory"  # CXCL9 CXCL10 CXCL11
merged.object$cell_state[merged.object$harm.graph_res.0.2%in%c(6)]="Lymphatic"
merged.object$cell_state[merged.object$harm.graph_res.0.2%in%c(1,3,5,7)]="Vascular"
pal=c("#f3c300","#875692","#f38400","#a1caf1","#be0032","#c2b280","#008856","#e68fac","#0067a5","#f99379","#604e97","#f6a600")
Endothelial=DimPlot(merged.object,reduction = "umap.subcluster.harmony",group.by = "cell_state",cols=pal,alph=0.7,label = T)+ ggplot2::theme(legend.position = "none")
Endothelial
ggsave(Endothelial,file="/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Endothelial_subcluster.pdf",width=5,height=5)

Endothelial_meta=merged.object@meta.data[,c("cell_type","cell_state")]
Endothelial_meta$cell.name=row.names(Endothelial_meta)
Endothelial_meta=Endothelial_meta[,c("cell.name","cell_type","cell_state")]
saveRDS(Endothelial_meta,file="/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Endothelial_meta.rds")

# Mesenchymal res0.3

# PVLs MCAM 4,8,9,10,12
# EMT-like CAF 3
# iCAF CXCL12/14  1,6
# myCAF 0,2,5,7,11
load("/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Mesenchymal_subcluster.Rdata")
features=c("MCAM","PDGFRB","CXCL12","CXCL16","CNN1","FAP","ACTA2")
features=c("FN1","FAP")
FeaturePlot(merged.object,reduction="umap.subcluster.harmony", features = features,
            min.cutoff = 0, max.cutoff = 10,)
DotPlot(merged.object, features = features,group.by ="harm.graph_res.0.2") + RotatedAxis()
DimPlot(merged.object,reduction = "umap.subcluster.harmony",group.by = "harm.graph_res.0.2",
        alph=0.7,label = T)+NoLegend()

merged.object <- JoinLayers(merged.object)
markers <- FindAllMarkers(
  merged.object,
  group.by = "harm.graph_res.0.2",
  test.use = 'wilcox',
  logfc.threshold = 0.25,
  only.pos = TRUE
)
merged.object@meta.data$cell_type=NA
merged.object$cell_type="Mesenchymal"
merged.object@meta.data$cell_state=NA
merged.object$cell_state[merged.object$harm.graph_res.0.2%in%c(5)]="PVLs"
merged.object$cell_state[merged.object$harm.graph_res.0.2%in%c(2)]="iCAF"
merged.object$cell_state[merged.object$harm.graph_res.0.2%in%c(1,4)]="myCAF" #FN1 FAP ACTA2
merged.object$cell_state[merged.object$harm.graph_res.0.2%in%c(3)]="ECM-CAF" # THBS1,CCN1,CCN2,CCN3
pal=c("#f3c300","#875692","#f38400","#a1caf1","#be0032","#c2b280","#008856","#e68fac","#0067a5","#f99379","#604e97","#f6a600")
Mesenchymal=DimPlot(merged.object,reduction = "umap.subcluster.harmony",group.by = "cell_state",cols=pal,alph=0.7,label = T)+ ggplot2::theme(legend.position = "none")
Mesenchymal
ggsave(Mesenchymal,file="/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Mesenchymal_subcluster.pdf",width=5,height=5)

Mesenchymal_meta=merged.object@meta.data[,c("cell_type","cell_state")]
Mesenchymal_meta$cell.name=row.names(Mesenchymal_meta)
Mesenchymal_meta=Mesenchymal_meta[,c("cell.name","cell_type","cell_state")]
saveRDS(Mesenchymal_meta,file="/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Mesenchymal_meta.rds")

#Myeloid res0.5
# Cycling myeloid MKI67 3,12
# Monocyte CD14 CD16 (FCGR3A)   6
# M1 macrophages CD68 CD86, CD80   4,5,9,13,
# M2 macrophages CD163,MRC1   2,11,15,18
# LAM TREM2,FABP4  10,17
# cDC CD1C, ITGAX, CLEC9A  1,16
# pDC CLEC4C  14
# MDSC CD16b CD33 0,7
# Neutrophil CXCR1 CXCR2 FCGR3B FPR2 8
# Mast cell CPA3 MS4A2

load("/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Myeloid_subcluster.Rdata")
#delete 11 4 9
table(merged.object$harm.graph_res.0.6)
Idents(merged.object)=merged.object$harm.graph_res.0.6
merged.object = subset(merged.object, idents = c(4,9,11), invert = TRUE)
DimPlot(merged.object,reduction = "umap.subcluster.harmony",group.by = "harm.graph_res.0.6",
        alph=0.7,label = T)+NoLegend()
features=c("MKI67","TOPA2","FCGR3A","CD14","CD16","CD68","CD86","CD163","TREM2","FABP4","ITGAX",
           "CLEC9A","CLEC4C","CD33","CXCR1","CXCR2","FCGR3B","FPR2","CSF3R","ADGRG3")
DotPlot(merged.object, features = features,group.by ="harm.graph_res.0.6") + RotatedAxis()

features=c("FCGR3A","CD14","CD68","MRC1","FCGR3B","CSF3R","CLEC9A","CLEC4C")
features=c("CD14","CD68")
FeaturePlot(merged.object,reduction="umap.subcluster.harmony", features = features,
            min.cutoff = 0, max.cutoff = 10)

merged.object <- JoinLayers(merged.object)
markers <- FindAllMarkers(
  merged.object,
  group.by = "harm.graph_res.0.6",
  test.use = 'wilcox',
  logfc.threshold = 0.25,
  only.pos = TRUE
)

merged.object@meta.data$cell_type=NA
merged.object$cell_type="Myeloid"
merged.object@meta.data$cell_state=NA
merged.object$cell_state[merged.object$harm.graph_res.0.6%in%c(1)]="M2_macrophages" #CD163
merged.object$cell_state[merged.object$harm.graph_res.0.6%in%c(12)]="M1_macrophages" # IL1B, TNF, IL6, CCL4, CCL8, CXCL2, OSM
merged.object$cell_state[merged.object$harm.graph_res.0.6%in%c(2)]="APOC1_macrophages" #CD163
merged.object$cell_state[merged.object$harm.graph_res.0.6%in%c(10)]="MARCO_macrophages" #CD163
merged.object$cell_state[merged.object$harm.graph_res.0.6%in%c(6)]="cDC1"  # CLEC9A, XCR1  
merged.object$cell_state[merged.object$harm.graph_res.0.6%in%c(3)]="cDC2"  # CD1C CD1D   
merged.object$cell_state[merged.object$harm.graph_res.0.6%in%c(5)]="pDC"  # CLEC4C⁺ IL3RA⁺   
merged.object$cell_state[merged.object$harm.graph_res.0.6%in%c(13)]="Langerhans"  # CD207 CD1A CD93
merged.object$cell_state[merged.object$harm.graph_res.0.6%in%c(7)]="Mast cells"
merged.object$cell_state[merged.object$harm.graph_res.0.6%in%c(8)]="LAM" # SPP1 CD36 LPL PLIN2
merged.object$cell_state[merged.object$harm.graph_res.0.6%in%c(14)]="ECM_macrophages" # CTSK⁺ MMP9⁺ ACP5⁺ ITGB3⁺ NFATC1

pal=c("#f3c300","#875692","#f38400","#a1caf1","#be0032","#c2b280","#008856","#e68fac","#0067a5","#f99379","#604e97","#f6a600")
Myeloid=DimPlot(merged.object,reduction = "umap.subcluster.harmony",group.by = "cell_state",cols=pal,alph=0.7,label = T)+ ggplot2::theme(legend.position = "none")
Myeloid
ggsave(Myeloid,file="/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Myeloid_subcluster.pdf",width=5,height=5)

DimPlot(merged.object,reduction = "umap.subcluster.harmony",group.by = "cell_state",
        alph=0.7,label = T)+NoLegend()


Myeloid_meta=merged.object@meta.data[,c("cell_type","cell_state")]
Myeloid_meta$cell.name=row.names(Myeloid_meta)
Myeloid_meta=Myeloid_meta[,c("cell.name","cell_type","cell_state")]
saveRDS(Myeloid_meta,file="/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Myeloid_meta.rds")


####################################
####################################
####################################
# T cell res0.5
# Exhausted T cell PDCD1 12 
# Cycling T cell MKI67 9 

# CD8T CD8A 
#### Naive CD8 T cell, CCR7 SELL XXX
# CTLs GZMK 2,11,17
# CD8 Tem CCR7,IL7R 5  

# CD4T 0,3,4,12,15
# CD4 naive SELL CCR7 IL7R  4,13
#Th1 IFNG CXCR3 CXCL9 CCL20 0
#Th2 GATA3,IL4R,STAT6,CCL17 13
#Th17 RORC 4
#Treg FOXP3 3
# NK cell KLRD1,NCAM1 10

load("/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Tcell_subcluster_final.Rdata")
gc()
DimPlot(merged.object,reduction = "umap.subcluster.harmony",group.by = "harm.graph_res.0.6",
        alph=0.7,label = T)+NoLegend()
DotPlot(merged.object, features = features,group.by ="harm.graph_res.0.8") + RotatedAxis()
merged.object <- JoinLayers(merged.object)
markers <- FindAllMarkers(
  merged.object,
  group.by = "harm.graph_res.0.8",
  test.use = 'wilcox',
  logfc.threshold = 0.25,
  only.pos = TRUE
)

# delete 4 6 8 9 12 13 14
table(merged.object$harm.graph_res.0.8)
merged.object@meta.data$cell_type=NA
merged.object@meta.data$cell_type="Tcell"
merged.object@meta.data$cell_state=NA         # CD3E, CD4, CD8A
merged.object$cell_state[merged.object$harm.graph_res.0.8%in%c(1)]="Th17"   # RORC CD40LG CCR6、KLRB1（CD161）、DPP4
merged.object$cell_state[merged.object$harm.graph_res.0.8%in%c(2,5)]="CTLs"   # GZMA GZMK KLRK1 CCL4 
merged.object$cell_state[merged.object$harm.graph_res.0.8%in%c(3)]="CD8_Trm"   # CD8A MIAT ITGA1 ITGAE 
merged.object$cell_state[merged.object$harm.graph_res.0.8%in%c(4)]="NK"  # FCGR3A "KLRK1","NCAM1"
merged.object$cell_state[merged.object$harm.graph_res.0.8%in%c(6,11)]="Treg"  # FOXP3 CD4 
merged.object$cell_state[merged.object$harm.graph_res.0.8%in%c(7)]="CXCR4_T"  # CD4 Naive
merged.object$cell_state[merged.object$harm.graph_res.0.8%in%c(8,12,13,14)]="Tfh"  # CXCL13 CXCR5 BCL6 ICOS
merged.object$cell_state[merged.object$harm.graph_res.0.8%in%c(9)]="Cycling T"  # TUBB TK1 MKI67  
merged.object$cell_state[merged.object$harm.graph_res.0.8%in%c(10)]="NKT"  # 
pal=c("#f3c300","#875692","#f38400","#a1caf1","#008856","grey","#e68fac","#0067a5","#f99379","#604e97","#f6a600","#be0032")
Tcell=DimPlot(merged.object,reduction = "umap.subcluster.harmony",group.by = "cell_state",cols=pal,alph=0.7,label = T)+ ggplot2::theme(legend.position = "none")
Tcell
ggsave(Tcell,file="/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Tcell_subcluster.pdf",width=5,height=5)



features=c("TRAC","TUBB","MKI67","TK1","CD4","CD8A","CD8B","GZMB","GZMK","PRF1","CCL5","IFNG",'MIAT',"ITGA1","ITGAE",  
           "PDCD1","LAG3","TIGIT","IL7R","CCR7","SELL","TCF7","LEF1","CD27","TBX21","STAT4","CXCR3",
           "GATA3","IL4","IL5","IL13","CCR4","RORC","IL17A","IL17F","CCR6","IL21","IL22","FOXP3",
           "IL2RA","CTLA4","IKZF2","CCR8","NCAM1","NCR1","KLRK1","FCGR3A","BCL6","CXCL13")
DotPlot(merged.object, features = features,group.by ="cell_state") + RotatedAxis()


Tcell_meta=merged.object@meta.data[,c("cell_type","cell_state")]
Tcell_meta$cell.name=row.names(Tcell_meta)
Tcell_meta=Tcell_meta[,c("cell.name","cell_type","cell_state")]
saveRDS(Tcell_meta,file="/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Tcell_meta.rds")



###########merge cell state with Seurat#############
data=readRDS("/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/data_harmony_labeled.rds")
nrow(data@meta.data)
str(data@meta.data)
Epi_meta=readRDS("/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Epi_meta.rds") 
Tcell_meta=readRDS("/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Tcell_meta.rds") 
Bcell_meta=readRDS("/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Bcell_meta.rds") 
Endothelial_meta=readRDS("/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Endothelial_meta.rds") 
Myeloid_meta=readRDS("/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Myeloid_meta.rds") 
Mesenchymal_meta=readRDS("/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Mesenchymal_meta.rds") 
meta=rbind(Epi_meta,Tcell_meta)%>%rbind(Bcell_meta)%>%rbind(Endothelial_meta)%>%rbind(Myeloid_meta)%>%rbind(Mesenchymal_meta)
data@meta.data$cell.name=rownames(data@meta.data)
a=data@meta.data%>%as.data.frame()
meta=left_join(a,meta,by="cell.name")
meta$cell_type[is.na(meta$cell_type)]="None"
data@meta.data$cell_type=NA
data$cell_type=meta$cell_type
data@meta.data$cell_state=NA
data$cell_state=meta$cell_state
Idents(data)=data$cell_type
table(Idents(data))
merged.object=subset(data,idents=c("Bcell","Endothelial","Mesenchymal","Myeloid","Normal epithelial","Tcell","Tumor epithelial"))

# change the dir of BPcells

merged.object@assays[["RNA"]]@layers[["counts"]]@matrix@matrix_list[[1]]@matrix@matrix@matrix@dir

# change BP dir
for (i in seq_along(merged.object@assays[["RNA"]]@layers[["counts"]]@matrix@matrix_list)) {
  old_dir <- merged.object@assays[["RNA"]]@layers[["counts"]]@matrix@matrix_list[[i]]@matrix@matrix@matrix@dir
  
  # replace dir
  new_dir <- sub(
    pattern = "^/proj/sens2022005/Xenium/PREDIX_HER2/data/BP_sep25/",
    replacement = "E:/Projects/PREDIX_HER2/Multimodal/Data/Xenium/BPcell/",
    x = old_dir
  )
  
  # updated dir
  merged.object@assays[["RNA"]]@layers[["counts"]]@matrix@matrix_list[[i]]@matrix@matrix@matrix@dir <- new_dir
}


merged.object@assays[["RNA"]]@layers[["counts"]]@matrix@matrix_list[[2]]@matrix@matrix@matrix@dir
saveRDS(merged.object,"/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Seurat_cell_state_curated_localBPcell_Oct2025.rds")



