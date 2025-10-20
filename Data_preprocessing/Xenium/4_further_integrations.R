library(ggplot2, lib.loc = '/sw/apps/R_packages/4.3.1/rackham')
library(Seurat);library(BPCells);library(future);library(dplyr);library(patchwork);library(tibble);library(tidyr)
options(future.globals.maxSize = 1000 * 1024^3)  

### First harmony run, second integration in general
merged.object <- readRDS('/proj/sens2022005/Xenium/PREDIX_HER2/result/integration/sketch_sep25/full_data_integrated_clustered.rds')
delete <- c('23', '51', '55', '64', '44', '48', '56', '58', '30', '15', '52', '34', '39', '67', '57', '70', '16', '65', '11', '35')
Idents(merged.object) <- merged.object$harmony_2.5
merged.object <- subset(merged.object, idents = delete, invert=TRUE)

DefaultAssay(merged.object) <- 'RNA'
merged.object <- SCTransform(merged.object, variable.features.n = 5099)
merged.object <- RunPCA(merged.object, npcs = 30, verbose = F)
merged.object <- IntegrateLayers(merged.object, method = HarmonyIntegration, assay = 'SCT',
                                 orig = "pca", new.reduction = "integrated.harmony", dims = 1:20, k.anchor = 20, verbose = F)
DefaultAssay(merged.object) <- "RNA"
merged.object[["SCT"]] <- NULL

merged.object <- NormalizeData(merged.object)
merged.object <- FindVariableFeatures(merged.object, nfeatures = 5099)
merged.object <- ScaleData(merged.object)

merged.object <- FindNeighbors(merged.object, dims = 1:20, reduction = "integrated.harmony",
                               graph.name = "harm.graph")
merged.object <- RunUMAP(merged.object,
                         reduction = "integrated.harmony",
                         dims = 1:20,       
                         reduction.name = "umap.harmony",
                         n.neighbors = 30,  
                         min.dist = 0.4)  
merged.object <- FindClusters(merged.object, 
                              algorithm = 4,
                              resolution = 4.0,
                              cluster.name = "harmony_4.0",
                              graph.name = "harm.graph")

merged.object <- JoinLayers(merged.object)
markers_list <- list(
  Adipocytes = c('FABP4', 'PLIN1', 'PLIN4', 'ADIPOQ', 'CD36', 'LPL', 'CIDEC', 'PPARG','PLIN2'),
  B_cells = c('CD79A', 'MS4A1', 'POU2AF1', 'CD19', 'BANK1', 'TNFRSF13C', 'CD79B', 'CD22'),
  Endothelial = c('PECAM1', 'VWF', 'PLVAP', 'ENG', 'CD34', 'KDR', 'FLT1', 'TIE1', 'EGFL7', 'CLEC14A', 'ROBO4','PROX1', 'ICAM2'),
  Epithelial = c('EPCAM', 'KRT19', 'CDH1', 'MUC1', 'FOXA1', 'ESR1', 'GATA3', 'ERBB2', 'ERBB3', 'GRB7'),
  Mast = c('CPA3', 'MS4A2', 'KIT', 'HDC', 'SIGLEC6', 'IL1RL1', 'CMA1','SLC18A2','GATA2'),
  T_cells = c('CD3E', 'CD3G', 'TRAC', 'IL7R', 'CD2', 'CD8A', 'CD4', 'CCL5', 'CD247', 'CCR7','GZMB','GZMK'),
  Plasma = c('MZB1', 'SLAMF7', 'CD38', 'XBP1', 'TENT5C', 'IRF4','PRDM1','DERL3'),
  Myeloid = c('CD14', 'FCGR3A', 'FCGR2A', 'MSR1', 'APOC1', 'CD68', 'ITGAX', 'FCN1', 'SLCO2B1'),
  Mesenchymal = c('FAP', 'PDPN', 'PDGFRB', 'ACTA2', 'MYH11', 'RGS5', 'MCAM', 'COL5A1','COL5A2','THY1','NOTCH3', 'COL11A1', 'THBS2')
)
for (celltype in names(markers_list)) {
  merged.object <- AddModuleScore(
    object = merged.object,
    features = list(markers_list[[celltype]]),
    name = celltype)}
score_cols <- paste0(names(markers_list), "1")
avg_scores <- FetchData(merged.object, vars = c("harmony_4.0", score_cols)) %>%
  group_by(harmony_4.0) %>%
  summarise(across(all_of(score_cols), mean))
cluster_assignments <- avg_scores %>%
  tidyr::pivot_longer(cols = -harmony_4.0, names_to = "celltype", values_to = "score") %>%
  group_by(harmony_4.0) %>%
  slice_max(score, n = 1) %>%
  ungroup()
cluster_ids <- setNames(gsub("1", "", cluster_assignments$celltype),
                        cluster_assignments$harmony_4.0)
Idents(merged.object) <- merged.object$harmony_4.0
merged.object <- RenameIdents(merged.object, cluster_ids)

DimPlot(merged.object, reduction = "umap.harmony", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(merged.object,reduction = "umap.harmony",group.by="harmony_4.0",alph=0.5, label = T)



### Second harmony run, third integration in total
delete <- c('48', '25', '53', '43', '72', '78')
Idents(merged.object) <- merged.object$harmony_4.0
merged.object <- subset(merged.object, idents = delete, invert = TRUE)

merged.object <- SCTransform(merged.object, variable.features.n = 5099)
print('SCT done')
merged.object <- RunPCA(merged.object, npcs = 30, verbose = F)
merged.object <- IntegrateLayers(merged.object, method = HarmonyIntegration, assay = 'SCT',
                                 orig = "pca", new.reduction = "integrated.harmony2", dims = 1:20, k.anchor = 20, verbose = F)
DefaultAssay(merged.object) <- "RNA"
merged.object[["SCT"]] <- NULL

merged.object <- NormalizeData(merged.object)
merged.object <- FindVariableFeatures(merged.object, nfeatures = 5099)
merged.object <- ScaleData(merged.object)

merged.object <- FindNeighbors(merged.object, dims = 1:20, reduction = "integrated.harmony2",
                               graph.name = "harm.graph2")
merged.object <- RunUMAP(merged.object,
                         reduction = "integrated.harmony2",
                         dims = 1:20,       
                         reduction.name = "umap.harmony2",
                         n.neighbors = 30,  
                         min.dist = 0.4)  
merged.object <- FindClusters(merged.object, 
                              algorithm = 4,
                              resolution = 4.0,
                              cluster.name = "harmony2_4.0",
                              graph.name = "harm.graph2")

Idents(merged.object) <- merged.object$harmony2_4.0
merged.object <- JoinLayers(merged.object)
markers_list <- list(
  Adipocytes = c('FABP4', 'PLIN1', 'PLIN4', 'ADIPOQ', 'CD36', 'LPL', 'CIDEC', 'PPARG','PLIN2'),
  B_cells = c('CD79A', 'MS4A1', 'POU2AF1', 'CD19', 'BANK1', 'TNFRSF13C', 'CD79B', 'CD22'),
  Endothelial = c('PECAM1', 'VWF', 'PLVAP', 'ENG', 'CD34', 'KDR', 'FLT1', 'TIE1', 'EGFL7', 'CLEC14A', 'ROBO4','PROX1', 'ICAM2'),
  Epithelial = c('EPCAM', 'KRT19', 'CDH1', 'MUC1', 'FOXA1', 'ESR1', 'GATA3', 'ERBB2', 'ERBB3', 'GRB7'),
  Mast = c('CPA3', 'MS4A2', 'KIT', 'HDC', 'SIGLEC6', 'IL1RL1', 'CMA1','SLC18A2','GATA2'),
  T_cells = c('CD3E', 'CD3G', 'TRAC', 'IL7R', 'CD2', 'CD8A', 'CD4', 'CCL5', 'CD247', 'CCR7','GZMB','GZMK'),
  Plasma = c('MZB1', 'SLAMF7', 'CD38', 'XBP1', 'TENT5C', 'IRF4','PRDM1','DERL3'),
  Myeloid = c('CD14', 'FCGR3A', 'FCGR2A', 'MSR1', 'APOC1', 'CD68', 'ITGAX', 'FCN1', 'SLCO2B1'),
  Mesenchymal = c('FAP', 'PDPN', 'PDGFRB', 'ACTA2', 'MYH11', 'RGS5', 'MCAM', 'COL5A1','COL5A2','THY1','NOTCH3', 'COL11A1', 'THBS2')
)
for (celltype in names(markers_list)) {
  merged.object <- AddModuleScore(
    object = merged.object,
    features = list(markers_list[[celltype]]),
    name = celltype)}
score_cols <- paste0(names(markers_list), "1")
avg_scores <- FetchData(merged.object, vars = c("harmony2_4.0", score_cols)) %>%
  group_by(harmony2_4.0) %>%
  summarise(across(all_of(score_cols), mean))
cluster_assignments <- avg_scores %>%
  tidyr::pivot_longer(cols = -harmony2_4.0, names_to = "celltype", values_to = "score") %>%
  group_by(harmony2_4.0) %>%
  slice_max(score, n = 1) %>%
  ungroup()
cluster_ids <- setNames(gsub("1", "", cluster_assignments$celltype),
                        cluster_assignments$harmony2_4.0)
Idents(merged.object) <- merged.object$harmony2_4.0
merged.object <- RenameIdents(merged.object, cluster_ids)

DimPlot(merged.object, reduction = "umap.harmony2", label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(merged.object, reduction = "umap.harmony2", group.by = 'harmony2_4.0', label = TRUE, pt.size = 0.5) + NoLegend()



### Third harmony run, fourth integration in total
delete <- c('83', '56')
Idents(merged.object) <- merged.object$harmony2_4.0
merged.object <- subset(merged.object, idents = delete, invert = TRUE)

merged.object <- SCTransform(merged.object, variable.features.n = 5099)
merged.object <- RunPCA(merged.object, npcs = 30, verbose = F)
merged.object <- IntegrateLayers(merged.object, method = HarmonyIntegration, assay = 'SCT',
                                 orig = "pca", new.reduction = "integrated.harmony3", dims = 1:20, k.anchor = 20, verbose = F)
DefaultAssay(merged.object) <- "RNA"
merged.object[["SCT"]] <- NULL

merged.object <- NormalizeData(merged.object)
merged.object <- FindVariableFeatures(merged.object, nfeatures = 5099)
merged.object <- ScaleData(merged.object)

merged.object <- FindNeighbors(merged.object, dims = 1:20, reduction = "integrated.harmony3",
                               graph.name = "harm.graph3")
merged.object <- RunUMAP(merged.object,
                         reduction = "integrated.harmony3",
                         dims = 1:20,       
                         reduction.name = "umap.harmony3",
                         n.neighbors = 30,  
                         min.dist = 0.5)  
merged.object <- FindClusters(merged.object, 
                              algorithm = 4,
                              resolution = 4.0,
                              cluster.name = "harmony3_4.0",
                              graph.name = "harm.graph3")

Idents(merged.object) <- merged.object$harmony3_4.0
merged.object <- JoinLayers(merged.object)
markers_list <- list(
  Adipocytes = c('FABP4', 'PLIN1', 'PLIN4', 'ADIPOQ', 'CD36', 'LPL', 'CIDEC', 'PPARG','PLIN2'),
  B_cells = c('CD79A', 'MS4A1', 'POU2AF1', 'CD19', 'BANK1', 'TNFRSF13C', 'CD79B', 'CD22'),
  Endothelial = c('PECAM1', 'VWF', 'PLVAP', 'ENG', 'CD34', 'KDR', 'FLT1', 'TIE1', 'EGFL7', 'CLEC14A', 'ROBO4','PROX1', 'ICAM2'),
  Epithelial = c('EPCAM', 'KRT19', 'CDH1', 'MUC1', 'FOXA1', 'ESR1', 'GATA3', 'ERBB2', 'ERBB3', 'GRB7'),
  Mast = c('CPA3', 'MS4A2', 'KIT', 'HDC', 'SIGLEC6', 'IL1RL1', 'CMA1','SLC18A2','GATA2'),
  T_cells = c('CD3E', 'CD3G', 'TRAC', 'IL7R', 'CD2', 'CD8A', 'CD4', 'CCL5', 'CD247', 'CCR7','GZMB','GZMK'),
  Plasma = c('MZB1', 'SLAMF7', 'CD38', 'XBP1', 'TENT5C', 'IRF4','PRDM1','DERL3'),
  Myeloid = c('CD14', 'FCGR3A', 'FCGR2A', 'MSR1', 'APOC1', 'CD68', 'ITGAX', 'FCN1', 'SLCO2B1'),
  Mesenchymal = c('FAP', 'PDPN', 'PDGFRB', 'ACTA2', 'MYH11', 'RGS5', 'MCAM', 'COL5A1','COL5A2','THY1','NOTCH3', 'COL11A1', 'THBS2')
)
for (celltype in names(markers_list)) {
  merged.object <- AddModuleScore(
    object = merged.object,
    features = list(markers_list[[celltype]]),
    name = celltype)}
score_cols <- paste0(names(markers_list), "1")
avg_scores <- FetchData(merged.object, vars = c("harmony3_4.0", score_cols)) %>%
  group_by(harmony3_4.0) %>%
  summarise(across(all_of(score_cols), mean))
cluster_assignments <- avg_scores %>%
  tidyr::pivot_longer(cols = -harmony3_4.0, names_to = "celltype", values_to = "score") %>%
  group_by(harmony3_4.0) %>%
  slice_max(score, n = 1) %>%
  ungroup()
cluster_ids <- setNames(gsub("1", "", cluster_assignments$celltype),
                        cluster_assignments$harmony3_4.0)
Idents(merged.object) <- merged.object$harmony3_4.0
merged.object <- RenameIdents(merged.object, cluster_ids)

DimPlot(merged.object, reduction = "umap.harmony3", group.by = 'harmony3_4.0', label = TRUE, pt.size = 0.5) + NoLegend()
DimPlot(merged.object, reduction = "umap.harmony3", label = TRUE, pt.size = 0.5) + NoLegend()

saveRDS(merged.object, '/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/data_harmony_labeled.rds')
