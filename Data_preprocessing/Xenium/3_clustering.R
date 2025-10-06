library(ggplot2, lib.loc = '/sw/apps/R_packages/4.3.1/rackham')
library(Seurat);library(BPCells);library(dplyr);library(Polychrome);library(tibble);library(tidyr);library(patchwork)
options(future.globals.maxSize = 300 * 1024^3)

# Label cell types with module score
merged.object <- readRDS('/proj/sens2022005/Xenium/PREDIX_HER2/result/integration/sketch_sep25/full_data_integrated_clustered.rds')
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
avg_scores <- FetchData(merged.object, vars = c("cca_recalc_leiden_4.0", score_cols)) %>%
  group_by(cca_recalc_leiden_4.0) %>%
  summarise(across(all_of(score_cols), mean))
cluster_assignments <- avg_scores %>%
  tidyr::pivot_longer(cols = -cca_recalc_leiden_4.0, names_to = "celltype", values_to = "score") %>%
  group_by(cca_recalc_leiden_4.0) %>%
  slice_max(score, n = 1) %>%
  ungroup()
cluster_ids <- setNames(gsub("1", "", cluster_assignments$celltype),
                        cluster_assignments$cca_recalc_leiden_4.0)
Idents(merged.object) <- merged.object$cca_recalc_leiden_4.0
merged.object <- RenameIdents(merged.object, cluster_ids)

features <- c("Adipocytes1", "B_cells1", "Endothelial1", "Epithelial1", "Mesenchymal1", "Mast1", "Myeloid1", "Plasma1", "T_cells1")  
common_min <- 0
common_max <- 0.5
color_palette <- c('lightgrey', 'blue')
plots <- FeaturePlot(
  object = merged.object,
  features = features,
  reduction = "umap.cca.full",
  pt.size = 0.5,
  combine = FALSE)
plots <- lapply(plots, function(p) {
  p + scale_color_gradientn(
    colors = color_palette,
    limits = c(common_min, common_max),
    oob = scales::squish)})
wrap_plots(plots)
saveRDS(merged.object, '/proj/sens2022005/Xenium/PREDIX_HER2/result/integration/sketch_sep25/full_data_integrated_clustered_labeled.rds')

### Now move on to isolate individual cell types and identify mixed clusters to filter
### Example of workflow with one cell type:

## Mesenchymal
mesench <- subset(merged.object, idents = 'Mesenchymal')
mesench_counts <- as.data.frame(table(mesench$cca_recalc_leiden_4.0))

DimPlot(mesench, reduction = "umap.cca.full", label = TRUE, pt.size = 0.5, group.by = 'cca_recalc_leiden_4.0') + NoLegend()
VlnPlot(mesench, features = 'nCount_Xenium', group.by = 'cca_recalc_leiden_4.0') + NoLegend()


# heatmap to see the values for module scores in each cell type
score_cols <- c("Adipocytes1", "B_cells1", "Endothelial1", "Epithelial1", "Mesenchymal1", "Mast1", "Myeloid1", "Plasma1", "T_cells1")
avg_scores_mesench <- mesench@meta.data %>%
  dplyr::select(cca_recalc_leiden_4.0, all_of(score_cols)) %>%
  group_by(cca_recalc_leiden_4.0) %>%
  summarise(across(all_of(score_cols), mean, .names = "{col}"), .groups = "drop") %>%
  arrange(desc(Mesenchymal1)) %>%
  as.data.frame()
avg_scores_long <- avg_scores_mesench %>%
  pivot_longer(cols = all_of(score_cols),
               names_to = "celltype",
               values_to = "score")
avg_scores_long$cca_recalc_leiden_4.0 <- factor(
  avg_scores_long$cca_recalc_leiden_4.0,
  levels = avg_scores_mesench$cca_recalc_leiden_4.0)
avg_scores_long <- avg_scores_long %>%
  group_by(cca_recalc_leiden_4.0) %>%
  mutate(mesench_score = score[celltype == "Mesenchymal1"]) %>%
  ungroup() %>%
  mutate(cca_recalc_leiden_4.0 = forcats::fct_reorder(cca_recalc_leiden_4.0,
                                                      mesench_score,
                                                      .desc = TRUE))
ggplot(avg_scores_long, aes(x = celltype, 
                            y = cca_recalc_leiden_4.0, 
                            fill = score)) +
  geom_tile() +
  geom_text(aes(label = round(score, 2)), size = 6) +  # show scores, 2 decimals
  scale_fill_gradientn(colours = c("white", "red")) +
  labs(x = "Module score", y = "Cluster", fill = "Avg score") +
  theme_minimal(base_size = 16) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())


# expression of module scores
features <- c("Adipocytes1", "B_cells1", "Endothelial1", "Epithelial1", "Mesenchymal1", "Mast1", "Myeloid1", "Plasma1", "T_cells1")  
common_min <- 0
common_max <- 0.5
color_palette <- c('lightgrey', 'blue')
plots <- FeaturePlot(
  object = mesench,
  features = features,
  pt.size = 0.5,
  combine = FALSE)
plots <- lapply(plots, function(p) {
  p + scale_color_gradientn(
    colors = color_palette,
    limits = c(common_min, common_max),
    oob = scales::squish)})
wrap_plots(plots)


# visualize expression of different markers
features <- c('MS4A1', 'CD79A', 'PECAM1', 'KDR', 'EPCAM', 'sct_ERBB2', 'KRT19', 'FAP', 'COL11A1', 'CD68', 'KIT', 'TENT5C', 'MZB1', 'PDGFRB', 'CD3E', 'TRAC')
common_min <- 0
common_max <- 1
color_palette <- c('lightgrey', 'blue')
plots <- FeaturePlot(
  object = mesench,
  features = features,
  pt.size = 0.5,
  combine = FALSE)
plots <- lapply(plots, function(p) {
  p + scale_color_gradientn(
    colors = color_palette,
    limits = c(common_min, common_max),
    oob = scales::squish)})
wrap_plots(plots)


# inspect location of each suspisious cluster
Idents(mesench) <- mesench$cca_recalc_leiden_4.0
DimPlot(mesench, group.by = "cca_recalc_leiden_4.0", 
        cells.highlight = WhichCells(mesench, idents = "64"), 
        cols.highlight = "red", cols = "lightgrey") +
  ggtitle('cluster 64') + NoLegend()
DimPlot(mesench, group.by = "cca_recalc_leiden_4.0", 
        cells.highlight = WhichCells(mesench, idents = "43"), 
        cols.highlight = "red", cols = "lightgrey") +
  ggtitle('cluster 43') + NoLegend()