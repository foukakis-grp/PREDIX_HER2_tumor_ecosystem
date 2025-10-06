library(ggplot2, lib.loc = '/sw/apps/R_packages/4.3.1/rackham')
library(Seurat);library(BPCells);library(future);library(dplyr);library(patchwork);library(tibble);library(tidyr)
options(future.globals.maxSize = 1000 * 1024^3)  
data=readRDS("/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/data_harmony_labeled.rds")
table(Idents(data))
data=subset(data,idents=c("B_cells"))
print("Subcluster Bcell")
gc()

## re-integration ##
DefaultAssay(data) <- 'RNA'
data[["RNA"]] <- split(data[["RNA"]], f = data$sampleID)
data <- SCTransform(data, variable.features.n = 5099)
data <- RunPCA(data, npcs = 30, verbose = F)
merged.object <- IntegrateLayers(data, method = HarmonyIntegration, assay = 'SCT',
        orig = "pca", new.reduction = "subcluster.harmony", dims = 1:20, k.anchor = 20, verbose = F)
rm(data)
DefaultAssay(merged.object) <- "RNA"
merged.object[["SCT"]] <- NULL
gc()
merged.object <- FindNeighbors(merged.object, dims = 1:20, reduction = "subcluster.harmony",
                               graph.name = "harm.graph")
merged.object <- RunUMAP(merged.object,
                         reduction = "subcluster.harmony",
                         dims = 1:20,       
                         reduction.name = "umap.subcluster.harmony",
                         n.neighbors = 30,  
                         min.dist = 0.4)  
print("re-integration done!")
## Leiden cluster ##
seq_res <- seq(0.2, 1.5, 0.1)

merged.object <- FindClusters(merged.object, 
                              algorithm = 4,
                              resolution = seq_res,
                              graph.name = "harm.graph")
clustree_plt <- clustree::clustree(merged.object,prefix = paste0("harm.graph_res."))

## Optimal cluster ##
cell_dists <- dist(merged.object@reductions$subcluster.harmony@cell.embeddings,method = "euclidean")
head(cell_dists)

cluster_info <- merged.object@meta.data[,grepl(paste0("harm.graph_res."),colnames(merged.object@meta.data))] %>%
                dplyr::mutate_all(as.character)%>%dplyr::mutate_all(as.numeric)
    
silhouette_res <- apply(cluster_info, 2, function(x){
  si <- cluster::silhouette(x, cell_dists)
  if(!any(is.na(si))) {
    mean(si[, 'sil_width'])
  } else {
    NA
  }
})     
head(silhouette_res)
merged.object[["opt_subclust_integrated"]] <- merged.object[[names(which.max(silhouette_res))]]

save(merged.object, silhouette_res, file = "/proj/sens2022005/Xenium/PREDIX_HER2/result/cell_typing/Bcell_subcluster.Rdata")  
