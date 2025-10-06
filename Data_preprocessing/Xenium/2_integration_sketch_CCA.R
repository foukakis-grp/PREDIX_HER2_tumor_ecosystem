library(ggplot2, lib.loc = '/sw/apps/R_packages/4.3.1/rackham')
library(Seurat);library(BPCells);library(dplyr)
options(future.globals.maxSize = 300 * 1024^3)  

setwd('/proj/sens2022005/Xenium/PREDIX_HER2/result/integration/sketch_sep25/')

# import metadata to obtain all BL samples
meta <- read.csv('/proj/sens2022005/Xenium/PREDIX_HER2/data/xenium_meta.csv')
selected_sampleIDs <- meta$sampleID[meta$tpt == "pre"]
sampleIDs <- paste0(selected_sampleIDs, "_seurat_QC.rds")

file.dir <- "/proj/sens2022005/Xenium/PREDIX_HER2/result/QC/seurat_object_sep25/"
bp.dir <- "/proj/sens2022005/Xenium/PREDIX_HER2/data/BP_sep25/"

# Lists to store BPCells matrices and metadata
data.list <- list()
metadata.list <- list()

# Loop through Xenium outputs
for (i in seq_along(sampleIDs)) {
  path <- file.path(file.dir, sampleIDs[i])
  
  # Load Xenium data as a Seurat object
  obj <- readRDS(path)
  
  # Convert and store counts as BPCells on-disk matrix
  bp_dir <- file.path(bp.dir, paste0(selected_sampleIDs[i], "_BP"))
  write_matrix_dir(mat = obj[["Xenium"]]$counts, dir = bp_dir)
  
  # Load BPCells matrix
  mat <- open_matrix_dir(bp_dir)
  data.list[[i]] <- mat
 
  # Extract metadata
  metadata <- obj@meta.data
  metadata$orig.ident <- selected_sampleIDs[i]
  metadata$sampleID <- selected_sampleIDs[i]
  metadata.list[[i]] <- metadata
}
names(data.list) <- selected_sampleIDs
metadata <- do.call(rbind, metadata.list)
merged.object <- CreateSeuratObject(counts = data.list, meta.data = metadata)

merged.object[["RNA"]] <- JoinLayers(merged.object[["RNA"]])
merged.object <- NormalizeData(merged.object)
merged.object <- FindVariableFeatures(merged.object, nfeatures = 4099)
merged.object[["RNA"]] <- split(merged.object[["RNA"]], f = merged.object$sampleID) 
merged.object <- NormalizeData(merged.object)
merged.object <- FindVariableFeatures(merged.object, nfeatures = 4099)

merged.object <- SketchData(
  object = merged.object,
  ncells = 10000,
  method = "LeverageScore",
  sketched.assay = "sketch")

DefaultAssay(merged.object) <- "sketch"
merged.object <- FindVariableFeatures(merged.object, nfeatures = 4099)
merged.object <- ScaleData(merged.object)
merged.object <- RunPCA(merged.object)

# test RPCA and CCA
merged.object <- IntegrateLayers(merged.object, method = CCAIntegration, 
                                 orig = "pca", new.reduction = "integrated.cca", dims = 1:20, k.anchor = 20, verbose = F)
merged.object <- IntegrateLayers(merged.object, method = RPCAIntegration, 
                                 orig = "pca", new.reduction = "integrated.rpca", dims = 1:20, k.anchor = 20, verbose = F)

# test different resolutions and algorithms to choose best
reses    <- c(1.5, 2.0, 2.5, 3.0)
algos    <- c(louvain_refined = 2, leiden = 4)

### RPCA branch
# 1) neighbors on integrated.rpca
merged.object <- FindNeighbors(merged.object, 
                               reduction = "integrated.rpca", 
                               dims = 1:20,
                               graph.name = "RPCA_snn", 
                               verbose = FALSE)

# 2) clustering grid (algorithms × resolutions)
for (an in names(algos)) {
  for (r in reses) {
    col <- sprintf("rpca_%s_r%.1f", an, r)
    merged.object <- FindClusters(
      merged.object,
      graph.name   = "RPCA_snn",
      algorithm    = algos[[an]],
      resolution   = r,
      cluster.name = col,
      verbose      = FALSE
    )
  }
}

# 3) UMAP on integrated.rpca
merged.object <- RunUMAP(
  merged.object,
  reduction = "integrated.rpca",
  dims = 1:20,
  n.neighbors = 30,
  min.dist = 0.2,
  return.model = TRUE,
  reduction.name = "umap.rpca",
  reduction.key  = "umaprpca_",
  verbose = FALSE)


### CCA branch

# 1) neighbors on integrated.cca
merged.object <- FindNeighbors(merged.object, 
                               reduction = "integrated.cca", 
                               dims = 1:20,
                               graph.name = "CCA_snn",
                               verbose = FALSE)

# 2) clustering grid
for (an in names(algos)) {
  for (r in reses) {
    col <- sprintf("cca_%s_r%.1f", an, r)
    merged.object <- FindClusters(
      merged.object,
      graph.name   = "CCA_snn",
      algorithm    = algos[[an]],
      resolution   = r,
      cluster.name = col,
      verbose      = FALSE
    )
  }
}

# 3) UMAP on integrated.cca
merged.object <- RunUMAP(
  merged.object,
  reduction = "integrated.cca",
  dims = 1:20,
  n.neighbors = 30,
  min.dist = 0.2,
  return.model = TRUE,
  reduction.name = "umap.cca",
  reduction.key  = "umapcca_",
  verbose = FALSE)


## Projection
merged.object <- ProjectIntegration(object = merged.object,
                                    sketched.assay = "sketch",
                                    assay = "RNA",
                                    reduction = "integrated.cca")

merged.object <- ProjectData(
  object = merged.object,
  assay = "RNA",
  full.reduction = "integrated.cca.full",
  sketched.assay = "sketch",
  sketched.reduction = "integrated.cca",
  dims = 1:30,
  refdata = list(leiden_2.5_full = "cca_leiden_r2.5"), 
  verbose = TRUE
)
DefaultAssay(merged.object)="RNA"
merged.object <- RunUMAP(merged.object, 
                         reduction = "integrated.cca.full", 
                         dims = 1:30, 
                         n.neighbors = 70,  
                         min.dist = 0.1,
                         reduction.name = "umap.cca.full",
                         reduction.key = "umapcca_full_")
merged.object <- FindClusters(merged.object,
                              graph.name = "CCA_full_snn",
                              algorithm = 4,  
                              resolution = 4.0,
                              cluster.name = "cca_recalc_leiden_4.0") 
saveRDS(merged.object, '/proj/sens2022005/Xenium/PREDIX_HER2/result/integration/sketch_sep25/full_data_integrated_clustered.rds')

# Continue with module score for cluster identification in 3_clustering.R
