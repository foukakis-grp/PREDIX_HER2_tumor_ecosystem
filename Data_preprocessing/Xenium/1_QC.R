library(Seurat);library(dplyr);library(arrow)
options(future.globals.maxSize = 500 * 1024^3)
# import metadata to obtain all BL samples
meta <- read.csv('/proj/sens2022005/Xenium/PREDIX_HER2/data/xenium_meta.csv')
selected_samples <- meta$sampleID[meta$tpt == "pre"]
data_path <- "/proj/sens2022005/Xenium/PREDIX_HER2/data/Xenium_ranger_output"
selected_files <- file.path(data_path, selected_samples)
# read Seurat
seurat_list <- lapply(selected_files, LoadXenium)
# name Seurat
names(seurat_list) <- selected_samples
# update `orig.ident`
for (i in seq_along(seurat_list)) {
  seurat_list[[i]]$orig.ident <- selected_samples[i]
  orig_id <- names(seurat_list)[i]
  seurat_list[[i]] <- RenameCells(seurat_list[[i]], new.names = paste0(orig_id, "_", Cells(seurat_list[[i]])))
}

output_dir <- "/proj/sens2022005/Xenium/PREDIX_HER2/result/QC/seurat_object_sep25/"
if (!dir.exists(output_dir)) dir.create(output_dir)

# QC loop
for (i in seq_along(seurat_list)) {
  seurat <- seurat_list[[i]]
  sample_name <- names(seurat_list)[i]
  seurat$badquality_count_per= (seurat$nCount_BlankCodeword+seurat$nCount_ControlCodeword+seurat$nCount_ControlProbe)/seurat$nCount_Xenium
  seurat <- subset(seurat, 
                   subset = nFeature_Xenium > 20 & 
                     nCount_Xenium > 50 & 
                     badquality_count_per < 0.05)
  saveRDS(seurat, file = file.path(output_dir, paste0(sample_name, "_seurat_QC.rds")))
  message("Finished QC and saved: ", sample_name)
}
print('done')

