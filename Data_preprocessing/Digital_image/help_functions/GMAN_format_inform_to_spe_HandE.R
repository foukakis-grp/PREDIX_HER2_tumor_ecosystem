GMAN_format_inform_to_spe_HandE <- function( count_i  ,spe_object , tissue_ref , pheno_ref , tissue_target , pheno_target , ref_advanced , target_advanced){

  spe_object <- data.frame(spe_object[[count_i]])

  if(is_empty(spe_object) | nrow(spe_object)==0 )
  {
    spe<- NA
  }else{
  
  rownames(spe_object) <- paste0('cell_',1:nrow(spe_object))
  
  if(tissue_ref=='Entire'){tissue_ref <- paste('Epithelial' , 'Stroma',sep='__')}
  if(tissue_target=='Entire'){tissue_target <- paste('Epithelial' , 'Stroma',sep='__')}
  
  build_ref <- data.frame()
  build_target <- data.frame()
  
  ###########################    REF
  #all phenotypes
  if(pheno_ref == 'ALL' ){
    build_ref <- spe_object
    if(ref_advanced){build_ref$inform_spatial <- build_ref$Cell.Type.Advanced
    }else{build_ref$inform_spatial <- build_ref$Cell.Type.Basic}
    #phenotype either basic or advanced
  }else{
    if(ref_advanced){
      build_ref <- spe_object[spe_object$Cell.Type.Advanced==pheno_ref,]
      build_ref$inform_spatial <- build_ref$Cell.Type.Advanced
    }else{
      build_ref <- spe_object[spe_object$Cell.Type.Basic==pheno_ref,]
      build_ref$inform_spatial <- build_ref$Cell.Type.Basic}
  }
  if(tissue_ref!='Epithelial__Stroma'){build_ref <- build_ref[ build_ref$Tissue.Category == tissue_ref,]}
  #
  ###########################    TARGET
  #all phenotypes
  if(pheno_target == 'ALL' ){
    build_target <- spe_object
    if(target_advanced){build_target$inform_spatial <- build_target$Cell.Type.Advanced
    }else{build_target$inform_spatial <- build_target$Cell.Type.Basic}
    #phenotype either basic or advanced
  }else{
    if(target_advanced){
      build_target <- spe_object[spe_object$Cell.Type.Advanced==pheno_target,]
      build_target$inform_spatial <- build_target$Cell.Type.Advanced
    }else{
      build_target <- spe_object[spe_object$Cell.Type.Basic==pheno_target,]
      build_target$inform_spatial <- build_target$Cell.Type.Basic}
  }
  if(tissue_target!='Epithelial__Stroma'){build_target <- build_target[ build_target$Tissue.Category == tissue_target,]} 
  
  duprows <- rownames(build_ref) %in% rownames(build_target)
  return_spe <- rbind(build_target, build_ref[!duprows,])
  
  ########## now do the spe!!!! ####################################
  
  if(nrow(return_spe)>0){
    
    #now create the spe object...
    #grab the intensity level, markers and cell IDs
    assay_data <- return_spe
    
    #transpose the matrix so every column is a cell and every row is a marker
    assay_data_matrix <- as.matrix(assay_data)
    colnames(assay_data_matrix) <- NULL
    rownames(assay_data_matrix) <- NULL
    assay_data_matrix_t <- t(assay_data_matrix)
    
    
    #Assign the phenotype, X and Y positions and cell property columns
    metadata_columns <- return_spe[,c("Phenotype", "Cell.X.Position",  "Cell.Y.Position",'imageName')]
    
    spe <- SpatialExperiment::SpatialExperiment(
      assay = assay_data_matrix_t,
      colData = metadata_columns,
      spatialCoordsNames = c("Cell.X.Position", "Cell.Y.Position"))#,
    #  sample_id = as.character(unique(return_spe$patient_id)))
    
    rownames(spe) <- colnames(return_spe)
    colnames(spe) <- 1:nrow(return_spe)#########################################return_spe$Cell.ID
    
  }else{spe<- NA}
  }
  return(spe)
}