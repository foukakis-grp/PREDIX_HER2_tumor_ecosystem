GMAN_main_summary_min_distance <- function(refer_target_names , my.cluster , spe_data , PopulationLabel){
  cell_types_of_interest = c(str_split_fixed(refer_target_names, "__", 2)[1], str_split_fixed(refer_target_names, "__", 2)[2])
  spatial_2_celldist <- parLapply(cl = my.cluster , spe_data  ,  GMAN_PREDIX_minimum_distances_between_celltypes,  
                                  cell_types_of_interest , feature_colname = "Phenotype")
  
  summary_distances_ALL <- data.frame()
  for(i in 1:length(spe_data) ){
    if(nrow(spatial_2_celldist[[i]])<1){
      tmp <- setNames(data.frame(matrix(ncol = 9, nrow = 2)), c('Pair','Mean' ,'Min', 'Max','Median','Std.Dev','Reference','Target','PatName'))
      tmp$Reference <- cell_types_of_interest
      tmp$Target <- rev(cell_types_of_interest)
      tmp$PatName <- unique(spe_data[[i]]$imageName)
    }else{
      tmp <- calculate_summary_distances_between_celltypes(spatial_2_celldist[[i]])
      tmp <- cbind.data.frame(  tmp  ,   PatName = unique(spe_data[[i]]$imageName) )    }
    
    tmp <- tmp[,-1]
    tmp_name <- unique(tmp$PatName)
    tmp <- tmp %>% dplyr::select(-contains(c('PatName')  ))
    tmp$overall <- paste0(  tmp$Reference , '__'  , tmp$Target  )
    tmp <- tmp %>% dplyr::select(-contains(c('Reference','Target')  ))
    tmp <- melt(tmp,'overall')
    tmp$overall <- paste0(  tmp$overall , '__'  , tmp$variable  )
    tmp <- tmp %>% dplyr::select(-contains(c('variable')  ))
    tmp <- data.table::transpose(tmp)
    colnames(tmp) <- tmp[1,]
    tmp <- tmp[-1,]
    rownames(tmp) <- paste0('pat_' , tmp_name)
    summary_distances_ALL <- rbind.data.frame(summary_distances_ALL , tmp)
  }
  
  colnames(summary_distances_ALL)  <- paste0( 'MinDist__' , PopulationLabel , '__'  ,  colnames(summary_distances_ALL)  )
  
  return(summary_distances_ALL)
}