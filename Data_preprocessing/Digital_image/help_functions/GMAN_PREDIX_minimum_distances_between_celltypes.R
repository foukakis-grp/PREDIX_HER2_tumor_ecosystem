GMAN_PREDIX_minimum_distances_between_celltypes <- function(spe_object, cell_types_of_interest, feature_colname ) {
  if(  sum(cell_types_of_interest%in%names(table(spe_object$Phenotype))) !=2   )
  {local_dist_mins <- setNames(data.frame(matrix(ncol = 9, nrow = 0)),c('RefCell','RefType','NearestCell','NearestType','Distance','Pair'))}else{
    local_dist_mins <- calculate_minimum_distances_between_celltypes(spe_object, feature_colname , cell_types_of_interest )}
  
  return(local_dist_mins)
}
    
    

