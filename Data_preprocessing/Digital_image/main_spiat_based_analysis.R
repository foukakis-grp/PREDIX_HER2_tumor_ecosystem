library("tiff")                   
library('tidyr')   
library('plyr')    
library('dplyr')    
library('data.table') 
library('stringr')    
library('writexl')   
library('openxlsx')
library('raster')   
library('splitstackshape')  
library('reshape2')    
library('dbscan')    
library('terra')
library('ggplot2')              
library('tibble')      
library('haven')        
library('SPIAT')      
library('ggthemes') 
library('spatstat')  
library('parallel')
library('foreach')
library(dplyr)
library(stringr)
#---------------------------------------------------------------------------------------------------------   
parallel::detectCores()
n.cores <- parallel::detectCores() - 2
my.cluster <- parallel::makeCluster(n.cores, type = "PSOCK")                          ##### PSOCK for clusters   FORK
clusterEvalQ(my.cluster, c(library(dplyr), library(SPIAT) , library(data.table) , library(stringr) , library(raster) ))
#---------------------------------------------------------------------------------------------------------
phenotypes_string <- tolower(c('ImmuneCells','Stroma'    ,'Tumor' ))
cols <- c( 'ImmuneCells' = 'cyan3' ,  'Stroma' = 'darkgreen'  ,   'Tumor'  = 'red' )
pair_combinations <- c( 'Tumor__ImmuneCells'  ,   'Tumor__Stroma')

#################################################################
data_root <- '/mimer/....'
path <- '/mimer/....'
subpath_to_save <- '.....'
dir.create(file.path(path,'spatial_results', subpath_to_save), showWarnings = FALSE)
setwd(path)
#---------------------------------------------------------------------------------------------------------

source(paste0(path,'help_functions/','GMAN_format_inform_to_spe_HandE.R'))
source(paste0(path,'help_functions/','GMAN_main_summary_min_distance.R'))
source(paste0(path,'help_functions/','GMAN_PREDIX_minimum_distances_between_celltypes.R'))
#---------------------------------------------------------------------------------------------------------

get_field_info = function(path) {
  info = readTIFFDirectory(path, all=FALSE)
  
  required_attributes = c('width', 'length', 'x.resolution', 'resolution.unit')
  missing_attributes = setdiff(required_attributes, names(info))
  if (length(missing_attributes) > 0) {
    missing = paste(missing_attributes, collapse=', ')
    stop(paste0('Image file is missing required attributes: ', missing))
  }
  
  if (info$resolution.unit != 'cm')
    stop(paste('Unsupported resolution unit:', info$resolution.unit))
  
  result = list()
  result$image_size = c(info$width, info$length)
  result$microns_per_pixel = as.numeric(10000/info$x.resolution)
  result$pixels_per_micron = 1/result$microns_per_pixel
  result$field_size = result$image_size * result$microns_per_pixel
  
  # Location directly from TIFF info
  result$location = c(info$x.position, info$y.position) * 10000
  result
}
#---------------------------------------------------------------------------------------------------------
quiet_basic <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) }
#---------------------------------------------------------------------------------------------------------
select_and_transpose <- function(input_df )
{
  input_df <- input_df[,colnames(input_df)%in%c('Cell_type','Proportion'), drop= FALSE]
  input_df <- data.table::transpose(input_df )
  colnames(input_df) <- input_df[1,]
  input_df<- input_df[-1,]
  
  input_df <- input_df%>%mutate_if(is.character, as.numeric)
  return(input_df)}
#---------------------------------------------------------------------------------------------------------
is_empty <- function(x) {
  if (length(x) == 0 & !is.null(x)) {
    TRUE
  } else {
    FALSE
  }}
#---------------------------------------------------------------------------------------------------------
get_colData <- function(spe_object){
  formatted_data <- data.frame(SummarizedExperiment::colData(spe_object))
  formatted_data <- cbind(formatted_data, 
                          data.frame(SpatialExperiment::spatialCoords(spe_object)))
  if (is.null(formatted_data$Cell.ID)){
    formatted_data <- formatted_data %>% tibble::rownames_to_column("Cell.ID")
  }
  
  # shouldn't delete column `sample_id`
  # formatted_data$sample_id <- NULL
  
  return(formatted_data)
}
#---------------------------------------------------------------------------------------------------------
return_files_info <- function(data_root )
{
  all_files <- list.files(data_root , pattern = "_classification.csv$", recursive = TRUE, full.names = T)
  
  pixel_size <-  list.files(data_root , pattern = "pixel_size.csv$", recursive = TRUE, full.names = T)
  
  data_list <- vector(mode = "list", length = length(all_files))
  name_list <- vector(mode = "list", length = length(all_files))
  
  for (i in 1:length(all_files))
  {
    pixel_size_i <- read.csv(pixel_size[i], header = TRUE)
    
    tmp <- read.csv(all_files[i], header = TRUE)
    
    if(nrow(tmp)>0){
      
      print(all_files[i])
      
      tmp <- tmp[!grepl("Other", tmp$class),]
      tmp <- tmp%>%drop_na(class)
      
      tmp <- tmp[!(is.na(tmp$class) | tmp$class==""), ]
      
      tmp <- tmp %>%
        mutate(across('class', str_replace, 'Immune cells', 'ImmuneCells'))

      names(tmp)[names(tmp)%in%c('cell_centroid_x')] <- "Cell.X.Position"
      names(tmp)[names(tmp)%in%c('cell_centroid_y')] <- "Cell.Y.Position"
      names(tmp)[names(tmp)%in%c('class')] <- "Cell.Type.Basic"
      tmp$Phenotype <- tmp$Cell.Type.Basic
      
      tmp$Cell.X.Position <- tmp$Cell.X.Position*pixel_size_i$pixelWidth 
      tmp$Cell.Y.Position <- tmp$Cell.Y.Position*pixel_size_i$pixelHeight
      data_list[[i]] <-  tmp
      name_list[[i]] <- all_files[i]
    }
  }
  
  data_list[sapply(data_list, is.null)] <- NULL
  name_list[sapply(name_list, is.null)] <- NULL
  
  return(list(data_list , name_list))
}
#---------------------------------------------------------------------------------------------------------
GMAN_PREDIX_average_percentage_of_cells_within_radius <- function(spe_object, reference_celltype, target_celltype, radius_vect , feature_colname,PopulationLabel ) {
  
  local_dist_mins <- setNames(data.frame(matrix(ncol = length(radius_vect), nrow = 1)),c(paste0(reference_celltype , '_to_' ,  target_celltype , '_Rad_',radius_vect)))
  rownames(local_dist_mins) <- paste0('pat_' , unique(spe_object$patient_id))
  
  if(  sum(    c(reference_celltype, target_celltype)%in%names(table(spe_object$Phenotype))    ) ==2   )
  {for(ii in 1:length(radius_vect)){local_dist_mins[ii] <- average_percentage_of_cells_within_radius(spe_object, reference_celltype, target_celltype, radius_vect[ii] , feature_colname )}}
  
  colnames(local_dist_mins)  <- paste0( 'PercentCells__' , PopulationLabel , '__'  ,  colnames(local_dist_mins)  )
  
  return(local_dist_mins)
}
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
GMAN_calculate_pairwise_distances_between_celltypes <- function(spe_data ,  refer_target_names , feature_colname ,  PopulationLabel){
  
  cell_types_of_interest <- c(str_split_fixed(refer_target_names, "__", 2)[1], str_split_fixed(refer_target_names, "__", 2)[2])
  tmp <- setNames(data.frame(matrix(ncol = 9, nrow = 4)), c('Pair','Mean' ,'Min', 'Max','Median','Std.Dev','Reference','Target','PatName'))
  
  if(  sum(    cell_types_of_interest%in%names(table(spe_data$Phenotype))    ) ==2   ){
    
    
    tmp2<-  calculate_pairwise_distances_between_celltypes(spe_data, cell_types_of_interest , feature_colname)
    tmp2 <- calculate_summary_distances_between_celltypes(tmp2)
    tmp2 <- cbind.data.frame(  tmp2  ,   PatName = unique(spe_data$imageName) )  
    
    
    #----- if calculate_summary_distances_between_celltypes return less than 4
    if (nrow(tmp2)<nrow(tmp)){tmp[1:nrow(tmp2),] <- tmp2}else{tmp <- tmp2}
    tmp$PatName <- unique(spe_data$imageName)
    
  }else{
    
    tmp$Reference <- c(  cell_types_of_interest[[1]]  ,  cell_types_of_interest[[2]]  , cell_types_of_interest[[1]]  , cell_types_of_interest[[2]] )
    tmp$Target <-    c(  cell_types_of_interest[[1]]  ,  cell_types_of_interest[[2]]  , cell_types_of_interest[[2]]  , cell_types_of_interest[[1]] )
    tmp$PatName <- unique(spe_data$imageName)
  }
  
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
  
  colnames(tmp)  <- paste0( 'PairDists__' , PopulationLabel , '__'  ,  colnames(tmp)  )
  
  return(tmp)
}
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
return_phenotypes_ifexist <- function(sample_metadata, phenotypes_string,cols)
{
  sample_metadata <- as.data.frame(sample_metadata)
  ret_cols <- data.frame()
  ret_pheno <- c()
  ret_sample_metadata <- data.frame()
  
  for (i in 1:length(phenotypes_string))
  {
    if(sum( tolower(sample_metadata$Phenotype) == phenotypes_string[i] ) > 0)
    {
      tmp <- sample_metadata[ tolower(sample_metadata$Phenotype)==phenotypes_string[i], ]
      tmp$Phenotype <- names(cols)[i]
      
      ret_sample_metadata <- rbind.data.frame( ret_sample_metadata  ,  tmp )
      
      ret_cols <- c( ret_cols , cols[i] )
      ret_pheno <- c( ret_pheno , names(cols)[i] )
    }
  }
  
  return( list(ret_sample_metadata = ret_sample_metadata , ret_cols = unlist(ret_cols)  ,  ret_pheno = ret_pheno , 
               patient_name = unique(sample_metadata$imageName)  , x_max_distance = (max(sample_metadata$Cell.X.Position) - min(sample_metadata$Cell.X.Position))  ,
               y_max_distance = (max(sample_metadata$Cell.Y.Position) - min(sample_metadata$Cell.Y.Position))    ) )
}
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
return_pairs_ifexist <- function(data_list, pair_combinations)
{
  dist_df <- data.frame(row.names=1)
  dist_df_names <- c()
  
  pheno_exist <- unique(data_list$Phenotype)
  
  for (i in 1:length(pair_combinations))
  {
    tmp <- str_split_1(pair_combinations[i], "__")
    
    #if pair exists do the math
    if(  tmp[1]%in%pheno_exist &  tmp[2]%in%pheno_exist  )
    {
      #first_part for initial pair_combination. second_part for reversed
      
      first_part  <- as.data.frame(data_list %>% filter(select_rows(data_list, tmp[1])))
      first_part  <- first_part %>% dplyr::select( paste0('Distance to ',tmp[2]) )
      
      dist_df <- cbind.data.frame( dist_df  , colMeans(first_part,na.rm = TRUE)   )
      dist_df <- cbind.data.frame( dist_df  , median(pull(first_part),na.rm = TRUE)   )
      
      second_part <- data_list %>% filter(select_rows(data_list, tmp[2]))
      second_part  <- second_part %>% dplyr::select( paste0('Distance to ',tmp[1]) )
      
      dist_df <- cbind.data.frame( dist_df  , colMeans(second_part,na.rm = TRUE)   )
      dist_df <- cbind.data.frame( dist_df  , median(pull(second_part),na.rm = TRUE)   )
    }
    #else return NA
    else{
      
      dist_df <- cbind.data.frame( dist_df  , NaN )
      dist_df <- cbind.data.frame( dist_df  , NaN )
      dist_df <- cbind.data.frame( dist_df  , NaN )
      dist_df <- cbind.data.frame( dist_df  , NaN )
      
      
    }
    #fill names
    dist_df_names <- c( dist_df_names , paste0(pair_combinations[i],'__mean') )
    dist_df_names <- c( dist_df_names , paste0(pair_combinations[i],'__median') )
    
    #------------------------------------------------------------------
    #reverse pair_combination name first
    splits <- strsplit(pair_combinations[i], "__")[[1]]
    reversed <- rev(splits)
    final_result <- paste(reversed, collapse = "__")
    
    dist_df_names <- c( dist_df_names , paste0(final_result,'__mean') )
    dist_df_names <- c( dist_df_names , paste0(final_result,'__median') )
    
  }
  
  colnames(dist_df) <- dist_df_names  
  return(dist_df)
}
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
remove_duplicates_per_tile <- function(sce_object){
  #split per images within the same patient
  exam_tmp <- split(sce_object , sce_object$sample)

  #for each image within the patient
  for (ii in 1:length(exam_tmp)){
    #-------------------------------------------------------------
    ######   SOS I removed duplicated rows containing SAME X,Y position !!!!!!!!!
    tmp <- exam_tmp[[ii]]
    exam_tmp[[ii]] <- tmp[    !duplicated(tmp[colnames(tmp)%in%c('Cell.X.Position','Cell.Y.Position') ]),]
    exam_tmp[[ii]]$Cell.ID <- 1:nrow(exam_tmp[[ii]])
    #-------------------------------------------------------------
  }
  
  new_sce_object <- plyr::ldply (exam_tmp, data.frame , .id=NULL  ,  .parallel = TRUE)
  new_sce_object$Cell.ID <- 1:nrow(new_sce_object)
  
  return(new_sce_object)
}
#-----------------------------------------------------------------------------------------------------------------------------------------------------------

clusterExport(my.cluster, c("gman_average_nearest_neighbor_index",'get_colData','is_empty','quiet_basic','calculate_cell_proportions','GMAN_main_summary_min_distance',
                            'GMAN_PREDIX_average_percentage_of_cells_within_radius'), envir=environment())
doParallel::registerDoParallel(cl = my.cluster)

#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#STEP1: get values from csv
super_info <- return_files_info(data_root)
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#STEP2: generate spe from df
super_info_spe_all <- parLapply(cl = my.cluster , seq_along(super_info[[1]])  ,  GMAN_format_inform_to_spe_HandE  ,  super_info[[1]] ,   tissue_ref =  'Entire'  ,  pheno_ref = 'ALL'  , 
                                tissue_target = 'Entire' , pheno_target = 'ALL' , ref_advanced = FALSE , target_advanced = FALSE)


unique_pheno_names <- lapply( seq_along(super_info_spe_all)  ,  function(x,super_info_spe_all ){if(!(sum(is.na(super_info_spe_all[[x]])))){unique(super_info_spe_all[[x]]$Phenotype)}} , super_info_spe_all)
unique_pheno_names <- unique(unlist(unique_pheno_names))

analysis_is <- 'HandE'
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#STEP3: Analysis
spatial_1_celldist <- NA

#--------------- ALL
spatial_1_celldist <- parLapply(cl = my.cluster , super_info_spe_all  ,  calculate_cell_proportions  ,  reference_celltypes = NULL,
                                feature_colname ="Phenotype",celltypes_to_exclude = "Others",plot.image = FALSE)

df_names <- lapply( seq_along(super_info_spe_all)  ,  function(x,super_info_spe_all ){unique(super_info_spe_all[[x]]$imageName)} , super_info_spe_all)

spatial_1_celldist <- parLapply(cl = my.cluster , spatial_1_celldist  ,  select_and_transpose )
spatial_1_celldist <- rbind.fill(spatial_1_celldist)
rownames(spatial_1_celldist) <- paste0('pat_' , unlist(df_names))
colnames(spatial_1_celldist) <- paste0( 'CellProps__ALL__'  , colnames(spatial_1_celldist) )
print('***************   calculate_cell_proportions ENDS **********')


for (i in 1:length(pair_combinations))
{all_mindist <- cbind.data.frame ( all_mindist  , GMAN_main_summary_min_distance(pair_combinations[i] , my.cluster , super_info_spe_all , 'ALL') )}
all_mindist<- all_mindist[,-1]
SPATIAL_2_MIN_DISTS <-   all_mindist
SPATIAL_2_MIN_DISTS <- SPATIAL_2_MIN_DISTS[,colSums(is.na(SPATIAL_2_MIN_DISTS))<nrow(SPATIAL_2_MIN_DISTS)]


