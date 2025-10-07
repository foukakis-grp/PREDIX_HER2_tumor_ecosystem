
#######  SOS to eortc_data NA VGAINEI META TO RUN TOU ----> TME_level_analysis_merged_data2.R
library('ggplot2')
library('tidyverse')
library('haven')
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
library(foreach)
library(doParallel)


is_empty <- function(x) (!(nrow(x)==0 || ncol(x) ==0))

cols <- c('Tumor'='red' ,'Stroma'='darkgreen' ,'Immune cells'='blue')

path_predix_hande_main <- '/mimer/.....'         #
path_predix_hande_specific <- 'predix-2/output'
new_path_to_write <- '/mimer/....'

#---------------------------------------------------------------------------------------------------------
parallel::detectCores()
n.cores <- parallel::detectCores() - 4
my.cluster <- parallel::makeCluster(n.cores, type = "FORK")                          ##### PSOCK for clusters   FORK
clusterEvalQ(my.cluster, c(library(dplyr), library(data.table) , library(stringr) ))
clusterExport(my.cluster, c("is_empty"), envir=environment())
#---------------------------------------------------------------------------------------------------------

TILE_SIZE <- c(1000,2000)

dir.create(file.path(new_path_to_write, strsplit(path_predix_hande_specific,'/')[[1]][1] ), showWarnings = FALSE)
dir.create(file.path(new_path_to_write, paste0(path_predix_hande_specific,'_READY')), showWarnings = FALSE)
patient_folders <- list.files(file.path(path_predix_hande_main , path_predix_hande_specific),full.names = TRUE)
options(digits=20)
#---------------------------------------------------------------------------------------------------------

predix_prepare_data <- function(  patient_folders_ ,  new_path_to_write  ,  path_predix_hande_specific  , TILE_SIZE  )
{
i <- patient_folders_
if(file.exists(   file.path(i,'tiles')    )){
  
dir.create(file.path(new_path_to_write, paste0(path_predix_hande_specific,'_READY'),basename(i)), showWarnings = FALSE)
patient_files <- list.files(i)#,pattern = "\\.csv$")
#-----------------------------------------------------------------------------------------------------------------------------------------------------------
#-------------------  classification file -------------------------------------------------------------------------------

classif_csv <- read.csv( file.path( i , grep('_classification', patient_files, value=TRUE) ) ,header = TRUE, sep = ",")


#-------------------  pixel size file -------------------------------------------------------------------------------
pixel_size_csv <- read.csv( file.path( i , grep('pixel_size', patient_files, value=TRUE) ) ,header = TRUE, sep = ",")

#-------------------  cell y coords file -------------------------------------------------------------------------------
Ncol <- max(unlist(lapply(strsplit(readLines(file.path( i , grep('cell_polygons_y', patient_files, value=TRUE) )), ","), length)))
cell_polygons_y_csv <- read.csv( file.path( i , grep('cell_polygons_y', patient_files, value=TRUE) ) ,header = F, 
                                 dec = "." , sep = ",", skip = 1, col.names=paste0("V", 1:Ncol))

#-------------------  cell y coords file CYTOPLASM-------------------------------------------------------------------------------
Ncol <- max(unlist(lapply(strsplit(readLines(file.path( i , grep('cytoplasm_polygons_y', patient_files, value=TRUE) )), ","), length)))
cell_polygons_y_CYTO_csv <- read.csv( file.path( i , grep('cytoplasm_polygons_y', patient_files, value=TRUE) ) ,header = F, 
                                      dec = "." , sep = ",", skip = 1, col.names=paste0("V", 1:Ncol))

#-------------------  cell x coords file -------------------------------------------------------------------------------
Ncol <- max(unlist(lapply(strsplit(readLines(file.path( i , grep('cell_polygons_x', patient_files, value=TRUE) )), ","), length)))
cell_polygons_x_csv <- read.csv( file.path( i , grep('cell_polygons_x', patient_files, value=TRUE) ) ,header = F, 
                                 dec = "." , sep = ",", skip = 1, col.names=paste0("V", 1:Ncol))

#-------------------  cell x coords file CYTOPLASM-------------------------------------------------------------------------------
Ncol <- max(unlist(lapply(strsplit(readLines(file.path( i , grep('cytoplasm_polygons_x', patient_files, value=TRUE) )), ","), length)))
cell_polygons_x_CYTO_csv <- read.csv( file.path( i , grep('cytoplasm_polygons_x', patient_files, value=TRUE) ) ,header = F, 
                                      dec = "." , sep = ",", skip = 1, col.names=paste0("V", 1:Ncol))

rownames(classif_csv) <- rownames(cell_polygons_y_csv) <- rownames(cell_polygons_x_csv) <- paste0('cell',1:nrow(cell_polygons_x_csv))
rownames(cell_polygons_y_CYTO_csv) <- rownames(cell_polygons_x_CYTO_csv) <- paste0('cell',1:nrow(cell_polygons_x_csv))

# remove "other" Cells
classif_csv <- classif_csv[!str_detect(tolower(classif_csv$class), "other"), ]  

#remove TILs having less than 20 cells!!
thres_cells <- 20

classif_csv_BU <- classif_csv
cell_polygons_y_csv_BU <- cell_polygons_y_csv
cell_polygons_x_csv_BU <- cell_polygons_x_csv 
cell_polygons_y_CYTO_csv_BU <- cell_polygons_y_CYTO_csv 
cell_polygons_x_CYTO_csv_BU <- cell_polygons_x_CYTO_csv

for (tile_num in TILE_SIZE){
  
  classif_csv <- classif_csv_BU
  cell_polygons_y_csv <- cell_polygons_y_csv_BU
  cell_polygons_x_csv <- cell_polygons_x_csv_BU 
  cell_polygons_y_CYTO_csv <- cell_polygons_y_CYTO_csv_BU 
  cell_polygons_x_CYTO_csv <- cell_polygons_x_CYTO_csv_BU
  
  required_column <- colnames(classif_csv)[grepl('tile_name' , colnames(classif_csv)) & grepl(tile_num , colnames(classif_csv))]
  classif_csv$fake_var <- classif_csv[[required_column]]
  
  classif_csv <- classif_csv[!with(classif_csv,is.na(classif_csv$fake_var)),]
  classif_csv <- classif_csv[!with(classif_csv,!nzchar(classif_csv$fake_var)),]
  
  
  stats_per_til_freq <- dcast(classif_csv, fake_var ~ annID, fun.aggregate = length)
  TILs_ID <- stats_per_til_freq$fake_var
  
  TILs_ID <- TILs_ID[nzchar(TILs_ID)]
  
  stats_per_til_freq <- stats_per_til_freq[,!colnames(stats_per_til_freq)%in%c('fake_var'), drop=FALSE]
  
  for (ii in 1:length(TILs_ID)){
    for (j in 1:ncol(stats_per_til_freq))
      if(stats_per_til_freq[ii,j]<thres_cells)
      {classif_csv <- classif_csv[ ! (classif_csv$fake_var==TILs_ID[ii]  &  classif_csv$annID%in%colnames(stats_per_til_freq)[j] ) ,  ]}}
  
  cell_polygons_x_csv <- cell_polygons_x_csv[ rownames(cell_polygons_x_csv)%in%rownames(classif_csv),  ]
  cell_polygons_x_CYTO_csv <- cell_polygons_x_CYTO_csv[ rownames(cell_polygons_x_CYTO_csv)%in%rownames(classif_csv),  ]
  
  cell_polygons_y_csv <- cell_polygons_y_csv[ rownames(cell_polygons_y_csv)%in%rownames(classif_csv),  ]
  cell_polygons_y_CYTO_csv <- cell_polygons_y_CYTO_csv[ rownames(cell_polygons_y_CYTO_csv)%in%rownames(classif_csv),  ]
  
  classif_csv$CellInfo <- cell_polygons_y_csv$CellInfo <- cell_polygons_x_csv$CellInfo <- rownames(classif_csv)
  cell_polygons_y_CYTO_csv$CellInfo <- cell_polygons_x_CYTO_csv$CellInfo <- rownames(classif_csv)
  
  #split into tils 
  classif_csv_split <- split(x = classif_csv, f = ~ annID + fake_var , sep = "___")
  classif_csv_split <- classif_csv_split[sapply(classif_csv_split, is_empty)] 
  
  dirs = list.dirs(file.path(i))
  #find correct paths for each zoom
  dirs  <- dirs[grepl(paste0(tile_num,'x',tile_num) , dirs)]
  
  for (ii in names(classif_csv_split)){
    
    print(ii)
    
    tiles_directory <- dirs[which(grepl(  strsplit(ii,'___')[[1]][1] ,dirs ))]
    
    tiles_files <- list.files(tiles_directory,pattern = "\\.png$")

    tile_to_move <- tiles_files[ which(   strsplit(ii,'___')[[1]][2]  == tiles_files   )]
    
    new_name_tile <- paste0(  'ZOOM' ,  tile_num , '_' ,  ii   )
    
    file.copy( file.path(tiles_directory ,tile_to_move ) ,  
               file.path(new_path_to_write, paste0(path_predix_hande_specific,'_READY'),basename(i)),
               overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)
    
    file.rename(file.path(new_path_to_write, paste0(path_predix_hande_specific,'_READY'),basename(i),tile_to_move )    ,
                file.path(new_path_to_write, paste0(path_predix_hande_specific,'_READY'),basename(i),new_name_tile) )
    
    #--------------------------------
    cell_polygons_x_csv_split <- cell_polygons_x_csv[ rownames(cell_polygons_x_csv)%in%rownames(classif_csv_split[[ii]]),  ] 
    cell_polygons_x_csv_split[is.na(cell_polygons_x_csv_split)] <- -9999
    
    cell_polygons_y_csv_split <- cell_polygons_y_csv[ rownames(cell_polygons_y_csv)%in%rownames(classif_csv_split[[ii]]),  ] 
    cell_polygons_y_csv_split[is.na(cell_polygons_y_csv_split)] <- -9999
    #--------------------------------
    cell_polygons_x_CYTO_csv_split <- cell_polygons_x_CYTO_csv[ rownames(cell_polygons_x_CYTO_csv)%in%rownames(classif_csv_split[[ii]]),  ] 
    cell_polygons_x_CYTO_csv_split[is.na(cell_polygons_x_CYTO_csv_split)] <- -9999
    
    cell_polygons_y_CYTO_csv_split <- cell_polygons_y_CYTO_csv[ rownames(cell_polygons_y_CYTO_csv)%in%rownames(classif_csv_split[[ii]]),  ] 
    cell_polygons_y_CYTO_csv_split[is.na(cell_polygons_y_CYTO_csv_split)] <- -9999
    
    new_name_tile2 <-  strsplit(new_name_tile,'.png')[[1]][1]
    
    
    write.table(classif_csv_split[[ii]], file.path(new_path_to_write, paste0(path_predix_hande_specific,'_READY'),basename(i), paste0(new_name_tile2,'.csv') ), sep='\t',row.names = FALSE, col.names = TRUE )
    
    write.table(cell_polygons_x_csv_split, file.path(new_path_to_write, paste0(path_predix_hande_specific,'_READY'),basename(i), paste0(new_name_tile2,'_PolygX.csv') ), sep='\t',row.names = FALSE, col.names = TRUE )
    write.table(cell_polygons_y_csv_split, file.path(new_path_to_write, paste0(path_predix_hande_specific,'_READY'),basename(i), paste0(new_name_tile2,'_PolygY.csv') ), sep='\t',row.names = FALSE, col.names = TRUE )
    
    write.table(cell_polygons_x_CYTO_csv_split, file.path(new_path_to_write, paste0(path_predix_hande_specific,'_READY'),basename(i), paste0(new_name_tile2,'_CYTOPolygX.csv') ), sep='\t',row.names = FALSE, col.names = TRUE )
    write.table(cell_polygons_y_CYTO_csv_split, file.path(new_path_to_write, paste0(path_predix_hande_specific,'_READY'),basename(i), paste0(new_name_tile2,'_CYTOPolygY.csv') ), sep='\t',row.names = FALSE, col.names = TRUE )
    
    write.table(pixel_size_csv, file.path(new_path_to_write, paste0(path_predix_hande_specific,'_READY'),basename(i), paste0('WSI_info.csv') ), sep='\t',row.names = FALSE, col.names = TRUE )
  }
}
}
}
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################
########################################################################################################################################

parLapply(cl = my.cluster , patient_folders  ,  predix_prepare_data, new_path_to_write  ,  path_predix_hande_specific   , TILE_SIZE)






  