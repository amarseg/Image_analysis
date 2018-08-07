load_CPA_data <- function(per_object_data, per_image_data, plate_n)
{
  require('tidyverse')
  object_data <- read_csv(per_object_data)
  well_data <- 
    read_csv(per_image_data) %>%
    select(ImageNumber, Image_Metadata_QCFlag,Image_FileName_DNA) %>%
    separate(Image_FileName_DNA, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
    transform(Well_n = str_sub(Well_n, start = 5, end = 7))
  
 strain_data <- read_csv('../Imaging_code/library_strains.csv') %>%
   select('Ver5.0 position','Systematic ID') %>%
   separate('Ver5.0 position', into = c('Version','Plate','Well')) %>%
   filter(Plate == plate_n)
  
  all_data <- left_join(object_data, well_data, by = c('ImageNumber' = 'ImageNumber'))
  all_data_strain <- left_join(all_data, strain_data, by = c('Well_n' = 'Well'))
  return(all_data_strain)
}

load_filtered_data <- function(csv_path, plate_n)
{
  require(tidyverse)
  
  strain_data <- read_csv('library_strains.csv') %>%
    select('Ver5.0 position','Systematic ID') %>%
    separate('Ver5.0 position', into = c('Version','Plate','Well')) %>%
    mutate(Plate = str_extract(Plate, pattern = '[:digit:]{2}')) %>%
    mutate(Plate = as.numeric(Plate), Well = as.numeric(Well)) %>%
    filter(Plate == plate_n)
  
  data <- 
    read_csv(csv_path) %>% 
    filter(Metadata_QCFlag != 1) %>%
    left_join(strain_data, by = c( 'Metadata_WellNumber' = 'Well'))
  
  return(data)
}