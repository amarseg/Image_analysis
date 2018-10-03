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
    select(`Ver5.0 position`,`Systematic ID`) %>%
    separate(`Ver5.0 position`, into = c('Version','Plate','Well')) %>%
    mutate(Plate = str_extract(Plate, pattern = '[:digit:]{2}')) %>%
    mutate(Plate = paste0('Plate',as.numeric(Plate)), Well = as.numeric(Well)) %>%
    filter(Plate == plate_n)
  
  data <- 
    read_csv(csv_path) %>% 
    mutate(File_name = basename(Metadata_FileLocation)) %>%
    separate(File_name, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
    filter(Metadata_QCFlag != 1) %>%
    left_join(strain_data, by = c( 'Metadata_WellNumber' = 'Well'))
  
  return(data)
}

plot_cell_proportion <- function(cell_count_data, bad_well_thr = 0.2)
{
  proportions <- cell_count_data %>%
    select(Count_Nuclei, Count_Nuclei_AR_filtered, Count_Nuclei_AR_Solidity_Filtered, FileName_DNA, Metadata_Plate_Name) %>%
    separate(FileName_DNA, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
    mutate(Proportion_AR_filter = Count_Nuclei_AR_filtered/(Count_Nuclei+Count_Nuclei_AR_filtered), 
           Proportion_solidity_filter = Count_Nuclei_AR_Solidity_Filtered/(Count_Nuclei_AR_filtered+ Count_Nuclei_AR_Solidity_Filtered))
  
  test <- gather(cell_number, key = feature, value = cell_count, -plate_name, -Well) %>%
    filter(feature == 'Proportion_AR_filter' | feature == 'Proportion_solidity_filter')
  
  ggplot(test, aes(x = Well, group = feature, fill = feature, y = cell_count)) +
    geom_bar( stat = 'identity') +
    facet_wrap(~plate_name)
  
  ggplot(test, aes(group = feature, fill = feature, x = cell_count)) +
    geom_histogram() +
    facet_wrap(~plate_name)
  
  bad_wells <- filter(proportions, Proportion_solidity_filter < 0.2) %>%
    select(plate_name, Well, Picture) %>%
    return()
  
}

draw_control_mean_dist <- function(plates_data, n, size, column)
{
  means <- rerun(n, sample_frac(plates_data[,column], size = size)) %>%
    sapply(function(x){mean(x[[1]])}) %>%
    as_data_frame() %>%
    return()
  
}

