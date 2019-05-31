rm(list = ls())
library(tidyverse)

#Create file ID using the count_Image.csv file
file_id <- read_csv('../live_imaging_second_round/output_cell_profiler/cell_countsImage.csv') 

file_list <- c('../live_imaging_second_round//output_cell_profiler/row_1_april_19_001/Nuclei_AR_Solidity_filtered.csv',
               '../live_imaging_second_round//output_cell_profiler/row_2_april_19_001/Nuclei_AR_Solidity_filtered.csv',
               '../live_imaging_second_round//output_cell_profiler/row_3_april_19_002/Nuclei_AR_Solidity_filtered.csv',
               '../live_imaging_second_round//output_cell_profiler/row_4_april_19_001/Nuclei_AR_Solidity_filtered.csv',
               '../live_imaging_second_round//output_cell_profiler/row_5_april_19_001/Nuclei_AR_Solidity_filtered.csv')



time_stamps <- c('../live_imaging_second_round/row_1_april_19_001/AcquisitionLog.dat',
                 '../live_imaging_second_round/row_2_april_19_001/AcquisitionLog.dat',
                 '../live_imaging_second_round/row_3_april_19_002/AcquisitionLog.dat',
                 '../live_imaging_second_round/row_4_april_19_001/AcquisitionLog.dat',
                 '../live_imaging_second_round/row_5_april_19_001/AcquisitionLog.dat')

time_thr <- 30
all_data <- list()
i <- 1
for(file in file_list)
{
  row <- gsub(".*[/]([^.]+)[/].*", "\\1", file)
  
  time_st <-read_delim(time_stamps[i], delim = '\t', skip = 1, col_names = c('IMAGEPOS','W','P','T','Z','IMAGEX','IMAGEY','IMAGEZ','TIME')) %>%
    select(-IMAGEPOS)
  
  time_st <- apply(time_st, c(1,2), function(x){str_split(x,pattern = '=')[[1]][2]})
  
  time_st <- apply(time_st, c(1,2), as.numeric) %>%
    as.data.frame()
    
  all_data[[i]] <- read_csv(file = file) %>%
    add_column(filename = row) %>%
    left_join(file_id, by = c('filename' = 'Metadata_Plate_Name', 'ImageNumber')) %>%
    separate(FileName_Outlines, into = c('Well','Well_n','Position','Z_axis','Time','Type'), sep = '--') %>%
    transform(Well_n = as.numeric(str_sub(Well_n, start = 5, end = 7))) %>%
    transform(Position = as.numeric(str_sub(Position, start = 5, end = 7))) %>%
    transform(Z_axis = as.numeric(str_sub(Z_axis, start = 5, end = 7))) %>%
    transform(Time = as.numeric(str_sub(Time, start = 5, end = 7))) %>% 
    left_join(time_st, by = c('Well_n' = 'W',
                                      'Position' = 'P',
                                      'Time' = 'T')) 
 
  i<- i+1 
}

strain_code <- read_csv('second_round_imaging/strain_code_ready_for_use.csv') %>%
  select(-original_well)

data <- bind_rows(all_data) %>%
  select(-Metadata_FileLocation,-Metadata_Folder, -Metadata_Position,
         -Metadata_Well, -Metadata_WellNumber) %>%
  left_join(strain_code, by = c('Well' = 'Real_well',
                               'Metadata_Plate_Name' = 'Plate_name')) %>%
  filter(Time < time_thr) %>%
  write_csv('second_round_imaging/small_tidy_results.csv')