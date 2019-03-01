library(tidyverse)

data_path <- 'Z:/Pers_Amalia/Screening/live_imaging/processed_data/plate1'
file_names <- list.files(data_path, recursive = T, full.names = T) %>%
  grep(pattern = '.csv', fixed = T, value = T)

join_data <- file_names %>%
  map(read_csv) %>%
  reduce(bind_rows) %>%
  add_column(plate = 'plate1')

path2 <- 'Z:/Pers_Amalia/Screening/live_imaging/processed_data/plate2'
file_names2 <- list.files(path2, recursive = T, full.names = T) %>%
  grep(pattern = '.csv', fixed = T, value = T)

join_data2 <- file_names2 %>%
  map(read_csv) %>%
  reduce(bind_rows) %>%
  add_column(plate = 'plate2')

#########Load time data############
time_stamp <- read_delim('AcquisitionLog.dat', delim = '\t', skip = 1, col_names = c('IMAGEPOS',
                                                                                     'W','P','T','Z',
                                                                                     'IMAGEX','IMAGEY','IMAGEZ','TIME')) %>%
  select(-IMAGEPOS)

time_stamp <- apply(time_stamp, c(1,2), function(x){str_split(x,pattern = '=')[[1]][2]})
time_stamp <- apply(time_stamp, c(1,2), as.numeric) %>%
  as.data.frame() %>%
  add_column(plate = 'plate1')

time_stamp2 <- read_delim('AcquisitionLog2.dat', delim = '\t', skip = 1, col_names = c('IMAGEPOS',
                                                                                     'W','P','T','Z',
                                                                                     'IMAGEX','IMAGEY','IMAGEZ','TIME')) %>%
  select(-IMAGEPOS) 

time_stamp2 <- apply(time_stamp2, c(1,2), function(x){str_split(x,pattern = '=')[[1]][2]})
time_stamp2<- apply(time_stamp2, c(1,2), as.numeric) %>%
  as.data.frame() %>%
  add_column(plate = 'plate2')

join_time_stamp <- bind_rows(time_stamp, time_stamp2)

################Load strain data#################