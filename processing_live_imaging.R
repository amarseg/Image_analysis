library(tidyverse)
library(gtools)

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

join_all_data <- bind_rows(join_data, join_data2) %>%
  select(-X1) %>%
  separate(Image_name, into = c('Well','Well_n','Position','Z','Time','File'), sep = '--') %>%
  select(-File) %>%
  mutate_at(vars(Well_n:Time), str_sub,start = 4, end = 6) %>%
  mutate_at(vars(Well_n:Time), as.numeric)

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

strain_data <- read_csv('live_cell_imaging_cp/out.csv') %>%
  select(`Systematic ID`,`96_well`) %>%
  add_column(plate = 'plate1') 

strain_data <- strain_data[mixedorder(strain_data$`96_well`),] %>%
  rownames_to_column()

strain_data2 <- read_csv('live_cell_imaging_cp/out2.csv') %>%
  select(`Systematic ID`,`96_well`) %>%
  add_column(plate = 'plate2')
  
strain_data2 <- strain_data2[mixedorder(strain_data2$`96_well`),] %>%
  rownames_to_column()

strain_data$rowname <- as.numeric(strain_data$rowname)
strain_data2$rowname <- as.numeric(strain_data2$rowname)

join_strain_data <- bind_rows(strain_data, strain_data2)

####################Join everything together#############################

strain_time <- right_join(join_strain_data, join_time_stamp, by = c('plate', 'rowname' = 'W') )

all_data <- left_join(join_all_data, strain_time, by = c('Well_n'='rowname',
                                                         'Position' = 'P',
                                                         'Time' = 'T',
                                                         'plate')) %>%
  mutate(aspect_ratio = Major_axis_length/Minor_axis_length)

all_data[rowSums(is.na(all_data))>1,]$`Systematic ID` <- 'wt'
write_csv(all_data, 'all_live_imaging_data.csv')




