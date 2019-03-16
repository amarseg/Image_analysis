library(tidyverse)
library(gtools)

file_list <- c('../olaf_cell_profiler/output/plate1/Nuclei_AR_Solidity_Filtered.csv',
               '../olaf_cell_profiler/output/plate2/Nuclei_AR_Solidity_Filtered.csv')
###Load strain data#########
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

####Load times#######
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


#############Load actual data###########
data <- map(file_list, read_csv) %>%
  reduce(bind_rows) %>%
  mutate(FileName = basename(Metadata_FileLocation))

image_data <- read_csv('../olaf_cell_profiler/output/cell_countsImage.csv') 

all_data <- left_join(data, image_data, by = c('FileName' = 'FileName_DNA',
                                               'Metadata_Plate_Name')) %>%
  filter(Metadata_QCFlag.y == 0) %>%
  mutate(Time = str_extract(FileName, pattern = 'T0000[0-6]')) %>%
  mutate(Time = as.numeric(str_sub(Time, start = 5, end = 6))) %>%
  separate(FileName, into = c('Well','Well_n','Position','Z_axis','Time','Type'), sep = '--') %>%
  transform(Well_n = as.numeric(str_sub(Well_n, start = 5, end = 7))) %>%
  transform(Position = as.numeric(str_sub(Position, start = 5, end = 7))) %>%
  transform(Z_axis = as.numeric(str_sub(Z_axis, start = 5, end = 7))) %>%
  transform(Time = as.numeric(str_sub(Time, start = 5, end = 7))) %>%
  left_join(join_time_stamp, by = c('Well_n' = 'W',
                                    'Position' = 'P',
                                    'Time' = 'T',
                                    'Metadata_Plate_Name' = 'plate')) %>%
  left_join(join_strain_data, by = c('Well' = '96_well',
                                     'Metadata_Plate_Name' = 'plate')) %>%
  filter()

all_data[rowSums(is.na(all_data))>1,]$`Systematic ID` <- 'wt'
write_csv(all_data, 'live_cell_imaging_cp/output_data_cleaned.csv')


###################Plot Image data##########################
proc_image_data <- image_data %>%
  separate(FileName_DNA, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
  select(Well, Picture, ImageQuality_Correlation_DNA_20:Metadata_Plate_Name) %>%
  gather(key = variable, value = parameter_value, - Well, -Metadata_Plate_Name, -Picture)


ggplot(proc_image_data, aes(x = parameter_value )) +
  geom_histogram() + 
  facet_wrap(~variable, scales = 'free') 
