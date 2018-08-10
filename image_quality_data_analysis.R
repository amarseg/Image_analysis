##Analysis of image features
library(tidyverse)
source('functions_CPA.R')

image_data <- read_csv('../Analysed_data/adding_quality_features/cell_countsImage.csv') %>%
  separate(FileName_DNA, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
  select(Well, Picture, ImageQuality_Correlation_DNA_20:Metadata_Plate_Name) %>%
  gather(key = variable, value = parameter_value, - Well, -Metadata_Plate_Name, -Picture)


ggplot(image_data, aes(x = parameter_value )) +
  geom_histogram() + 
  facet_wrap(~variable, scales = 'free') 

ggsave('figs/Distribution_of_quality_features.pdf')

low_power_log <- filter(image_data, variable == 'ImageQuality_PowerLogLogSlope_DNA' & parameter_value < -2.5)
ggplot(low_power_log, aes(x = Metadata_Plate_Name)) +
  geom_bar()

high_correlation <- filter(image_data, variable == 'ImageQuality_Correlation_DNA_20' & parameter_value > 0.5)
ggplot(high_correlation, aes(x = Metadata_Plate_Name)) +
  geom_bar()


high_correlation <- filter(image_data, variable == 'ImageQuality_Correlation_DNA_20' & parameter_value > 0.5
                           | variable == 'ImageQuality_PowerLogLogSlope_DNA' & parameter_value < -2.5 )
ggplot(high_correlation, aes(x = Metadata_Plate_Name)) +
  geom_bar()

