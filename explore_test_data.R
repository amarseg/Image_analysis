library('tidyverse')
source('functions_CPA.R')

plate1_test <- 
  load_CPA_data('../Plate1_20f/test_table','../Plate1_20f/image_well.csv', plate_n = 1) %>%
  filter(Image_Metadata_QCFlag == 0 & Type == 'DAPI.tif')

p <- ggplot(plate1_test, aes(x = Well, y = Nuclei_AR_Solidity_Filtered_AreaShape_Area))
p + geom_boxplot()

p <- ggplot(plate1_test, aes(x = Well))
p + geom_bar()

