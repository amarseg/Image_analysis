library('tidyverse')
source('functions_CPA.R')

all <- 
  load_CPA_data('../test_table','../image_well.csv') %>%
  filter(Image_Metadata_QCFlag == 0 & Type == 'DAPI.tif')

p <- ggplot(all, aes(x = Well, y = Nuclei_AR_Solidity_Filtered_AreaShape_Area))
p + geom_boxplot()

p <- ggplot(all, aes(x = Well))
p + geom_bar()
