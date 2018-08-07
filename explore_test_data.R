library('tidyverse')
source('functions_CPA.R')


plate1 <- load_filtered_data('../Analysed_data/Plate1_20f/Nuclei_AR_Solidity_Filtered.csv', plate_n = 1)
plate2 <- load_filtered_data('../Analysed_data/Plate2_20f/Nuclei_AR_Solidity_Filtered.csv', plate_n = 2)
plate3 <- load_filtered_data('../Analysed_data/Plate3_20f/Nuclei_AR_Solidity_Filtered.csv', plate_n = 3)
plate4 <- load_filtered_data('../Analysed_data/Plate4_wPI/Nuclei_AR_Solidity_Filtered.csv', plate_n = 4)

all_plates <- bind_rows(plate1, plate2, plate3, plate4)
p <- ggplot(all_plates, aes(x = Metadata_Well, y = AreaShape_Area))
p + geom_boxplot()+ facet_wrap(~Metadata_Plate_Name)

p <- ggplot(all_plates, aes(x = Metadata_Well))
p + geom_bar()+ facet_wrap(~Metadata_Plate_Name)

p <- ggplot(all_plates, aes(x = Metadata_Well, y = AreaShape_MinorAxisLength))
p + geom_boxplot()+ facet_wrap(~Metadata_Plate_Name)
