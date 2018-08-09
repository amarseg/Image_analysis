library('tidyverse')
library('clusterProfiler')
source('functions_CPA.R')

#analyse plate data
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

#Load number of cells per filter
cell_number <- read_csv('../Analysed_data/cell_countsImage.csv') %>%
  select(Count_Nuclei, Count_Nuclei_AR_filtered, Count_Nuclei_AR_Solidity_Filtered, FileName_DNA) %>%
  separate(FileName_DNA, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
  mutate(Proportion_AR_filter = Count_Nuclei_AR_filtered/(Count_Nuclei+Count_Nuclei_AR_filtered), 
         Proportion_solidity_filter = Count_Nuclei_AR_Solidity_Filtered/(Count_Nuclei_AR_filtered+ Count_Nuclei_AR_Solidity_Filtered))

plate_name <- c(rep('Plate1',96*9), rep('Plate2',96*9), rep('Plate3',96*9), rep('Plate4',96*8))
cell_number$plate_name <- plate_name

test <- gather(cell_number, key = feature, value = cell_count, -plate_name, -Well) %>%
  filter(feature == 'Proportion_AR_filter' | feature == 'Proportion_solidity_filter')

ggplot(test, aes(group = feature, fill = feature, x = cell_count)) +
  geom_histogram() +
  facet_wrap(~plate_name)

bad_wells <- filter(cell_number, Proportion_solidity_filter < 0.2) %>%
  select(plate_name, Well,Picture) 

#Add mean and sd of randomised control distribution to area plot
median_dist <- draw_control_dist(all_plates, n = 1000, size = 0.25, column = 'AreaShape_Area')

p <- ggplot(all_plates, aes(x = Metadata_Well, y = AreaShape_Area))
p + geom_boxplot() + 
  facet_wrap(~Metadata_Plate_Name) +
  geom_hline(yintercept = mean(median_dist$value)) +
  geom_hline(yintercept = mean(median_dist$value) + 2*sd(median_dist$value), color = 'red') +
  geom_hline(yintercept = mean(median_dist$value) - 2*sd(median_dist$value), color = 'red')

#Filter mutants higher than the mean
hits <- all_plates %>%
  group_by(Metadata_Plate_Name, Metadata_Well) %>%
  summarise(mean_area = mean(AreaShape_Area), id = unique(`Systematic ID`)) %>%
  filter(mean_area > mean(median_dist$value)) %>%
  write_csv('hits.csv')
