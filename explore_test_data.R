library('tidyverse')
library('clusterProfiler')
source('functions_CPA.R')

#analyse plate data
plate1 <- load_filtered_data('../Analysed_data/adding_quality_features/Plate1_20f/Nuclei_AR_Solidity_Filtered.csv', plate_n = 1)
plate2 <- load_filtered_data('../Analysed_data/adding_quality_features/Plate2_20f/Nuclei_AR_Solidity_Filtered.csv', plate_n = 2)
plate3 <- load_filtered_data('../Analysed_data/adding_quality_features/Plate3_20f/Nuclei_AR_Solidity_Filtered.csv', plate_n = 3)
plate4 <- load_filtered_data('../Analysed_data/adding_quality_features/Plate4_wPI/Nuclei_AR_Solidity_Filtered.csv', plate_n = 4)

all_plates <- bind_rows(plate1, plate2, plate3, plate4)
p <- ggplot(all_plates, aes(x = Metadata_Well, y = AreaShape_Area))
p + geom_boxplot()+ facet_wrap(~Metadata_Plate_Name)

p <- ggplot(all_plates, aes(x = Metadata_Well))
p + geom_bar()+ facet_wrap(~Metadata_Plate_Name)

p <- ggplot(all_plates, aes(x = Metadata_Well, y = AreaShape_MinorAxisLength))
p + geom_boxplot()+ facet_wrap(~Metadata_Plate_Name)

#Load number of cells per filter
cell_number <- read_csv('../Analysed_data/adding_quality_features/cell_countsImage.csv') %>%
  separate(FileName_DNA, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
  mutate(Proportion_AR_filter = Count_Nuclei_AR_filtered/(Count_Nuclei+Count_Nuclei_AR_filtered), 
         Proportion_solidity_filter = Count_Nuclei_AR_Solidity_Filtered/(Count_Nuclei_AR_filtered+ Count_Nuclei_AR_Solidity_Filtered))

test <- gather(cell_number, key = feature, value = cell_count, -Metadata_Plate_Name, -Well) %>%
  filter(feature == 'Proportion_AR_filter' | feature == 'Proportion_solidity_filter')

ggplot(test, aes(group = feature, fill = feature, x = as.numeric(cell_count))) +
  geom_histogram() +
  facet_wrap(~Metadata_Plate_Name)

bad_wells <- filter(cell_number, Proportion_solidity_filter < 0.2) %>%
  select(Metadata_Plate_Name, Well,Picture) 

#Add mean and sd of randomised control distribution to area plot
median_dist <- draw_control_dist(all_plates, n = 1000, size = 0.25, column = 'AreaShape_Area')

p <- ggplot(all_plates, aes(x = Metadata_Well, y = AreaShape_Area))
p + geom_boxplot() + 
  facet_wrap(~Metadata_Plate_Name) +
  geom_hline(yintercept = mean(median_dist$value)) +
  geom_hline(yintercept = mean(median_dist$value) + 2*sd(median_dist$value), color = 'red') +
  geom_hline(yintercept = mean(median_dist$value) - 2*sd(median_dist$value), color = 'red')


#Compare with image quality measures
image_data <- read_csv('../Analysed_data/adding_quality_features/cell_countsImage.csv') %>%
  separate(FileName_DNA, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
  select(Well, Picture, ImageQuality_Correlation_DNA_20:Metadata_Plate_Name) %>%
  right_join(all_plates, by = c('Well' = 'Metadata_Well','Metadata_Plate_Name') ) %>%
  filter(ImageQuality_Correlation_DNA_20 < 0.5 & ImageQuality_PowerLogLogSlope_DNA > -2.5)

median_dist <- draw_control_mean_dist(image_data, n = 1000, size = 0.25, column = 'AreaShape_Area')

p <- ggplot(image_data, aes(x = Metadata_Well, y = AreaShape_Area))
p + geom_boxplot()+ facet_wrap(~Metadata_Plate_Name) +
  geom_hline(yintercept = mean(median_dist$value)) +
  geom_hline(yintercept = mean(median_dist$value) + 2*sd(median_dist$value), color = 'red') +
  geom_hline(yintercept = mean(median_dist$value) - 2*sd(median_dist$value), color = 'red')

hits <- image_data %>%
  group_by(Metadata_Plate_Name, Well) %>%
  summarise(mean_area = mean(AreaShape_Area), id = unique(`Systematic ID`)) %>%
  filter(mean_area > mean(median_dist$value) + 2*sd(median_dist$value)) %>%
  write_csv('hits.csv')

all_plates_hits <- image_data %>%
  mutate(hits = ifelse(image_data$`Systematic ID` %in% hits$id, 'Hits','No Hits'))


p <- ggplot(all_plates_hits, aes(x = Well, y = AreaShape_Area, fill = hits))
p + geom_boxplot()+ facet_wrap(~Metadata_Plate_Name) +
  geom_hline(yintercept = mean(median_dist$value)) +
  geom_hline(yintercept = mean(median_dist$value) + 2*sd(median_dist$value), color = 'red') +
  geom_hline(yintercept = mean(median_dist$value) - 2*sd(median_dist$value), color = 'red')

#Let's do a t.test with the backgroudn distribution I made! Then correct p-values using fdr
#or bonferroni

control_dist <- image_data$AreaShape_Area

t_test <- image_data %>%
  group_by(Well, Metadata_Plate_Name) %>%
  summarise(p.value = t.test(x = control_dist, y = AreaShape_Area)$p.value,
            id = unique(`Systematic ID`)) %>%
  mutate(adj.pval = p.adjust(p.value, method = 'bonferroni'))
