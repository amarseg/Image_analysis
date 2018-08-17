#Comparison between the two replicates
#one after 12 hours of analogue, and the other one after 11
library(tidyverse)
source('functions_CPA.R')
dir.create('figs/replicate_comparison')

##
#let's load image quality data to be able to remove shitty images
##

qc_rep1 <- read_csv('../Analysed_data/adding_quality_features/cell_countsImage.csv') %>%
  separate(FileName_DNA, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
  select(Well, Picture, ImageQuality_Correlation_DNA_20:Metadata_Plate_Name) %>%
  gather(key = feature, value = feature_value, -Well, -Picture, -Metadata_Plate_Name)

qc_rep2 <- read_csv('../Analysed_data/11hr_exp/cell_countsImage.csv') %>%
  separate(FileName_DNA, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
  select(Well, Picture, ImageQuality_Correlation_DNA_20:Metadata_Plate_Name) %>%
  gather(key = feature, value = feature_value, -Well, -Picture, -Metadata_Plate_Name)


ggplot(qc_rep1, aes(x = feature_value)) +
  geom_histogram() +
  facet_wrap(~feature, scales = 'free') %>% ggsave('figs/replicate_comparison/features_rep1_pdf', scale = 2)

ggplot(qc_rep2, aes(x = feature_value)) +
  geom_histogram() +
  facet_wrap(~feature, scales = 'free') %>% ggsave('figs/replicate_comparison/features_rep2_pdf', scale = 2)

qc_data <- bind_rows(qc_rep1, qc_rep2, .id = 'id' )
ggplot(qc_data, aes(x = feature_value, fill = id)) +
  geom_histogram(alpha = 0.8) +
  facet_wrap(~feature, scales = 'free') %>% ggsave('figs/replicate_comparison/comparison_image_quality.pdf')

##
#Load image data
##
plate1 <- load_filtered_data('../Analysed_data/adding_quality_features/Plate1_20f/Nuclei_AR_Solidity_Filtered.csv', plate_n = 1)
plate2 <- load_filtered_data('../Analysed_data/adding_quality_features/Plate2_20f/Nuclei_AR_Solidity_Filtered.csv', plate_n = 2)
plate3 <- load_filtered_data('../Analysed_data/adding_quality_features/Plate3_20f/Nuclei_AR_Solidity_Filtered.csv', plate_n = 3)
plate4 <- load_filtered_data('../Analysed_data/adding_quality_features/Plate4_wPI/Nuclei_AR_Solidity_Filtered.csv', plate_n = 4)

rep1 <- bind_rows(plate1, plate2, plate3, plate4) %>%
  group_by(Metadata_Plate_Name, Metadata_Well) %>%
  summarise(avg_area = mean(AreaShape_Area)) %>%
  separate(Metadata_Plate_Name, into = c('Metadata_Plate_Name','bad_bit'), sep = '_') %>%
  select(-bad_bit)

plate1 <- load_filtered_data('../Analysed_data/11hr_exp/Plate1/Nuclei_AR_Solidity_Filtered.csv', plate_n = 1)
plate2 <- load_filtered_data('../Analysed_data/11hr_exp/Plate2/Nuclei_AR_Solidity_Filtered.csv', plate_n = 2)
plate3 <- load_filtered_data('../Analysed_data/11hr_exp/Plate3/Nuclei_AR_Solidity_Filtered.csv', plate_n = 3)
plate4 <- load_filtered_data('../Analysed_data/11hr_exp/Plate4/Nuclei_AR_Solidity_Filtered.csv', plate_n = 4)

rep2 <- bind_rows(plate1, plate2, plate3, plate4) %>% 
  group_by(Metadata_Plate_Name, Metadata_Well) %>%
  summarise(avg_area = mean(AreaShape_Area))

all_reps <- inner_join(rep1, rep2, suffix = c('_rep1','_rep2'), by = c('Metadata_Plate_Name', 'Metadata_Well')) %>%
  mutate(ratio = avg_area_rep1/avg_area_rep2)

ggplot(all_reps, aes(x = avg_area_rep1, y = avg_area_rep2, label = Metadata_Well)) +
  geom_point() +
  facet_wrap(~Metadata_Plate_Name, scales = 'free') +
  geom_label() %>% ggsave('figs/replicate_comparison/mean_area_plot_comparison.pdf')

ggplot(all_reps, aes(x = ratio)) +
  geom_histogram() +
  facet_wrap(~Metadata_Plate_Name, scales = 'free') +
  scale_x_log10() + 
  geom_rug() %>% ggsave('figs/replicate_comparison/ratio_between_replicates_histogram.pdf')




