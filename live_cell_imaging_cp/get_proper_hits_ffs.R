library(tidyverse)

rep1 <- read_csv('output_rep1/statistics_rep1_good.csv') %>%
  mutate(all_area = mean(mean_area), sd_all = sd(mean_area)) %>%
  mutate(z_score = (mean_area - all_area)/sd_all) %>%
  mutate(pvalue = 1-pnorm(abs(z_score))) %>% 
  mutate(hits = ifelse(pvalue < 0.15, 'hit','no hit')) 

rep2 <- read_csv('output_rep2/statistics_rep2.csv') %>%
  mutate(all_area = mean(mean_area), sd_all = sd(mean_area)) %>%
  mutate(z_score = (mean_area - all_area)/sd_all) %>%
  mutate(pvalue = 1-pnorm(abs(z_score))) %>% 
  mutate(hits = ifelse(pvalue < 0.15, 'hit','no hit')) 


all_things <- inner_join(rep1, rep2, by   = 'Systematic ID', suffix = c('_rep1','_rep2'))

strain_data <- read_csv('library_strains.csv') %>%
  dplyr::select('Ver5.0 position','Systematic ID') %>%
  separate('Ver5.0 position', into = c('Version','Plate','Well')) 

hits_i_guess <- all_things %>%
  filter(hits_rep1 == 'hit' & hits_rep2 == 'hit') %>%
  inner_join(strain_data, by = 'Systematic ID')
