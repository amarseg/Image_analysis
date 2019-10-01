library(tidyverse)
source('functions_CPA.R')

rep1 <- c('Z:/Pers_Amalia/Screening/Analysed_data/aug_sept_2018/')
rep2 <- c('Z:/Pers_Amalia/Screening/Analysed_data/replicate_2/', 
          'Z:/Pers_Amalia/Screening/Analysed_data/replicate_2_batch_2/',
          'Z:/Pers_Amalia/Screening/Analysed_data/batch_3/',
          'Z:/Pers_Amalia/Screening/Analysed_data/batch_4/')

plates <- c(1:36)
plates <- paste0('Plate',plates)
plates <- paste(plates, collapse = '|')

areas1 <- list()
i = 1
for(file in rep1_extract)
{
  plate_n <- str_extract(file, pattern = 'Plate[0-9]{1,2}')
  data <- load_filtered_data(file, plate_n = plate_n) %>%
    left_join(cell_number, by = c('Well' = 'Well','Metadata_Plate_Name','ImageNumber') ) %>%
    filter(ImageQuality_Correlation_DNA_20 < 0.5 & ImageQuality_PowerLogLogSlope_DNA < -2 & AreaShape_Solidity > 0.9) 
  
  areas1[[i]] <- data %>%
    add_count(`Systematic ID`) %>%
    filter(n >=50) %>%
    select(AreaShape_Area,`Systematic ID`,Metadata_Plate_Name, Well)
  i=i+1
}

areas2 <- list()
i = 1
for(file in rep2_folders)
{
  plate_n <- str_extract(file, pattern = 'Plate[0-9]{1,2}')
  data <- load_filtered_data(file, plate_n = plate_n) %>%
    left_join(cell_number_2, by = c('Well' = 'Well','Metadata_Plate_Name','ImageNumber') ) %>%
    filter(ImageQuality_Correlation_DNA_20 < 0.5 & ImageQuality_PowerLogLogSlope_DNA < -2 & AreaShape_Solidity > 0.9) 
  
  areas2[[i]] <- data %>%
    add_count(`Systematic ID`) %>%
    filter(n >=50) %>%
    select(AreaShape_Area,`Systematic ID`,Metadata_Plate_Name, Well)
  
  
  i=i+1
}

areas1_tib <- bind_rows(areas1) %>%
  write_csv('clean_code/areas_rep1.csv') %>%
  add_column(rep = 'replicate1')

areas2_tib <- bind_rows(areas2) %>%
  write_csv('clean_code/areas_rep2.csv') %>%
  add_column(rep = 'replicate2')


summary_stats1 <- areas1_tib %>%
  group_by(`Systematic ID`) %>%
  summarise(mean_area = mean(AreaShape_Area,trim = 0.1), sd_area = sd(AreaShape_Area)) %>%
  add_column(all_mean = mean(areas1_tib$AreaShape_Area, trim = 0.1), all_sd = sd(areas1_tib$AreaShape_Area)) %>%
  mutate(z_score = (mean_area - all_mean)/all_sd) %>%
  mutate(pval = pnorm(abs(z_score)))
  

summary_stats2 <- areas2_tib %>%
  group_by(`Systematic ID`) %>%
  summarise(mean_area = mean(AreaShape_Area, trim = 0.1), sd_area = sd(AreaShape_Area)) %>%
  add_column(all_mean = mean(areas2_tib$AreaShape_Area, trim = 0.1), all_sd = sd(areas2_tib$AreaShape_Area)) %>%
  mutate(z_score = (mean_area - all_mean)/all_sd) %>%
  mutate(pval = pnorm(abs(z_score)))

wild_type_path <- '../Analysed_data/wt_data/Nuclei_AR_Solidity_Filtered.csv'

wild_type <- read_csv(wild_type_path) %>%
  filter(Metadata_QCFlag == 0) %>%
  separate(Metadata_FileLocation, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
  transform(Well_n = str_sub(Well_n, start = 5, end = 7)) %>%
  transform(Well = str_extract(Well,pattern = '[A-Z][0-9]{1,2}'))

cell_number <- read_csv('../Analysed_data/aug_sept_2018/cell_countsImage.csv') %>%
  separate(FileName_DNA, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
  mutate(Proportion_AR_filter = Count_Nuclei_AR_filtered/(Count_Nuclei+Count_Nuclei_AR_filtered), 
         Proportion_solidity_filter = Count_Nuclei_AR_Solidity_Filtered/(Count_Nuclei_AR_filtered+ Count_Nuclei_AR_Solidity_Filtered)) 

wild_type <- left_join(wild_type, cell_number, by = 'Well') %>%
  filter(ImageQuality_Correlation_DNA_20< 0.5 & ImageQuality_PowerLogLogSlope_DNA > -2.5 & AreaShape_Solidity > 0.9) 

id <- 'SPAC6C3.08'

interesting <- areas1_tib %>%
  filter(`Systematic ID` == id)

ggplot(wild_type, aes(AreaShape_Area)) +
  facet_wrap(~Metadata_Well) +
  geom_histogram() +
  geom_histogram(data = interesting, aes(AreaShape_Area), colour = 'red') +
  geom_vline(xintercept=1000, colour = 'red') +
  coord_cartesian(ylim=c(0, 1000))



all_areas <- bind_rows(areas1_tib, areas2_tib)
for(id in unique(all_areas$`Systematic ID`))
{
  todo <- filter(all_areas, `Systematic ID` == id)
  
  ggplot(todo, aes(AreaShape_Area, fill = rep)) +
    geom_histogram(alpha = 0.75) +
    ggtitle(id)
  
  ggsave(paste0('clean_code/',id,'.jpg'))
}
