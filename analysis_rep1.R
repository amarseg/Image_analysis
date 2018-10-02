###This will take the whole data and make plots (hopefully)
rm(list = ls())
library('tidyverse')
library('broom')
source('functions_CPA.R')
dir.create('plots')

wild_type_path <- '../Analysed_data/wt_data/Nuclei_AR_Solidity_Filtered.csv'

wild_type <- read_csv(wild_type_path) %>%
  filter(Metadata_QCFlag == 0) %>%
  separate(Metadata_FileLocation, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
  transform(Well_n = str_sub(Well_n, start = 5, end = 7))

cell_number <- read_csv('../Analysed_data/aug_sept_2018/cell_countsImage.csv') %>%
  separate(FileName_DNA, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
  mutate(Proportion_AR_filter = Count_Nuclei_AR_filtered/(Count_Nuclei+Count_Nuclei_AR_filtered), 
         Proportion_solidity_filter = Count_Nuclei_AR_Solidity_Filtered/(Count_Nuclei_AR_filtered+ Count_Nuclei_AR_Solidity_Filtered))

ggplot(wild_type, aes(x = Well_n, y = AreaShape_Area)) +
  geom_boxplot()

mean_wt <- mean(wild_type$AreaShape_Area, trim = 0.1)
sd_wt <- sd(wild_type$AreaShape_Area)

test <- gather(cell_number, key = feature, value = cell_count, -Metadata_Plate_Name, -Well) %>%
  filter(feature == 'Proportion_AR_filter' | feature == 'Proportion_solidity_filter')

all_files <- list.files('../Analysed_data/aug_sept_2018/', full.names = T, recursive = T, pattern = 'Nuclei_AR_Solidity_Filtered.csv')
all_files <- all_files[1:36]

pval <- list()
z_scores <- list()
areas <- list()
i = 1
for(file in all_files)
{
  plate_n <- str_extract(file, pattern = 'Plate[0-9]{1,2}')
  data <- load_filtered_data(file, plate_n = plate_n) %>%
    left_join(cell_number, by = c('Well' = 'Well','Metadata_Plate_Name') ) %>%
    filter(ImageQuality_Correlation_DNA_20 < 0.5 & ImageQuality_PowerLogLogSlope_DNA > -2.5 & AreaShape_Solidity > 0.9) 
  
  ggplot(data, aes(x = Well, y = AreaShape_Area)) +
    geom_boxplot(notch = T) + ylim(0,2500) +
    geom_hline(yintercept = mean_wt, color = 'red') +
    geom_hline(yintercept = mean_wt + sd_wt, color = 'red',linetype="dashed") +
    geom_hline(yintercept = mean_wt - sd_wt, color = 'red', linetype="dashed")
  
  filename <- paste0('plots/',plate_n,'.pdf')
  ggsave(filename = filename, scale = 2)
  
  test <- data %>% 
    select(Well, Metadata_Plate_Name,AreaShape_Area) %>%
    group_by(Well,Metadata_Plate_Name) %>%
    add_count(Well) %>%
    filter(n > 50) %>%
    summarise(pvalue = t.test(AreaShape_Area, y = wild_type$AreaShape_Area)$p.value, Plate = plate_n)
  
  pval[[i]] <- test
  names(pval)[i] <- plate_n
  
  z_scores[[i]] <- data %>%
    group_by(Well, Metadata_Plate_Name) %>%
    add_count(Well) %>%
    filter(n > 50) %>%
    summarise(mean_area = mean(AreaShape_Area, trim = 0.1), sd_area = sd(AreaShape_Area)) %>%
    mutate(z_score = (mean_area - mean_wt)/sd_wt)

  
  areas[[i]] <- data %>%
    add_count(Well) %>%
    filter(n >50) %>%
    select(AreaShape_Area,`Systematic ID`,Metadata_Plate_Name)
  
  i = i+1
  rm(data)
  
}

pval_df <-as.data.frame(pval)


