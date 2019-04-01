###This will take the whole data and make plots (hopefully)
rm(list = ls())
library('tidyverse')
library('broom')
source('functions_CPA.R')
detach(package:plyr)
dir.create('plots')
dir.create('output_rep1')

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

ggplot(wild_type, aes(x = Well, y = AreaShape_Area)) +
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
    left_join(cell_number, by = c('Well' = 'Well','Metadata_Plate_Name','ImageNumber') ) %>%
    filter(ImageQuality_Correlation_DNA_20 < 0.5 & ImageQuality_PowerLogLogSlope_DNA > -2.5 & AreaShape_Solidity > 0.9) 
  
  ggplot(data, aes(x = Well, y = AreaShape_Area)) +
    geom_boxplot(notch = T) + ylim(0,2500) +
    geom_hline(yintercept = mean_wt, color = 'red') +
    geom_hline(yintercept = mean_wt + sd_wt, color = 'red',linetype="dashed") +
    geom_hline(yintercept = mean_wt - sd_wt, color = 'red', linetype="dashed")
  
  filename <- paste0('plots/',plate_n,'.pdf')
  ggsave(filename = filename, scale = 2)
  
  test <- data %>% 
    select(Well, Metadata_Plate_Name,AreaShape_Area,`Systematic ID`) %>%
    group_by(`Systematic ID`) %>%
    add_count(Well) %>%
    filter(n > 50) %>%
    summarise(pvalue = t.test(AreaShape_Area, y = wild_type$AreaShape_Area)$p.value)
  
  pval[[i]] <- test
  names(pval)[i] <- plate_n
  
  z_scores[[i]] <- data %>%
    add_count(`Systematic ID`) %>%
    filter(n > 50) %>%
    group_by(`Systematic ID`) %>%
    summarise(mean_area = mean(AreaShape_Area, trim = 0.1), sd_area = sd(AreaShape_Area)) %>%
    mutate(cv = sd_area/mean_area)

  
  areas[[i]] <- data %>%
    add_count(`Systematic ID`) %>%
    filter(n >50) %>%
    select(AreaShape_Area,`Systematic ID`,Metadata_Plate_Name, Well)
  
  i = i+1
  rm(data)
  
}

load('lists.011018.rda')
gene_list <- stack(out)

area_tib <- bind_rows(areas)
pval_tib <- bind_rows(pval)

wt_values <- c('wt',mean_wt,sd_wt, 0, sd_wt/mean_wt, 'wt' )

z_score_tib <- bind_rows(z_scores) %>%
  mutate(cv = sd_area/mean_area) %>%
  left_join(gene_list, by = c('Systematic ID' = 'values')) %>%
  add_row(`Systematic ID` = 'wt',
          mean_area = mean_wt,
          sd_area = sd_wt, 
          cv=sd_wt/mean_wt, 
          ind = 'wt') %>%
  write_csv('output_rep1/summary_stats_rep1.csv')

dup_ids <- z_score_tib[duplicated(z_score_tib$`Systematic ID`),]$`Systematic ID`
t <- filter(z_score_tib, `Systematic ID` %in% dup_ids)


wt_areas <- wild_type %>%
  select(AreaShape_Area, Well) 
wt_areas$Metadata_Plate_Name <- 'Wt'
wt_areas$`Systematic ID` <- 'Wt'

all_areas <- bind_rows(area_tib, wt_areas) %>%
  left_join(gene_list, by = c('Systematic ID' = 'values'))

all_areas$ind <- as.character(all_areas$ind)
all_areas[which(all_areas$`Systematic ID` == 'Wt'),]$ind <- 'wt'

write.table(all_areas, 'output_rep1/cell_areas.csv', sep = ',')

ggplot(z_score_tib, aes(x = mean_area, y = cv, color = ind)) +
  geom_point()
ggsave('figs/mean_vs_cv.pdf')


only_imp <-all_areas %>%
  filter(!is.na(ind))

ggplot(subset(only_imp, ind != 'var'), aes(x=reorder(`Systematic ID`, AreaShape_Area, mean), y = AreaShape_Area, color = ind)) +
  geom_boxplot() +
  ylim(100, 1000) +
  theme(axis.text.x  = element_text(angle=90, vjust=0.5))
ggsave('figs/hits_from_sam_avg_area.pdf', scale = 2, width = 7, heigh = 7)

