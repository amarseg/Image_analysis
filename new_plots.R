library(tidyverse)

wild_type_path <- '../Analysed_data/wt_data/Nuclei_AR_Solidity_Filtered.csv'

wild_type <- read_csv(wild_type_path) %>%
  filter(Metadata_QCFlag == 0) %>%
  separate(Metadata_FileLocation, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
  transform(Well_n = str_sub(Well_n, start = 5, end = 7)) %>%
  transform(Well = str_extract(Well,pattern = '[A-Z][0-9]{1,2}')) 

cell_number <- read_csv('../Analysed_data/wt_data/cell_countsImage.csv') %>%
  separate(FileName_DNA, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
  mutate(Proportion_AR_filter = Count_Nuclei_AR_filtered/(Count_Nuclei+Count_Nuclei_AR_filtered), 
         Proportion_solidity_filter = Count_Nuclei_AR_Solidity_Filtered/(Count_Nuclei_AR_filtered+ Count_Nuclei_AR_Solidity_Filtered)) 

wild_type <- left_join(wild_type, cell_number, by = c('Well','Picture')) %>%
  filter(ImageQuality_Correlation_DNA_20< 0.5 & ImageQuality_PowerLogLogSlope_DNA > -2.5 & AreaShape_Solidity > 0.9) %>%
  select(AreaShape_Area) 

wild_type$Metadata_Plate_Name <- 'Wt'
wild_type$`Systematic ID` <- 'Wt'
wild_type$rep <- 'Wt'


area_1 <- read_csv('clean_code/areas_rep1.csv') %>%
  add_column(rep = 'replicate1')

area_2 <- read_csv('clean_code/areas_rep2.csv') %>%
  add_column(rep = 'replicate2')


all_areas <- bind_rows(area_1, area_2, wild_type)

hits <- read_csv('live_cell_imaging_cp/hits_i_guess.csv')

only_hits <- all_areas %>%
  filter(`Systematic ID` %in% hits$`Systematic ID` | str_detect(`Systematic ID`,'Wt'))


wild_type_median <- all_areas %>%
  filter(str_detect(`Systematic ID`, 'Wt')) %>%
  summarise(median_area = median(AreaShape_Area))

library(cowplot)

gene_ids <- read_tsv('gene_IDs_names.tsv', skip = 1, col_names = c('Systematic ID','name','extra_names'))

gene_ids[is.na(gene_ids$name),]$name <- gene_ids[is.na(gene_ids$name),]$`Systematic ID`


annot_hits <- only_hits %>%
  left_join(gene_ids, by = c('Systematic ID'))

annot_hits[is.na(annot_hits$name),]$name <- 'wt'

ggplot(annot_hits, aes(x = reorder(name ,AreaShape_Area, median), y = AreaShape_Area)) +
  geom_boxplot() +
  geom_boxplot(data = filter(annot_hits,str_detect(`Systematic ID`, 'Wt')), colour = 'red') +
  coord_cartesian(xlim = c(0,2000)) +
  theme_cowplot(12) +
  coord_flip(ylim = c(0,2000)) +
  ylab('Area (px^2)') +
  xlab('Systematic ID') +
  geom_hline(yintercept = 544, colour = 'red', linetype = 'dashed', lwd = 1.5)
