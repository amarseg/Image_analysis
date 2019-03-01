library('tidyverse')

cell_number <- read_csv('../Analysed_data/aug_sept_2018/cell_countsImage.csv') %>%
  separate(FileName_DNA, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
  mutate(Proportion_AR_filter = Count_Nuclei_AR_filtered/(Count_Nuclei+Count_Nuclei_AR_filtered), 
         Proportion_solidity_filter = Count_Nuclei_AR_Solidity_Filtered/(Count_Nuclei_AR_filtered+ Count_Nuclei_AR_Solidity_Filtered)) 

cell_count_files <- list.files(rep2, full.names = T, recursive = T, pattern = 'cell_countsImage.csv') %>%
  map(read_csv) %>%
  purrr::reduce(bind_rows)

cell_number_2 <- cell_count_files %>%
  separate(FileName_DNA, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--') %>%
  mutate(Proportion_AR_filter = Count_Nuclei_AR_filtered/(Count_Nuclei+Count_Nuclei_AR_filtered), 
         Proportion_solidity_filter = Count_Nuclei_AR_Solidity_Filtered/(Count_Nuclei_AR_filtered+ Count_Nuclei_AR_Solidity_Filtered))

areas_1 <- read_csv('output_rep1/cell_areas.csv')
areas_2 <- read_csv('output_rep2/areas_rep2.csv')

t <- sum(cell_number_2$Count_Nuclei) + sum(cell_number$Count_Nuclei)
