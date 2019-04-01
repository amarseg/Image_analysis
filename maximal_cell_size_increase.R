library(tidyverse)

all_data <- read_csv('live_cell_imaging_cp/output_data_cleaned.csv') %>%
  rename('ID' = `Systematic ID`) %>%
  #filter(Time < 4) %>%
  group_by(ID, Time) %>%
  summarise(avg_area = median(AreaShape_Area))

annot_data <- read_csv('final_hit_list.csv')

tp0 <- filter(all_data, Time == 0)

max_size <- all_data %>%
  group_by(ID) %>%
  summarise(max_area = max(avg_area)) %>%
  left_join(annot_data, by = c('ID' = 'Systematic ID')) %>%
  left_join(tp0, by = 'ID', suffix = c('_max','_tp0')) %>%
  mutate(max_ratio = max_area/avg_area)

max_size[max_size$ID =='wt',]$size <- 'wt'


ggplot(max_size, aes(max_area, fill = size)) +
  geom_histogram(position = 'dodge')

ggplot(max_size, aes(max_ratio, fill = size)) +
  geom_histogram(position = 'dodge')
