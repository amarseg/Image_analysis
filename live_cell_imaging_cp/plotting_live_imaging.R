library(tidyverse)

all_data <- read_csv('live_cell_imaging_cp/output_data_cleaned.csv') %>%
  filter(`Systematic ID` == 'wt')

ggplot(all_data, aes(x = Time, y = AreaShape_Area, group = `Systematic ID`)) +
  geom_point() +
  stat_summary(fun.y = 'mean', color = 'red', geom = 'line') +
  facet_grid(~Well)


