source('functions_CPA.R')
library(tidyverse)
library(plyr)

area <- read_csv('output_rep1/cell_areas.csv')
summary <- read_csv('output_rep1/summary_stats_rep1.csv')

omics_lists <- load_gene_lists()s

complete_summary <- left_join(summary, omics_lists, by = 'Systematic ID')

ggplot(complete_summary, aes(x = mean_area, y = cv, color = type)) +
  geom_point()
ggsave('figs/omics_summary.pdf')


complete_area <- left_join(area, omics_lists, by = 'Systematic ID') %>%
  filter(!is.na(type) | ind == 'wt') 

complete_area[which(complete_area$ind == 'wt'),]$type <- 'wt'

ggplot(subset(complete_area,type == 'wt' | type == 'Up RNA'), 
       aes(x=reorder(`Systematic ID`, AreaShape_Area, mean), 
           y = AreaShape_Area, 
           color = type)) +
  geom_boxplot() +
  ylim(100, 1000) 


ggplot(subset(complete_area,type == 'wt' | type == 'Down RNA'), 
       aes(x=reorder(`Systematic ID`, AreaShape_Area, mean), 
           y = AreaShape_Area, 
           color = type)) +
  geom_boxplot() +
  ylim(100, 1000) 

ggplot(subset(complete_area,type == 'wt' | type == 'Up prot'), 
       aes(x=reorder(`Systematic ID`, AreaShape_Area, mean), 
           y = AreaShape_Area, 
           color = type)) +
  geom_boxplot() +
  ylim(100, 1000) 

ggplot(subset(complete_area,type == 'wt' | type == 'Down prot'), 
       aes(x=reorder(`Systematic ID`, AreaShape_Area, mean), 
           y = AreaShape_Area, 
           color = type)) +
  geom_boxplot() +
  ylim(100, 1000) 



ggsave('figs/omics_avg_area.pdf')


