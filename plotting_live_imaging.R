rm(list = ls())
library(tidyverse)

all_data <- read_csv('all_live_imaging_data.csv')

ggplot(all_data, aes(x = aspect_ratio)) +
  geom_histogram()

ggplot(all_data, aes(x = Solidity)) +
  geom_histogram()

as_filter <- all_data %>%
  filter(aspect_ratio > 1.5 & Solidity > 0.9) %>%
  group_by(`Systematic ID`,Time) %>%
  summarise(mean_area = mean(Area),
            median_area = median(Area),
            sd_area = sd(Area))

test_gene <- all_data %>%
  filter((`Systematic ID` == 'SPCC285.17' | `Systematic ID` == 'SPCC18B5.09c') & Time <12)


ggplot(test_gene, aes(x = Time, y = log2(Area), group = interaction(Time, `Systematic ID`), color = `Systematic ID`)) +
  geom_boxplot() 


ggplot(all_data, aes(x = TIME/60, y = Area)) +
  geom_hex()

only_wt <- as_filter %>%
  filter( `Systematic` == 'wt')

ggplot(only_wt, aes(x = Time, y = median_area)) +
  geom_point() + 
  geom_line() + 
  geom_smooth()


###################Add data from previous screening##################
hits <- read_csv('z_score_hits.csv') %>%
  select(`Systematic ID`,hits, size)

all_data_annotated <- left_join(all_data, hits, by = 'Systematic ID') %>%
  write_csv('all_data_live_imaging_annotated.csv')

all_data_annotated[all_data_annotated$`Systematic ID`=='wt',]$size <- 'wt'

ggplot(subset(all_data_annotated, Time < 12), aes(x = Time, y = log2(Area))) +
  geom_boxplot(aes(group = Time))+
  facet_wrap(~size)
