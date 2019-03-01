rm(list = ls())
library(tidyverse)

all_data <- read_csv('all_live_imaging_data.csv')

ggplot(all_data, aes(x = aspect_ratio)) +
  geom_histogram()

ggplot(all_data, aes(x = Solidity)) +
  geom_histogram()

as_filter <- all_data %>%
  filter(aspect_ratio > 1.5 & Solidity > 0.9 & Time < 10) %>%
  group_by(`Systematic ID`,Time) %>%
  summarise(mean_area = mean(Area),
            median_area = median(Area),
            sd_area = sd(Area))

ggplot(as_filter, aes(x = Time, y = Area)) +
  geom_point() + 
  geom_line(aes(group = `Systematic ID`)) + 
  geom_smooth()


ggplot(all_data, aes(x = TIME/60, y = Area)) +
  geom_hex()
