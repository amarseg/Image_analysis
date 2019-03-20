library(tidyverse)

second_list <- read_csv('second_list_of_hits.csv')

ggplot(second_list, aes(x = mean_area_rep1, y = mean_area_rep2)) +
  geom_point()

second_list$size = NA
second_list[which(second_list$mean_area_rep1 < 500 & second_list$mean_area_rep2 < 600),]$size <- 'small'
second_list[which(second_list$mean_area_rep1 > 500 & second_list$mean_area_rep2 > 600),]$size <- 'large'

second_list_clean <- second_list %>%
  select(`Systematic ID`, size)

######create master list##########
#load 1rst list#
first_list <- read_csv('z_score_hits.csv') %>%
  select(`Systematic ID`, size)

all_list <- bind_rows(second_list_clean, first_list) %>%
  write_csv('final_hit_list.csv')
