library(tidyverse)
library(broom)
library('reticulate')

all_data <- read_csv('live_cell_imaging_cp/output_data_cleaned.csv') %>%
  rename('ID' = `Systematic ID`) %>%
  filter(Time < 4) %>%
  group_by(ID, Time) %>%
  summarise(avg_area = median(AreaShape_Area))

area_lms <- all_data %>%
  group_by(ID) %>%
  do(area_lms = lm(avg_area ~ Time, data = .))

slopes <- tidy(area_lms, area_lms)

coeffs <- glance(area_lms, area_lms) %>%
  add_column(p_adj = p.adjust(.$p.value, method = 'BY'))

ggplot(coeffs, aes(p_adj)) +
  geom_histogram()

ggplot(coeffs, aes(r.squared)) +
  geom_histogram()

ggplot(filter(slopes, term == 'Time'), aes(estimate)) +
  geom_histogram()

dir.create('live_cell_imaging_cp/figures/regression/')

for (id in unique(all_data$ID)) 
{
  sub_id <- filter(all_data, ID == id)
  ggplot(sub_id, aes(x = Time, y = avg_area)) +
    geom_point()+
    geom_line() +
    geom_smooth(method = 'lm')
  
  file_name = paste0('live_cell_imaging_cp/figures/regression/',id,'.jpg')
  ggsave(file_name)
}

############Is the screening confirmed? ######

annot_data <- read_csv('final_hit_list.csv')

annot_slopes <- slopes %>%
  left_join(annot_data, by = c('ID'= 'Systematic ID'))

annot_slopes[annot_slopes$ID =='wt',]$size <- 'wt'


ggplot(filter(annot_slopes, term == 'Time'), aes(estimate, fill = size)) +
  geom_histogram(position = 'dodge')

annot_slopes$slope_res <- NA 
annot_slopes[which(annot_slopes$estimate > 25 & annot_slopes$term == 'Time'),]$slope_res <- 'positive'
annot_slopes[which(annot_slopes$estimate < 25 & annot_slopes$term == 'Time'),]$slope_res <- 'negative'


output_res <- annot_slopes %>%
  filter(term == 'Time') %>%
  select(ID, size, slope_res) %>%
  write_csv('live_cell_imaging_cp/movie_results.csv')

############Do something similar with the tracked data, app. not good##########

track_data <- read_csv('live_cell_imaging_cp/track_cell_data.csv') %>%
  rename('ID' = `Systematic ID`)

track_lms <- track_data %>%
  group_by(ID) %>%
  do(track_lms = lm(log2(AreaShape_Area) ~ TIME, data = .))

slopes_track <- tidy(track_lms, track_lms)

coeffs_track <- glance(track_lms, track_lms) %>%
  add_column(p_adj = p.adjust(.$p.value, method = 'BY'))

ggplot(coeffs_track, aes(p_adj)) +
  geom_histogram()

ggplot(coeffs_track, aes(r.squared)) +
  geom_histogram()

ggplot(filter(slopes_track, term == 'TIME'), aes(estimate)) +
  geom_histogram()

dir.create('live_cell_imaging_cp/figures/track_regression/')

for (id in unique(track_data$ID)) 
{
  sub_id <- filter(track_data, ID == id)
  ggplot(sub_id, aes(x = Time, y = AreaShape_Area)) +
    geom_point() +
    geom_smooth(method = 'lm')
  
  file_name = paste0('live_cell_imaging_cp/figures/track_regression/',id,'.jpg')
  ggsave(file_name)
}

