library(tidyverse)

input_path <- 'Z:/Pers_Amalia/Screening/live_imaging_second_round_8_bit/outlines/'
fig_path <- 'analysis_live_imaging/figures/'
csv_paths <- list.files(path = input_path,
                        pattern = 'yeastObjects',
                        full.names = TRUE,
                        recursive = TRUE)

data <- csv_paths %>%
  map(read_csv) %>%    
  reduce(rbind)

#######Make plot strain function

strain_data <- read_csv('analysis_live_imaging/annotation/strain_code_ready_for_use.csv') %>%
  mutate(Plate_name = str_extract(Plate_name, 'row_[0-9]{1}'))

mini_data <- data %>% 
  inner_join(strain_data,by = c('Metadata_Well' = 'Real_well', 'Metadata_Row' = 'Plate_name'))%>%
  filter(AreaShape_Solidity > 0.9) %>%
  mutate(Metadata_Time = as.numeric(Metadata_Time)) %>%
  group_by(Metadata_Time, Metadata_Row, Metadata_Well, `Systematic ID`) %>%
  summarise(median_area = median(AreaShape_Area), sd_area = sd(AreaShape_Area)) %>%
  write_csv('analysis_live_imaging/output_data/summary_data.csv')

ggplot(mini_data,aes(x = Metadata_Time, 
                     y = median_area, 
                     group = interaction(Metadata_Row, Metadata_Well))) +
  geom_line() +
  geom_line(data = subset(mini_data, `Systematic ID` == 'wt'), 
            color = 'red') +
  theme_bw()

ggsave('traces_all_strains.pdf',
       path = fig_path)

ggplot(subset(mini_data, `Systematic ID` == 'wt'),aes(x = Metadata_Time, 
                     y = median_area, 
                     group = interaction(Metadata_Row, Metadata_Well),
                     colour = Metadata_Row)) +
  geom_line() +
  theme_bw()

ggsave('only_wild_types.pdf',
       path = fig_path)


##########Are the sizes the same as we saw before?##############
primary_screen <- read_csv('analysis_live_imaging/annotation/summary_all_hits.csv')

annot_mini_data <- mini_data %>%
  inner_join(primary_screen, by = 'Systematic ID')

ggplot(annot_mini_data,aes(x = Metadata_Time, 
                           y = median_area, 
                           group = interaction(Metadata_Row, Metadata_Well),
                           colour = size)) +
  geom_line() +
  theme_bw() +
  facet_wrap(~size)

annot_mini_data %>% filter(size == 'large') %>% arrange(desc(median_area))

