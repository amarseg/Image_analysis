library(tidyverse)

wt_data <- read_csv('live_cell_imaging_cp/output_data_cleaned.csv') %>%
  filter(ImageQuality_PowerLogLogSlope_DNA > -1 & ImageQuality_Correlation_DNA_20 > 0.25) %>%
  filter(`Systematic ID` == 'wt' & Time < 4 )

ggplot(wt_data, aes(x = Time, y = AreaShape_Area, group = `Systematic ID`)) +
  geom_point() +
  stat_summary(fun.y = 'mean', color = 'red', geom = 'line') +
  facet_grid(~Well)


dir.create('live_cell_imaging_cp/figures')


all_data <- read_csv('live_cell_imaging_cp/output_data_cleaned.csv') %>%
  filter(ImageQuality_PowerLogLogSlope_DNA > -1 & ImageQuality_Correlation_DNA_20 > 0.25) %>%
  group_by(Time, Position, `Systematic ID`) %>%
  mutate(avg_time = mean(TIME)) %>%
  ungroup()

hits <- read_csv('z_score_hits.csv') %>%
  select(`Systematic ID`,hits, size)

all_data_annotated <- left_join(all_data, hits, by = 'Systematic ID') %>%
  write_csv('live_cell_imaging_cp/all_data_live_imaging_annotated.csv')

all_data_annotated[all_data_annotated$`Systematic ID`=='wt',]$size <- 'wt'

all_data_annotated$color = NA

all_data_annotated[which(all_data_annotated$size == 'small'),]$color <-'blue' 
all_data_annotated[which(all_data_annotated$size == 'large'),]$color <-'red' 
all_data_annotated[which(all_data_annotated$size == 'wt'),]$color <- 'brown' 


for(id in unique(all_data_annotated$`Systematic ID`))
{
  sub_id <- all_data_annotated %>%
    filter(`Systematic ID` == id)
  
  file_name = paste0('live_cell_imaging_cp/figures/',id, '.jpg')
  
  boxplot(log2(AreaShape_Area) ~ Time, data = sub_id, col = unique(sub_id$color))
  dev.copy(jpeg,filename = file_name)
  dev.off()
  
  ggplot(sub_id, aes(x = Time, y = log2(AreaShape_Area), group = Time)) +
    geom_violin() +
    stat_summary(fun.y = 'median', geom = 'point') +
    stat_summary(fun.y = 'median', geom = 'line')
    
  file_name = paste0('live_cell_imaging_cp/figures/',id,'_violin.jpg')
  ggsave(file_name)
}

##############Plot track data#################

track_data <- read_csv('live_cell_imaging_cp/track_cell_data.csv') %>%
  filter(cell_n != 0) %>%
  left_join( hits, by = 'Systematic ID') 

track_data[track_data$`Systematic ID`=='wt',]$size <- 'wt'


ggplot(track_data, aes(x = TIME, y = AreaShape_Area, group = interaction(Position, cell_n))) +
  geom_line()+
  facet_wrap(~size) +
  geom_smooth()

for (id in unique(track_data$`Systematic ID`))
{
  sub_id <- filter(track_data, `Systematic ID` == id)
  
  ggplot(track_data, aes(x = TIME, y = AreaShape_Area, group = interaction(Position, cell_n))) +
    geom_line()+
    stat_summary(fun.y = 'median', geom = 'line', color = 'red')
    geom_smooth()
    
  file_name = paste0('live_cell_imaging_cp/figures/',id,'_track.jpg')
  ggsave(file_name)
}
#############Plot summary of different strains################
sum_track <- track_data %>%
  group_by(`Systematic ID`, Time) %>%
  summarise(avg_area = median(AreaShape_Area), avg_time = mean(TIME), size = unique(size))

ggplot(sum_track, aes(x = avg_time, y = log2(avg_area), group = `Systematic ID`)) +
  geom_line() +
  facet_grid(~size) 
