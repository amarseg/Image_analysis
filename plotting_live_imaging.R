rm(list = ls())
library(tidyverse)
library('nlme')
library('broom')
all_data <- read_csv('all_live_imaging_data.csv') %>%
  filter(Area > 500)

ggplot(all_data, aes(x = aspect_ratio)) +
  geom_histogram()

ggplot(all_data, aes(x = Solidity)) +
  geom_histogram()

as_filter <- all_data %>%
  #filter(aspect_ratio > 1.5 & Solidity > 0.9) %>%
  group_by(`Systematic ID`,Time, Well) %>%
  summarise(mean_area = mean(Area),
            median_area = median(Area),
            sd_area = sd(Area), 
            TIME = mean(TIME), 
            plate = unique(plate)) 

test_gene <- all_data %>%
  filter((`Systematic ID` == 'SPCC285.17' | `Systematic ID` == 'SPCC18B5.09c') & Time <12)


ggplot(test_gene, aes(x = Time, y = log2(Area), group = interaction(Time, `Systematic ID`), color = `Systematic ID`)) +
  geom_boxplot() 


ggplot(all_data, aes(x = TIME/60, y = Area)) +
  geom_hex()

only_wt <- as_filter %>%
  filter( `Systematic ID` == 'wt')

ggplot(subset(only_wt, Time < 5) , aes(x = TIME/60, y = median_area, color = plate)) +
  geom_point() + 
  geom_line(aes(group = Well)) 


###################Add data from previous screening##################
hits <- read_csv('z_score_hits.csv') %>%
  select(`Systematic ID`,hits, size)

all_data_annotated <- left_join(all_data, hits, by = 'Systematic ID') %>%
  write_csv('all_data_live_imaging_annotated.csv')

all_data_annotated[all_data_annotated$`Systematic ID`=='wt',]$size <- 'wt'

ggplot(subset(all_data_annotated, Time < 12), aes(x = Time, y = log2(Area))) +
  geom_boxplot(aes(group = Time))+
  facet_wrap(~size)


#############Fit line to data and create plots for everybody###############
dir.create('output_live_imaging')

by_id <- as_filter %>%
  group_by(`Systematic ID`)

res <-do(by_id, glance(lm(median_area ~ TIME, data = .))) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(.$p.value, method = 'BY'))

slopes <-do(by_id, tidy(lm(median_area ~ TIME, data = .))) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(.$p.value, method = 'BY')) %>%
  filter(term == 'TIME') 


ggplot(res, aes(p.value)) +
  geom_histogram()

ggplot(res, aes(p_adj)) +
  geom_histogram()

ggplot(res, aes(r.squared)) +
  geom_histogram()

p_val_cutoff <- 0.05
r_squared_cutoff <- 0.5

good_hits <- slopes %>%
  filter(p_adj < p_val_cutoff & estimate > r_squared_cutoff)

########Plots all mutants#####################
all_data_annotated$color = NA

all_data_annotated[which(all_data_annotated$size == 'small'),]$color <-'blue' 
all_data_annotated[which(all_data_annotated$size == 'large'),]$color <-'red' 
all_data_annotated[which(all_data_annotated$size == 'wt'),]$color <- 'brown' 

all_data_annotated <- all_data_annotated %>%
  filter(Time < 12)

for (id in unique(as_filter$`Systematic ID`))
{
  data <- subset(as_filter, `Systematic ID` == id)
  
  file_name = paste0('output_live_imaging/',id,'.jpg')
  
  boxplot(log2(Area) ~ Time, data, col = unique(data$color), outline = F)
  dev.copy(jpeg,filename = file_name)
  dev.off()
}

