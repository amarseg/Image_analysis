library(tidyverse)
library(pheatmap)
#Load data#
all_data <- read_csv('live_cell_imaging_cp/output_data_cleaned.csv') %>%
  rename('ID' = `Systematic ID`) %>%
  #filter(Time < 4) %>%
  group_by(ID, Time) %>%
  summarise(avg_area = median(AreaShape_Area))

annot_data <- read_csv('final_hit_list.csv')

tp0 <- filter(all_data, Time == 0)

wt_max <- all_data %>%
  filter(ID == 'wt')

wt_max <- max(wt_max$avg_area)

max_size <- all_data %>%
  group_by(ID) %>%
  summarise(max_area = max(avg_area)) %>%
  left_join(annot_data, by = c('ID' = 'Systematic ID')) %>%
  left_join(tp0, by = 'ID', suffix = c('_max','_tp0')) %>%
  mutate(max_ratio = max_area/avg_area) %>%
  mutate(wt_ratio = max_area/wt_max)


hit_list <- read_csv('output_rep1/summary_rep1_pval_0.15.csv') %>%
  filter(hits == 'hit')
second_list <- read_csv('output_rep2/statistics_rep2_pvalue_0.15.csv') %>%
  filter(hits == 'hit')

all_hits <- inner_join(hit_list, second_list, by = 'Systematic ID', suffix = c('_rep1','_rep2')) %>%
  select(`Systematic ID`,z_score_rep1, z_score_rep2)  %>%
  left_join(max_size, by = c('Systematic ID' = 'ID')) %>%
  write_csv('live_cell_imaging_cp/summary_all_hits.csv')

int <- select(all_hits, z_score_rep1, z_score_rep2, max_ratio, wt_ratio, `Systematic ID`)
int$max_ratio <- log2(int$max_ratio)
int$wt_ratio <- log2(int$wt_ratio)

col <- colorRampPalette(c('blue','gray','red'))
cl <- pheatmap(int[,-5], cluster_cols = F, labels_row = int$`Systematic ID`, cutree_rows = 10, breaks= seq(-3.5,3.5,0.5),
               color = col(14))

hcl <- data.frame(cluster = cutree(cl$tree_row, k = 10), id = int[,5]) %>%
  write_csv('live_cell_imaging_cp/result_clustering.csv')

########Lines lines lines############

all_hits_complete <- inner_join(hit_list, second_list, by = 'Systematic ID', suffix = c('_rep1','_rep2'))
all_hits_complete$size = NA
all_hits_complete[which(all_hits_complete$mean_area_rep1 < 500 & all_hits_complete$mean_area_rep2 < 600),]$size <- 'small'
all_hits_complete[which(all_hits_complete$mean_area_rep1 > 500 & all_hits_complete$mean_area_rep2 > 600),]$size <- 'large'

all_hits_complete <- all_hits_complete %>%
  select(`Systematic ID`, size)

thing <- all_data %>%
  left_join(all_hits_complete, by = c('ID' = 'Systematic ID'))

dir.create('live_cell_imaging_cp/figures/line_plots')

wt <- all_data %>%
  filter(ID == 'wt')

for(id in unique(thing$ID))
{
  t <- filter(thing, ID == id)
  plot(x = t$Time, y = t$avg_area, col = 'black',
       main = paste0(id, '-', unique(t$size)), type = 'b', ylim = c(0,1000))
  points(x = wt$Time, y = wt$avg_area, col= 'red')
  lines(x = wt$Time, y = wt$avg_area, col = 'red')
  
  file_name = paste0('live_cell_imaging_cp/figures/line_plots/', id, '.jpg')
  
  dev.copy(device = jpeg, file_name)
  dev.off()
}
