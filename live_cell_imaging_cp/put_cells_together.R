library(tidyverse)

####Create lineages######
center_thr <- 10

tp0 <- read_csv('live_cell_imaging_cp/output_data_cleaned.csv') %>%
  #group_by(Metadata_Plate_Name, Position, `Systematic ID`) %>% 
  filter(Time == 0) %>%
  mutate(cell_n = ObjectNumber)

data <- read_csv('live_cell_imaging_cp/output_data_cleaned.csv') %>%
  #group_by(Metadata_Plate_Name, Position, `Systematic ID`) %>% 
  filter(Time > 0) %>%
  add_column(cell_n = 0)

for(i in 1:nrow(tp0))
{
  tp_cell <- tp0[i,]
  
  if (nrow(data[which(data$Metadata_Plate_Name == tp_cell$Metadata_Plate_Name &
               data$Position == tp_cell$Position &
               data$`Systematic ID` == tp_cell$`Systematic ID`& 
               data$Location_Center_X > tp_cell$Location_Center_X - center_thr & 
               data$Location_Center_X < tp_cell$Location_Center_X + center_thr &
               data$Location_Center_Y > tp_cell$Location_Center_Y - center_thr &
               data$Location_Center_Y < tp_cell$Location_Center_Y + center_thr),]) >= 1)
  {
    
    
    n = nrow(data[which(data$Metadata_Plate_Name == tp_cell$Metadata_Plate_Name &
                          data$Position == tp_cell$Position &
                          data$`Systematic ID` == tp_cell$`Systematic ID`& 
                          data$Location_Center_X > tp_cell$Location_Center_X - center_thr & 
                          data$Location_Center_X < tp_cell$Location_Center_X + center_thr &
                          data$Location_Center_Y > tp_cell$Location_Center_Y - center_thr &
                          data$Location_Center_Y < tp_cell$Location_Center_Y + center_thr),])
    print(n)
    data[which(data$Metadata_Plate_Name == tp_cell$Metadata_Plate_Name &
                 data$Position == tp_cell$Position &
                 data$`Systematic ID` == tp_cell$`Systematic ID`& 
                 data$Location_Center_X > tp_cell$Location_Center_X - center_thr & 
                 data$Location_Center_X < tp_cell$Location_Center_X + center_thr &
                 data$Location_Center_Y > tp_cell$Location_Center_Y - center_thr &
                 data$Location_Center_Y < tp_cell$Location_Center_Y + center_thr),]$cell_n <- rep(tp_cell$ObjectNumber, n)
  }
    
    
}

track_data <- bind_rows(data, tp0) %>%
  write_csv('live_cell_imaging_cp/track_cell_data.csv')
