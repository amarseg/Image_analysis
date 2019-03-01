library('tidyverse')

 hits <-read_csv('z_score_hits.csv') %>%
   filter(!is.na(Plate))
 