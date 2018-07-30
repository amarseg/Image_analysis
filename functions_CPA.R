load_CPA_data <- function(per_object_data, per_image_data)
{
  object_data <- read_csv(per_object_data)
  well_data <- 
    read_csv(per_image_data) %>%
    select(ImageNumber, Image_Metadata_QCFlag,Image_FileName_DNA) %>%
    separate(Image_FileName_DNA, into = c('Well','Well_n','Picture','Z_axis','Time','Type'), sep = '--')
  
  all_data <- left_join(object_data, well_data, by = c('ImageNumber' = 'ImageNumber'))
  return(all_data)
}