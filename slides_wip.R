library('tidyverse')
source('functions_CPA.R')
library('pheatmap')
library('gtools')

omics_dataset <- load_omics_data()

tp0 <- filter(omics_dataset, time_point == 0)

normalised_omics <- omics_dataset %>%
  left_join(tp0, by = c('replicate','ID','molecule')) %>%
  transform(read_number.y = read_number.y + 0.001,
            read_number.x = read_number.x +0.001) %>%
  mutate(read_number.y = read_number.x / read_number.y) %>%
  mutate(log2foldchange = log2(read_number.y)) %>%
  rename('read_number.y', 'fold_change') %>%
  transform(time_point.x = as.numeric(time_point.x))

gene_lists <- load_gene_lists()

omics_annot <- normalised_omics %>%
  inner_join(gene_lists, by = c('ID' = 'Systematic ID')) %>%
  filter(type == 'Up RNA' & molecule == 'RNA')

omics_annot$sample_name <- paste(omics_annot$replicate,omics_annot$time_point.x, omics_annot$molecule,sep = '_')


test <- dplyr::select(omics_annot, log2foldchange, ID, sample_name)

hits_for_heatmap <- test %>%
  spread(sample_name, log2foldchange, fill = 0)
hits_for_heatmap <- hits_for_heatmap[, mixedsort(colnames(hits_for_heatmap))]

col <- colorRampPalette(c('blue','grey','yellow'))
t<-pheatmap(dplyr::select(hits_for_heatmap,-ID), cluster_cols = F, labels_row = hits_for_heatmap$ID, color = col(20),
         breaks = seq(-20,20,length.out = 20))

cl<-cutree(t$tree_row, k=3)

clusters <- data.frame(cl = cl, names = hits_for_heatmap$ID)

write_csv(clusters,'upregulated_cl.csv')


omics_annot <- normalised_omics %>%
  inner_join(gene_lists, by = c('ID' = 'Systematic ID')) %>%
  filter(type == 'Down RNA' & molecule == 'RNA')

omics_annot$sample_name <- paste(omics_annot$replicate,omics_annot$time_point.x, omics_annot$molecule,sep = '_')


test <- dplyr::select(omics_annot, log2foldchange, ID, sample_name)

hits_for_heatmap <- test %>%
  spread(sample_name, log2foldchange, fill = 0)
hits_for_heatmap <- hits_for_heatmap[, mixedsort(colnames(hits_for_heatmap))]

col <- colorRampPalette(c('blue','grey','yellow'))
t<-pheatmap(dplyr::select(hits_for_heatmap,-ID), cluster_cols = F, labels_row = hits_for_heatmap$ID, color = col(20),
            breaks = seq(-5,5,length.out = 20), cutree_rows = )

cl<-cutree(t$tree_row, k=3)

clusters <- data.frame(cl = cl, names = hits_for_heatmap$ID)

write_csv(clusters,'upregulated_cl.csv')
