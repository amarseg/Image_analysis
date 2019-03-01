library('tidyverse')
source('functions_CPA.R')
library('pheatmap')
library('gtools')

common_maha_hits <- read_csv('common_hits_z_Score.csv')

########################Hits in common with omics data######################

omics_data <- load_gene_lists()

common_hits <- inner_join(common_maha_hits, omics_data, by = c('x' = 'Systematic ID'))

##########################Plot hits in both replicates######################

summ_1 <- read_csv('output_rep1/summary_stats_rep1.csv')
summ_2 <- read_csv('output_rep2/statistics_rep2.csv')

all_summ <- inner_join(summ_1, summ_2, by = 'Systematic ID', suffix = c('_rep1','_rep2')) %>%
  add_column(hits = ifelse(.$`Systematic ID` %in% common_maha_hits$x,
                           'hit',
                           NA))

ggplot(all_summ, aes(x = mean_area_rep1, y = mean_area_rep2, color = hits)) +
  geom_point() +
  geom_smooth(method = 'lm')

ggplot(filter(all_summ,!is.na(hits)), aes(x = cv_rep1, y = cv_rep2, colour = hits)) +
  geom_point()


ggplot(all_summ, aes(x = mean_area_rep1, y = cv_rep1, colour = hits)) +
  geom_point()

ggplot(all_summ, aes(x = mean_area_rep2, y = cv_rep2, colour = hits)) +
  geom_point() +
  ylim(0.25,2)

ggplot(filter(all_summ, !is.na(hits)), aes(x = mean_area_rep1,colour = hits)) +
  geom_histogram() 

ggplot(filter(all_summ, !is.na(hits)), aes(x = mean_area_rep2,colour = hits)) +
  geom_histogram() 
ggplot(filter(all_summ, !is.na(hits)), aes(x = mean_area_rep1,colour = hits)) +
  geom_histogram() 

ggplot(filter(all_summ, !is.na(hits)), aes(x = mean_area_rep2,colour = hits)) +
  geom_histogram() 
#######################How are these guys behaving in the omics dataset?##################
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



hits_omics <- inner_join(normalised_omics, common_maha_hits, by = c('ID' = 'x') ) %>%
  filter(!is.na(log2foldchange) & !is.infinite(log2foldchange) & !is.nan(log2foldchange))


ggplot(hits_omics, aes(x = time_point.x, y = log2foldchange)) +
  geom_point(color = 'gray')+
  geom_line(color = 'gray', aes(group = interaction(ID,replicate)))+
  stat_summary( fun.y=mean, geom="line", colour="red", size = 2) +
  facet_grid(~molecule)



hits_omics$sample_name <- paste(hits_omics$replicate,hits_omics$time_point.x, hits_omics$molecule,sep = '_')

test <- select(hits_omics, log2foldchange, ID, sample_name)

hits_for_heatmap <- test %>%
  spread(sample_name, log2foldchange, fill = 0)

hits_for_heatmap <- hits_for_heatmap[, mixedsort(colnames(hits_for_heatmap))]

prot <- hits_for_heatmap[, grep(colnames(hits_for_heatmap), pattern = 'Prot|ID')] 
prot <- prot[,!is.na(prot)]

rna <- hits_for_heatmap[, grep(colnames(hits_for_heatmap), pattern = 'RNA|ID')]

pheatmap(select(prot,-ID), cluster_cols = F, labels_row = prot$ID)

rna_cl <- pheatmap(select(rna,-ID), cluster_cols = F, labels_row = rna$ID,
         kmeans_k = 2)
rna_clusters <- data.frame(id = rna$ID, cl = rna_cl$kmeans$cluster) %>%
  write.csv('rna_clusters.csv')

##########################These are my arbitrary thresholds#####################
size_thr <- 500

all_summ$size <- NA
all_summ[which(all_summ$mean_area_rep1 > 500 & all_summ$mean_area_rep2 >700),]$size <- 'large'
all_summ[which(all_summ$mean_area_rep1 < 500 & all_summ$mean_area_rep2 <700),]$size <- 'small'

ggplot(filter(all_summ, !is.na(hits)), aes(x = mean_area_rep1, y = mean_area_rep2, colour = size)) +
  geom_point() 

write_csv(filter(all_summ, !is.na(hits)),'z_score_hits.csv')

size_hits <- filter(all_summ, !is.na(hits))


annot_hits_omics <- left_join(hits_omics, all_summ, by = c('ID'= 'Systematic ID'))

ggplot(filter(annot_hits_omics, !is.na(size)), aes(x = time_point.x, y = log2foldchange)) +
  geom_point(color = 'gray')+
  geom_line(color = 'gray', aes(group = interaction(ID,replicate)))+
  stat_summary( fun.y=mean, geom="line", colour="red", size = 2) +
  facet_grid(~molecule*size)
ggsave('hits_in_omics.pdf')
##########################Use jackies dataset#######################################
jackie_pheno <- read_csv('deletion_phenotypes.csv') %>%
  inner_join(size_hits, by = 'Systematic ID')

ggplot(filter(jackie_pheno,!is.na(size)) , aes(x = size, fill = `Phenotypic classification used for analysis`))+
  geom_bar() +
  theme_bw()
ggsave('figs/jackie_summary_maha.pdf')

ggplot(jackie_pheno , aes(x = size, fill = `Deletion mutant phenotype description`))+
  geom_bar() +
  theme_bw()
ggsave('figs/deletion_description_phenotype_maha.pdf')

