###Calculate mahalanobis distance, comparison with z_score predictions
library(tidyverse)
source('functions_CPA.R')
library('clusterProfiler')

areas <- read_csv('output_rep1/cell_areas.csv')
summary_stats <- read_csv('output_rep1/summary_stats_rep1.csv') %>%
  mutate(z_score = (mean_area - mean(mean_area))/sd(mean_area))

summary_stats <- add_column(summary_stats, maha_dist = mahalanobis(summary_stats[,c('mean_area','cv')], 
                                     center = colMeans(summary_stats[,c('mean_area','cv')]),
                                     cov = cov(summary_stats[,c('mean_area','cv')])))

pval_thr <- 0.95
t <- sqrt(-2*log(1-pval_thr))

summary_stats <- add_column(summary_stats, outlier = ifelse(summary_stats$maha_dist > t, 'yes','no'))

ggplot(summary_stats, aes(x = mean_area, y = cv, colour = outlier)) +
  geom_point()
ggsave('figs/mahalanohibs_hits_rep1.pdf')
###Do k_means to separate in 4 clusters
maha_hits <- filter(summary_stats, outlier == 'yes')

cl <- kmeans(x = select(maha_hits, mean_area, cv), centers  = 3)

maha_hits$cluster=factor(cl$cluster)
centers=as.data.frame(cl$centers)

ggplot(maha_hits, aes(x = mean_area, y = cv, colour = cluster)) +
  geom_point()#does not work, keep the code to show
ggsave('figs/clusters_maha_rep1.pdf')

ggplot(maha_hits, aes(x = mean_area)) +
  geom_histogram(aes(y = ..density..)) +
  geom_density()

ggplot(maha_hits, aes(x =cv)) +
  geom_histogram(aes(y = ..density..)) +
  geom_density()

write_csv(maha_hits,'output_rep1/maha_hits.csv')

enrichKEGG(maha_hits[which(maha_hits$cluster == 1),]$`Systematic ID`, organism = 'spo')
enrichKEGG(maha_hits[which(maha_hits$cluster == 2),]$`Systematic ID`, organism = 'spo')
enrichKEGG(maha_hits[which(maha_hits$cluster == 3),]$`Systematic ID`, organism = 'spo')

source('Z:/Pers_Amalia/SGA_analysis/enrichment_functions.R')

go_db<-load_go()
go_enrich <- enricher(maha_hits[which(maha_hits$cluster == 1),]$`Systematic ID`, TERM2GENE = go_db$term2gene, TERM2NAME = go_db$term2name)
dotplot(go_enrich)

go_enrich <- enricher(maha_hits[which(maha_hits$cluster == 2),]$`Systematic ID`, TERM2GENE = go_db$term2gene, TERM2NAME = go_db$term2name)
dotplot(go_enrich)

go_enrich <- enricher(maha_hits[which(maha_hits$cluster == 3),]$`Systematic ID`, TERM2GENE = go_db$term2gene, TERM2NAME = go_db$term2name)
dotplot(go_enrich)
###We will classify hits in two axis, size and variability

#################
#Let's do the same with the second replicate and see what they have in common
areas2 <- read_csv('output_rep2/areas_rep2.csv')
summary_stats2 <- read_csv('output_rep2/statistics_rep2.csv') %>%
  mutate(z_score = (mean_area - mean(mean_area))/sd(mean_area))

summary_stats2 <- add_column(summary_stats2, maha_dist = mahalanobis(summary_stats2[,c('mean_area','cv')], 
                                                                   center = colMeans(summary_stats2[,c('mean_area','cv')]),
                                                                   cov = cov(summary_stats2[,c('mean_area','cv')])))
summary_stats2 <- add_column(summary_stats2, outlier = ifelse(summary_stats2$maha_dist > t, 'yes','no'))
ggplot(summary_stats2, aes(x = mean_area, y = cv, colour = outlier)) +
  geom_point()
ggsave('figs/mahalanohibs_hits_rep2.pdf')

maha_hits2 <- filter(summary_stats2, outlier == 'yes')
common_hits<- intersect(maha_hits2$`Systematic ID`, maha_hits$`Systematic ID`)
write.csv(common_hits,'maha_common_hits.csv')
