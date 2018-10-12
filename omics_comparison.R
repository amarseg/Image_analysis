source('functions_CPA.R')
library(plyr)
library(tidyverse)
library(gridExtra)
load('lists.011018.rda')
library(pheatmap)
library(gtools)
hits_summary <- stack(out)

###plot omics hits in the screening data
area <- read_csv('output_rep1/cell_areas.csv')
summary <- read_csv('output_rep1/summary_stats_rep1.csv')

omics_lists <- load_gene_lists()

complete_summary <- left_join(summary, omics_lists, by = 'Systematic ID')

ggplot(complete_summary, aes(x = mean_area, y = cv, color = type)) +
  geom_point()
ggsave('figs/omics_summary.pdf')


complete_area <- left_join(area, omics_lists, by = 'Systematic ID') %>%
  filter(!is.na(type) | ind == 'wt') 

complete_area[which(complete_area$ind == 'wt'),]$type <- 'wt'

p1<-ggplot(subset(complete_area,type == 'wt' | type == 'Up RNA'), 
       aes(x=reorder(`Systematic ID`, AreaShape_Area, mean), 
           y = AreaShape_Area, 
           color = type)) +
  geom_boxplot() +
  ylim(100, 1000) 


p2<-ggplot(subset(complete_area,type == 'wt' | type == 'Down RNA'), 
       aes(x=reorder(`Systematic ID`, AreaShape_Area, mean), 
           y = AreaShape_Area, 
           color = type)) +
  geom_boxplot() +
  ylim(100, 1000) 

p3<-ggplot(subset(complete_area,type == 'wt' | type == 'Up prot'), 
       aes(x=reorder(`Systematic ID`, AreaShape_Area, mean), 
           y = AreaShape_Area, 
           color = type)) +
  geom_boxplot() +
  ylim(100, 1000) 

p4<-ggplot(subset(complete_area,type == 'wt' | type == 'Down prot'), 
       aes(x=reorder(`Systematic ID`, AreaShape_Area, mean), 
           y = AreaShape_Area, 
           color = type)) +
  geom_boxplot() +
  ylim(100, 1000) 

ggsave('figs/omics_avg_area.pdf', marrangeGrob(p1,p2,p3,p4, ncol = 1, nrow = 2))

###plot screening hits in the omics data
omics_dataset <- load_omics_data()

tp0 <- filter(omics_dataset, time_point == 0)

normalised_omics <- omics_dataset %>%
  left_join(tp0, by = c('replicate','ID','molecule')) %>%
  mutate(read_number.y = read_number.x / read_number.y) %>%
  mutate(log2foldchange = log2(read_number.y)) %>%
  rename('read_number.y', 'fold_change') %>%
  transform(time_point.x = as.numeric(time_point.x))



hits_omics <- inner_join(normalised_omics, hits_summary, by = c('ID' = 'values') ) %>%
  filter(!is.na(log2foldchange) & !is.infinite(log2foldchange) & !is.nan(log2foldchange))


ggplot(hits_omics, aes(x = time_point.x, y = log2foldchange)) +
  geom_point(color = 'gray')+
  geom_line(color = 'gray', aes(group = interaction(ID,replicate)))+
  stat_summary( fun.y=mean, geom="line", colour="red", size = 2) +
  facet_grid(~ind*molecule)


all_omics_hits <- left_join(normalised_omics, hits_summary, by = c('ID' = 'values') ) %>%
  filter(!is.na(log2foldchange) & !is.infinite(log2foldchange) & !is.nan(log2foldchange))

all_omics_hits$transparency <- ifelse(is.na(all_omics_hits$ind),yes = 0.2, no = 1)

ggplot(all_omics_hits, aes(x = time_point.x, y = log2foldchange, color = ind, alpha = transparency)) +
  geom_line(data = subset(all_omics_hits, !is.na(ind)),aes(group = interaction(ID,replicate))) +
  stat_summary( fun.y=mean, geom="line", colour="black", size = 1.25, linetype = 'dashed') +
  facet_wrap(~replicate*molecule)
ggsave('figs/omics_hits.pdf', scale = 2)



##Let's do some clustering?
set.seed(25)
cl_data <- all_omics_hits %>%
  select(- read_number.x,- read_number.y, -time_point.y, -transparency) %>%
  unite(temp, replicate, time_point.x) %>%
  spread(key = temp, value = log2foldchange) %>%
  filter(!is.na(ind))

rna <- subset(cl_data, molecule == 'RNA')

rna<-rna[,mixedorder(colnames(rna)),]

rna_cl <-pheatmap(select(rna, -ID,-molecule, -ind), 
                  cluster_cols = F, 
                  kmeans_k = 4,
                  annotation_row = select(rna, ind))

out_rna <- rna %>%
  add_column(cluster = rna_cl$kmeans$cluster) %>%
  write_csv('output_rep1/rna_clusters.csv')

prot <- subset(cl_data, molecule == 'Protein')
prot<-prot[,mixedorder(colnames(prot)),]

prot_cl <- pheatmap(select(prot, -ID, -molecule, -ind), 
                    cluster_cols = F, 
                    kmeans_k = 2, 
                    annotation_row = select(prot, ind))

out_pro <- prot %>%
  add_column(cluster = prot_cl$kmeans$cluster) %>%
  write_csv('output_rep1/prot_clusters.csv')

#plot of pfk1
ggplot(subset(all_omics_hits, ID == 'SPBC16H5.02'),aes(x = time_point.x, y = log2foldchange)) +
  facet_wrap(~molecule) +
  geom_line(aes(group = interaction(ID,replicate)))
ggsave('figs/pfk1_omics.pdf')


ggplot(subset(all_omics_hits, ID == 'SPAC10F6.16'),aes(x = time_point.x, y = log2foldchange)) +
  facet_wrap(~molecule) +
  geom_line(aes(group = interaction(ID,replicate)))
ggsave('figs/igo1_omics.pdf')


ggplot(subset(all_omics_hits, ID == 'SPCC1322.08'),aes(x = time_point.x, y = log2foldchange)) +
  facet_wrap(~molecule) +
  geom_line(aes(group = interaction(ID,replicate)))
ggsave('figs/srk1_omics.pdf')

ggplot(subset(all_omics_hits, ID == 'SPBC31F10.13c'),aes(x = time_point.x, y = log2foldchange)) +
  facet_wrap(~molecule) +
  geom_line(aes(group = interaction(ID,replicate)))
ggsave('figs/hip1_omics.pdf')


###phenotype summary
fypo <- read_csv('deletion_phenotypes.csv') %>%
  inner_join(hits_summary, by = c(`Systematic ID`= 'values'))

ggplot(fypo, aes(x = ind, fill = `Phenotypic classification used for analysis`))+
  geom_bar() 

ggplot(fypo, aes(x = ind, fill = `Deletion mutant phenotype description`))+
  geom_bar() 


###

