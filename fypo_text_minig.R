rm(list = ls())
library(tidyverse)
library(tm)
library(igraph)
library(SnowballC)
load('lists.011018.rda')
library(pheatmap)
library(gtools)
source('functions_CPA.R')
hits_summary <- stack(out)
library(wordcloud)

fypo_db <- fypo_database_loading() %>%
  inner_join(hits_summary, by = c(`Gene systematic ID` = 'values'))

corpus <- Corpus(VectorSource(fypo_db$Definition))
corpus <- tm_map(corpus, removeWords, c("cell", "vegetative","increased",'decreased','during','population','viable','normal','growth',
                                        'sensitivity','morphology','cellular','level','sensitive','with')) 

corpus <- tm_map(corpus, removeWords, stopwords("english"))

dtm <- TermDocumentMatrix(corpus)
dtm <- removeSparseTerms(dtm, sparse=0.99)

m <- as.matrix(dtm)
v <- sort(rowSums(m),decreasing=TRUE)
d <- data.frame(word = names(v),freq=v)
set.seed(1234)
wordcloud(words = d$word, freq = d$freq, min.freq = 1,
          max.words=200, random.order=FALSE, rot.per=0.35, 
          colors=brewer.pal(8, "Dark2"))


findAssocs(dtm, terms = "abnormal", corlimit = 0.1)
findAssocs(dtm, terms = "resistance", corlimit = 0.1)
findAssocs(dtm, terms = "carbon", corlimit = 0.1)



# change it to a Boolean matrix
m[m>=1] <- 1
# transform into a term-term adjacency matrix
termMatrix <- m %*% t(m)

g <- graph.adjacency(termMatrix, weighted = T, mode = 'undirected')
# remove loops
g <- igraph::simplify(g)
# set labels and degrees of vertices
V(g)$label <- V(g)$name
V(g)$degree <- degree(g)

# set seed to make the layout reproducible
layout1 <- layout.fruchterman.reingold(g)
plot(g, layout=layout1, vertex.size=20, 
vertex.label.color="darkred")


source('Z:/Pers_Amalia/SGA_analysis/enrichment_functions.R')

####enrichment######

go_db<-load_go()
go_enrich <- enricher(out$var, TERM2GENE = go_db$term2gene, TERM2NAME = go_db$term2name)
dotplot(go_enrich)
write_csv(as.tibble(go_enrich),path = 'data/enrichment_GO.csv')