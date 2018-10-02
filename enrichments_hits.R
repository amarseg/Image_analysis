library(tidyverse)
library('clusterProfiler')

load('lists.011018.rda')

cls <- compareCluster(out, fun = 'enrichKEGG', org = 'spo')
write.table(out$up, file = 'large.csv', sep = ',')
write.table(out$dn, file = 'small.csv', sep = ',')
write.table(out$var, file = 'variable.csv', sep = ',')
