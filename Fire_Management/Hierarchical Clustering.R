setwd("C:/Users/ags10001111/OneDrive/2020 Winter/Statistical Consulting/Major project/Clustering")
load('Hierarchical Clustering.RData')
library(tidyverse)

#-------------------------------------------------------
# Read the data
NWR_fire = read.csv('Clean_data1.csv')
ggplot(NWR_fire, aes(OBJECTIVE)) + geom_bar() + ggtitle('Fire cases count by OBJECTIVE')

table(NWR_fire[NWR_fire$GENERAL_CAUSE == 'LTG',]$OBJECTIVE)
# Extract the variables to consider
cluster_variable = c('OBJECTIVE', 'INIT_REP_GRP', 'INIT_ATT_TYPE',
                     'LOC_ATTACK', 'GROUND_FORCE', 'AIR_TANKERS',
                     'AIRCRAFT')
NWR_fire_cvar = NWR_fire[, which(names(NWR_fire) %in% cluster_variable)]
#---------------------------------------------------------------
# Do clustering using the cluster package
## Compute the Gower's distance for mixed data
library(StatMatch)
gdist = as.dist(gower.dist(NWR_fire_cvar))
gdist[is.na(gdist)] = max(gdist, na.rm = T)

library(cluster)
## methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

## function to compute coefficient
ac <- function(x) {
  agnes(gdist, diss = T, method = x)$ac
}

## Find the best linkage method
ac_methods = map_dbl(m, ac)
#-----------------------------------------------------------------
# Do clustring using hclust function
## Hierarchical Clustering
NWR_cluster = hclust(gdist, method = 'ward.D')

## Plot the dendrogram
library(ggplot2)
library(ggdendro)
ggdendrogram(NWR_cluster, rotate = TRUE, theme_dendro = FALSE) + 
  xlim('') + xlab('Observations') + ylab("Gower's Distance")

## Get the membership of each observation
group3 = cutree(NWR_cluster, 3)
group4 = cutree(NWR_cluster, 4)

write.csv(group3, 'result.csv')
group31 = read.csv('cluster.csv')
group32 = read.csv('result4.csv')

group31$x1 = ifelse(group31[, 2] == 1, 2, ifelse(group31[, 2] == 2, 3, 1))
mean(group31[, 3] == group3)

group32$x1 = ifelse(group32[, 4] == 1, 1, ifelse(group32[, 4] == 2, 2, 3))

hc_kproto = mean(group31[, 3] == group3)
hc_pam = mean(group32[, 5] == group3)
pam_kproto = mean(group32$x1 == group31$x1)

data.frame(hc_kproto, hc_pam, pam_kproto)



