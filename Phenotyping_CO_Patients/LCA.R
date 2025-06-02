library(tidyverse)
#load('LCA.RData')

###Data Preparation###-------------------------------------------------------------------------------------------
AFP_encode <- read_csv('Encode_data.csv')[, -1] # 3-class data
AFP_clean <- read_csv('AFP_clean.csv')[, -1]
AFP_dis <- read_csv('Discrete_Data.csv')[, -1]
con_var <- colnames(AFP_encode)[1:9] # Continuous variables
pain_var <- colnames(AFP_encode)[c(2:4, 7:9)]
cat_var <- colnames(AFP_clean)[c(1, 3, 9:15)]

df_lca <- as.data.frame(AFP_dis)
df_lca[df_lca == 'Unknown'] <- -1
df_lca[cat_var] <- sapply(df_lca[cat_var], as.numeric)

###Latent Class Analysis (poLCA package)###----------------------------------------------------------------------
library(poLCA)
LCA_res = poLCA(cbind(Sex, Age = Age + 1,
                   `Employment Status` = `Employment Status` + 2, 
                   `Body Diagram` = `Body Diagram` + 1,
                   `Pain severity (Mean)`, `Pain Interference (Mean)`,
                   `# of pain medications`, `# of opioids`,
                   Smoke = Smoke + 2, `Massage Therapy` = `Massage Therapy` + 2,
                   `Acupuncture` = `Acupuncture` + 2,
                   `Chiropractic` = `Chiropractic` + 2,
                   `Physiotherapy` = `Physiotherapy` + 2,
                   `Mental Health` = `Mental Health` + 2,
                   `Pain Stage` = `Pain Stage` + 2,
                   `DEP (T-Score)`, `ANX (T-Score)`, `SOM (T-Score)`
                   ) ~ 1, 
             maxiter=50000, nclass=2, graphs = T, 
             nrep=50, data=df_lca)

#BIC(n=2): 7504.483
#BIC(n=3): 7470.251 Best Model!
#BIC(n=4): 7565.053
#BIC(n=5): 7696.731

# Summarize the cluster profile
df_lca1 <- cbind(df_lca, LCA_res$predclass)
colnames(df_lca1)[19] <- 'Cluster'

library(data.table)
tbl_list1 <- apply(df_lca1[df_lca1$Cluster == 1, ], 2, table)
tmp1 <- rbindlist(lapply(tbl_list1, function(x) as.data.frame.list(x)), fill=TRUE)[-19,]
colnames(tmp1) <- c('1', '2', '0', '-1', '3', '4')
tmp1$Variable <- colnames(df_lca1)[-19]
tmp1$Cluster <- '1'

tbl_list2 <- apply(df_lca1[df_lca1$Cluster == 2, ], 2, table)
tmp2 <- rbindlist(lapply(tbl_list2, function(x) as.data.frame.list(x)), fill=TRUE)[-19,]
colnames(tmp2) <- c('1', '2', '0', '-1', '3', '4')
tmp2$Variable <- colnames(df_lca1)[-19]
tmp2$Cluster <- '2'

tbl_list3 <- apply(df_lca1[df_lca1$Cluster == 3, ], 2, table)
tmp3 <- rbindlist(lapply(tbl_list3, function(x) as.data.frame.list(x)), fill=TRUE)[-19,]
colnames(tmp3) <- c('1', '2', '0', '-1', '3', '4')
tmp3$Variable <- colnames(df_lca1)[-19]
tmp3$Cluster <- '3'

tmp <- cbind(tmp1, tmp2, tmp3)

# Another plot
lcmodel <- reshape2::melt(LCA_res$probs, level=2)
zp1 <- ggplot(lcmodel,aes(x = L2, y = value, fill = Var2))
zp1 <- zp1 + geom_bar(stat = "identity", position = "stack")
zp1 <- zp1 + facet_grid(Var1 ~ .) 
zp1 <- zp1 + scale_fill_brewer(type="seq", palette="Greys") +theme_bw()
zp1 <- zp1 + labs(x = "Fragebogenitems",y="Anteil der Item-\nAntwortkategorien", fill ="Antwortkategorien")
zp1 <- zp1 + theme( axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),                    
                    panel.grid.major.y=element_blank())
zp1 <- zp1 + guides(fill = guide_legend(reverse=TRUE))
print(zp1)

###Latent Class Analysis (Exclude treatment variables)###------------------------------------------------------
df_lca_wt <- df_lca[ , -which(names(df_lca) %in% c("Massage Therapy", "Acupuncture", "Chiropractic", 
                                               "Physiotherapy", "Mental Health"))]
library(poLCA)
LCA_res1 = poLCA(cbind(Sex, Age = Age + 1,
                      `Employment Status` = `Employment Status` + 2, 
                      `Body Diagram` = `Body Diagram` + 1,
                      `Pain severity (Mean)`, `Pain Interference (Mean)`,
                      `# of pain medications`, `# of opioids`,
                      Smoke = Smoke + 2,
                      `Pain Stage` = `Pain Stage` + 2,
                      `DEP (T-Score)`, `ANX (T-Score)`, `SOM (T-Score)`
) ~ 1, 
maxiter=50000, nclass=3, graphs = F, 
nrep=50, data=df_lca_wt)

#BIC(n=2): 4748.843 / 4554.835(AIC) Best Model in BIC!
#BIC(n=3): 4797.977 / 4505.321
#BIC(n=4): 4869.584 / 4478.28
#BIC(n=5): 4949.175 / 4459.223  Best Model in AIC!
#BIC(n=6): 5051.743 / 4463.143
#BIC(n=7): 5172.713 / 4471.506
#BIC(n=8): 5291.777 / 4491.607
#BIC(n=9): 5402.482 / 4517.938
#BIC(n=10): 5526.395 / 4543.203

# Summarize the cluster profile
df_lca1 <- cbind(df_lca_wt, LCA_res1$predclass)
colnames(df_lca1)[14] <- 'Cluster'
write.csv(df_lca1, '234.csv')

library(data.table)
tbl_list1 <- apply(df_lca1[df_lca1$Cluster == 1, ], 2, table)
tmp1 <- rbindlist(lapply(tbl_list1, function(x) as.data.frame.list(x)), fill=TRUE)[-14,]
colnames(tmp1) <- c('1', '2', '0', '-1', '3', '4')
tmp1$Variable <- colnames(df_lca1)[-14]
tmp1$Cluster <- '1'

tbl_list2 <- apply(df_lca1[df_lca1$Cluster == 2, ], 2, table)
tmp2 <- rbindlist(lapply(tbl_list2, function(x) as.data.frame.list(x)), fill=TRUE)[-14,]
colnames(tmp2) <- c('1', '2', '0', '-1', '3', '4')
tmp2$Variable <- colnames(df_lca1)[-14]
tmp2$Cluster <- '2'

tbl_list3 <- apply(df_lca1[df_lca1$Cluster == 3, ], 2, table)
tmp3 <- rbindlist(lapply(tbl_list3, function(x) as.data.frame.list(x)), fill=TRUE)[-14,]
colnames(tmp3) <- c('1', '2', '0', '-1', '3', '4')
tmp3$Variable <- colnames(df_lca1)[-14]
tmp3$Cluster <- '3'

tmp <- cbind(tmp1, tmp2, tmp3)
write.csv(tmp, 'LCA_result(WT).csv')

###Cluster profile (OF HC)###----------------------------------------------------------------------------------
library(readxl)
df_of <- read_xlsx('Cluster Result(OF).xlsx')[, -1]
df_of[df_of$Cluster_HC == 0, 19] <- 3

library(data.table)
tbl_list1 <- apply(df_of[df_of$Cluster_HC == 1, ], 2, table)
tmp1 <- rbindlist(lapply(tbl_list1, function(x) as.data.frame.list(x)), fill=TRUE)[-19,]
colnames(tmp1) <- c('1', '2', '0', '-1', '3', '4')
tmp1$Cluster <- '1'

tbl_list3 <- apply(df_of[df_of$Cluster_HC == 3, ], 2, table)
tmp3 <- rbindlist(lapply(tbl_list3, function(x) as.data.frame.list(x)), fill=TRUE)[-19,]
colnames(tmp3) <- c('1', '2', '0', '-1', '3', '4')
tmp3$Cluster <- '3'

tbl_list2 <- apply(df_of[df_of$Cluster_HC == 2, ], 2, table)
tmp2 <- rbindlist(lapply(tbl_list2, function(x) as.data.frame.list(x)), fill=TRUE)[-19,]
colnames(tmp2) <- c('1', '2', '0', '-1', '3', '4')
tmp2$Cluster <- '2'

tmp <- cbind(tmp1, tmp2, tmp3)

write.csv(tmp, '111.csv')

table(df_of$Cluster_HC, df_lca1$Cluster)

###Cluster profile (OF HC)[Without treatment variables]###-----------------------------------------------------
library(readxl)
df_of <- read_xlsx('Cluster Result_WT(OF)3.xlsx')[, -1]
#df_of[df_of$Cluster_HC == 1, 19] <- 2
#df_of[df_of$Cluster_HC == 0, 19] <- 1

library(data.table)
tbl_list1 <- apply(df_of[df_of$Cluster_HC == 1, ], 2, table)
tmp1 <- rbindlist(lapply(tbl_list1, function(x) as.data.frame.list(x)), fill=TRUE)[-19,]
colnames(tmp1) <- c('1', '2', '0', '-1', '3', '4')
tmp1$Cluster <- '1'

tbl_list2 <- apply(df_of[df_of$Cluster_HC == 2, ], 2, table)
tmp2 <- rbindlist(lapply(tbl_list2, function(x) as.data.frame.list(x)), fill=TRUE)[-19,]
colnames(tmp2) <- c('1', '2', '0', '-1', '3', '4')
tmp2$Cluster <- '2'

tbl_list3 <- apply(df_of[df_of$Cluster_HC == 3, ], 2, table)
tmp3 <- rbindlist(lapply(tbl_list3, function(x) as.data.frame.list(x)), fill=TRUE)[-19,]
colnames(tmp3) <- c('1', '2', '0', '-1', '3', '4')
tmp3$Cluster <- '3'

tmp <- cbind(tmp1, tmp2, tmp3)
write.csv(tmp, 'HC_result(WT).csv')

###Profile of mismatched patienets in two clusters###----------------------------------------------------------
df_of$Cluster <- df_lca1$Cluster
tbl_list12 <- apply(df_of[df_of$Cluster_HC == 1 & df_of$Cluster == 2, ], 2, table)
tmp12 <- rbindlist(lapply(tbl_list12, function(x) as.data.frame.list(x)), fill=TRUE)[-(19:20),]
colnames(tmp12) <- c('1', '2', '0', '-1', '3', '4')
tmp12$Variable <- colnames(df_of)[-(19:20)]

###Cluster profile (Correlation Distance)[Without treatment variables]###-----------------------------------------------------
df_of <- read_csv('corr_cluster.csv')[, -1]
#df_of[df_of$Cluster == 1, 19] <- 2
#df_of[df_of$Cluster == 0, 19] <- 1

library(data.table)
tbl_list1 <- apply(df_of[df_of$Cluster == 1, ], 2, table)
tmp1 <- rbindlist(lapply(tbl_list1, function(x) as.data.frame.list(x)), fill=TRUE)[-19,]
colnames(tmp1) <- c('1', '2', '0', '-1', '3', '4')
tmp1$Cluster <- '1'

tbl_list2 <- apply(df_of[df_of$Cluster == 2, ], 2, table)
tmp2 <- rbindlist(lapply(tbl_list2, function(x) as.data.frame.list(x)), fill=TRUE)[-19,]
colnames(tmp2) <- c('1', '2', '0', '-1', '3', '4')
tmp2$Cluster <- '2'

tbl_list3 <- apply(df_of[df_of$Cluster == 3, ], 2, table)
tmp3 <- rbindlist(lapply(tbl_list3, function(x) as.data.frame.list(x)), fill=TRUE)[-19,]
colnames(tmp3) <- c('1', '2', '0', '3')
tmp3$Cluster <- '3'

tmp <- cbind(tmp1, tmp2, tmp3)
write.csv(tmp, 'HC_result(correlation).csv')
