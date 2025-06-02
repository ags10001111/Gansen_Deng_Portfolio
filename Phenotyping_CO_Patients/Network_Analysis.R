library(tidyverse)
load('AFP_NA.RData')

AFP_encode0 <- read_csv('Encode_data.csv')[, -1] # 3-class data
AFP_clean <- read_csv('AFP_clean.csv')[, -1]

con_var <- colnames(AFP_encode)[1:9] # Continuous variables
pain_var <- colnames(AFP_encode)[c(2:4, 7:9)]
cat_var <- colnames(AFP_clean)[c(1, 3, 9:15)]

###Normality Checking###---------------------------------------------------------------------------------------
library(MVN)
mvn(AFP_encode[con_var], mvnTest = 'mardia') # Only Age is normally distributed

# Transform data to let it follows normal distribution
library(bestNormalize)
## Box-Cox's Transformation
boxcox_transform <- function(x){
  x = x + 0.0000001
  result <- boxcox(x)
  return(result$x.t)
}

## Yeo-Johnson's Transformation
yj_transform <- function(x){
  result <- yeojohnson(x)
  return(result$x.t)
}

## Pick up the best transformation automatically
bs_transform <- function(x){
  result <- bestNormalize(x, allow_lambert_s = T, allow_lambert_h = T)
  return(result$x.t)
}

# Transform the data
AFP_transform <- apply(AFP_encode[con_var], 2, yj_transform)

# Use the Mardia's test to check if the data is normally distributed
mvn(AFP_transform, mvnTest = 'mardia')

# Combine continuous and categorical variables
df_na <- cbind(AFP_transform, AFP_clean[cat_var])
df_na[df_na == 'Unknown'] <- -1
df_na[cat_var] <- sapply(df_na[cat_var], as.ordered)

###Estimation and Plotting###-----------------------------------------------------------------------------------
library(bootnet)
library(qgraph)

# The network with only continuous variables
Network_con <-estimateNetwork(AFP_transform, default = "EBICglasso", 
                              threshold = FALSE, corMethod = "cor_auto")
plot(Network_con, layout = 'spring', labels = T)

centralityPlot(Network_con, include = c("Strength","ExpectedInfluence"),
               orderBy = "Strength")

# The network with all variables
Network_full <- estimateNetwork(df_na, default = "EBICglasso",
                                 threshold = FALSE, corMethod = "cor_auto")
plot(Network_full, layout = 'spring', labels = T, label.cex = 2, esize = 10)
print(Network_full)

centralityPlot(Network_full, include = c("Strength", "Closeness", "Betweenness"),
               orderBy = "Strength")

# The network without treatment variables
df_na_wt <- df_na[ , -which(names(df_na) %in% c("Massage Therapy", "Acupuncture", "Chiropractic", 
                                                  "Physiotherapy", "Mental Health"))]
Network_wt <- estimateNetwork(df_na_wt, default = "EBICglasso",
                                threshold = FALSE, corMethod = "cor_auto")
plot(Network_wt, layout = 'spring', labels = T, label.cex = 2, esize = 10)
print(Network_wt)

centralityPlot(Network_wt, include = c("Strength", "Closeness", "Betweenness"),
               orderBy = "Strength")

###Network Accuracy###------------------------------------------------------------------------------------------
# Edge-weight accuracy
boot1 <- bootnet(Network_wt, nBoots = 1000, nCores = 8)
plot(boot1, labels = F, order = "sample")
summary(boot1)

# Centrality stability
boot2 <- bootnet(Network_wt, nBoots = 1000, type = "case", nCores = 8)
plot(boot2)
corStability(boot2)

###Mixed Graphical Model###----------------------------------------------------------------------------------------
library(mgm)
df_na[cat_var] <- sapply(df_na[cat_var], as.numeric)
type <- c(rep('g', 9), rep('c', 8),'g')
cat_num <- c(rep(1, 9), 2, 4, 3, 5, 5, 5, 5, 5, 1)
Network_mix <- mgm(df_na, type, cat_num, lamda.sel="EBIC")

library(qgraph)
qgraph(Network_mix, 
       vsize=3, layout="spring", 
       edge.color = rgb(33,33,33,100, 
                        maxColorValue = 255), 
       border.width=1.5,
       border.color="black",
       nodeNames=colnames(df_na),
       legend=TRUE, 
       legend.mode="groups",
       legend.cex=.75)
?mgmfit




###Including interactions###------------------------------------------------------------------------------------
# Get the top correlations (greater than 0.25)
cor_strength <- Network_wt$results$optnet
table(cor_strength)
which(cor_strength>0.25, arr.ind = T)

# Fit a multinomial regression model
AFP_encode <- read_csv('Discrete_Data.csv')[, -1] # Discretized data
HC_result <- read.csv('HC_result.csv')
LCA_result <- read.csv('LCA_result.csv')

AFP_encode <- AFP_encode0
AFP_encode$Diagnose <- AFP_encode0$Diagnose
#AFP_encode$Diagnose <- HC_result$Cluster_HC
#AFP_encode$Diagnose <- LCA_result$Cluster
#AFP_encode$Diagnose <- ifelse(HC_result$Cluster_HC != LCA_result$Cluster, HC_result$Cluster_HC, 3)

trt = FALSE
cont = TRUE
if (trt & cont){
  library(nnet)
  model_mn <- multinom(Diagnose ~ . + `Employment Status_2.0`*Age + `Employment Status_3.0`*Age + 
                         `Employment Status_Unknown`*Age + `Pain Interference (Mean)`*`Pain severity (Mean)` +
                         `ANX (T-Score)`*`SOM (T-Score)` + `ANX (T-Score)`*`DEP (T-Score)` +
                         `DEP (T-Score)`*`SOM (T-Score)`, data = AFP_encode)
} else if (cont){
  AFP_encode <- AFP_encode[ , -which(names(AFP_encode) %in% c("Massage Therapy", "Acupuncture", "Chiropractic", 
                                                              "Physiotherapy", "Mental Health", "Massage Therapy_0", 
                                                              "Acupuncture_0", "Chiropractic_0", 
                                                              "Physiotherapy_0", "Mental Health_0",
                                                              "Massage Therapy_Unknown", "Acupuncture_Unknown",
                                                              "Chiropractic_Unknown", 
                                                              "Physiotherapy_Unknown", "Mental Health_Unknown"))]
  model_mn <- multinom(Diagnose ~ . + `Employment Status_2.0`*Age + `Employment Status_3.0`*Age + 
                         `Employment Status_Unknown`*Age + `Pain Interference (Mean)`*`Pain severity (Mean)` +
                         `ANX (T-Score)`*`SOM (T-Score)` + `ANX (T-Score)`*`DEP (T-Score)` +
                         `DEP (T-Score)`*`SOM (T-Score)` + (`Pain Interference (Mean)`)^3 +
                         (`Pain severity (Mean)`)^3 + (Age)^3, data = AFP_encode)
} else{
  model_mn <- multinom(Diagnose ~ . + `Employment Status_1.0`*Age_1 + `Employment Status_2.0`*Age_1 + 
                         `Employment Status_3.0`*Age_1 + `Pain Interference (Mean)_2`*`Pain severity (Mean)_2` +
                         `Pain Interference (Mean)_2`*`Pain severity (Mean)_3` + `Pain Interference (Mean)_3`*`Pain severity (Mean)_2` +
                         `Pain Interference (Mean)_3`*`Pain severity (Mean)_3` +
                         `ANX (T-Score)_2`*`SOM (T-Score)_2` +
                         `ANX (T-Score)_2`*`SOM (T-Score)_3` + `ANX (T-Score)_3`*`SOM (T-Score)_2` +
                         `ANX (T-Score)_3`*`SOM (T-Score)_3` +
                         `ANX (T-Score)_2`*`DEP (T-Score)_2` +
                         `ANX (T-Score)_2`*`DEP (T-Score)_3` + `ANX (T-Score)_3`*`DEP (T-Score)_2` +
                         `ANX (T-Score)_3`*`DEP (T-Score)_3` +
                         `SOM (T-Score)_2`*`DEP (T-Score)_2` +
                         `SOM (T-Score)_2`*`DEP (T-Score)_3` + `SOM (T-Score)_3`*`DEP (T-Score)_2` +
                         `SOM (T-Score)_3`*`DEP (T-Score)_3`, data = AFP_encode)
}

# Stepwise regression
model_mn_step <- step(model_mn, k=log(nrow(AFP_encode)))
model_mn_step$coefnames

### Type III ANOVA test
library(afex)
set_sum_contrasts()
library(car)
aov_test <- Anova(model_mn_step, type="III",contrasts=c("contr.sum", "contr.poly"))

### Print out the result in a friendly form
three_class <- T
if (three_class){
  coef_df <- data.frame(cbind(t(coef(model_mn_step)), c(NA, aov_test$`Pr(>Chisq)`)))
  colnames(coef_df) <- c('ceof (2 VS 1)', 'coef (3 VS 1)', 'p-value (LRT)')
} else{
  coef_df <- data.frame(cbind(coef(model_mn_step), c(NA, aov_test$`Pr(>Chisq)`)))
  colnames(coef_df) <- c('coef', 'p-value (LRT)')
}

#coef_combine <- coef_df
rownames(coef_combine) <- coef_combine$Row.names
coef_combine <- merge(coef_combine[,-1], coef_df, by = "row.names", all.x = TRUE)

colnames(coef_df) <- c('coe (2 VS 1)', 'coe (3 VS 1)', 'p-values (LRT)')
rownames(coef_combine1) <- coef_combine1$Row.names
coef_combine <- merge(coef_combine1[,-1], coef_df, by = "row.names", all.x = TRUE)


write.csv(coef_combine, 'Multinomal coeffecients and p-values.csv', row.names = T)



