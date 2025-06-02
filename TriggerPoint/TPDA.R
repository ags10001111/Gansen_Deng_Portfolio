library(ggplot2)
library(GGally)
library(CCA)
library(CCP)
library(tidyverse)
library(zoo)


df_PD <- read.csv('PD_new.csv')
df_PPT <- read.csv('PPT_new.csv',na.strings='')

df_PD <- df_PD %>% mutate_if(is.character, as.numeric)
df_PD[is.na(df_PD)] <- 0

df_PPT <- df_PPT %>% mutate_if(is.character, as.numeric)
df_PPT <- na.aggregate(df_PPT, FUN = max)

cc1 <- cc(df_PD, df_PPT)
cc1$cor

cc1[3:4]

cc2 <- comput(df_PD, df_PPT, cc1)

# display canonical loadings
cc2[3:6]

# tests of canonical dimensions
rho <- cc1$cor
## Define number of observations, number of variables in first set, and number of variables in the second set.
n <- dim(df_PD)[1]
p <- length(df_PD)
q <- length(df_PPT)

## Calculate p-values using the F-approximations of different test statistics:
p.asym(rho, n, p, q, tstat = "Wilks")
p.asym(rho, n, p, q, tstat = "Hotelling")
p.asym(rho, n, p, q, tstat = "Roy")

### Correlation Matrix between two dataframes-----------------------------------------------------------
cor_mat <- as.data.frame(cor(df_PD, df_PPT))
t_mat <- as.matrix(cor_mat * sqrt(nrow(df_PD) - 2) / sqrt(1 - cor_mat^2))
pv_mat <- 2 * pt(t_mat, df=nrow(df_PD) - 2, lower.tail = TRUE)
pv_df <- as.data.frame(pv_mat)

which(pv_df < 0.05, arr.ind=TRUE)

result_df <- pv_df %>%
  rownames_to_column() %>%
  pivot_longer(cols = -rowname) %>%
  unite(name, rowname, name, sep = " & ")

coef_df <- cor_mat %>%
  rownames_to_column() %>%
  pivot_longer(cols = -rowname) %>%
  unite(name, rowname, name, sep = " & ")

result_df$coef <- coef_df$value
colnames(result_df) <- c('Variable Pair', 'p_value', 'Correlation Coefficient')
result_df1 <- subset(result_df, p_value < 0.05)
write.csv(result_df1, 'Variable Pairs with significant correlation.csv', row.names = F)


# CCA between MTrP Criteria and PPT -----------------------------------------------------
f_PD <- read.csv('PD_new.csv')
df_PPT <- read.csv('PPT_new.csv',na.strings='')

df_PD <- df_PD %>% mutate_if(is.character, as.numeric)
df_PD[is.na(df_PD)] <- 0

df_PPT <- df_PPT %>% mutate_if(is.character, as.numeric)
df_PPT <- na.aggregate(df_PPT, FUN = max)

cc1 <- cc(df_PD, df_PPT)
cc1$cor

cc1[3:4]

cc2 <- comput(df_PD, df_PPT, cc1)

# display canonical loadings
cc2[3:6]

# tests of canonical dimensions
rho <- cc1$cor
## Define number of observations, number of variables in first set, and number of variables in the second set.
n <- dim(df_PD)[1]
p <- length(df_PD)
q <- length(df_PPT)

## Calculate p-values using the F-approximations of different test statistics:
p.asym(rho, n, p, q, tstat = "Wilks")
p.asym(rho, n, p, q, tstat = "Hotelling")
p.asym(rho, n, p, q, tstat = "Roy")

