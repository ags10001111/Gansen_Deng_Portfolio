setwd("C:/Users/ags10001111/OneDrive/2020 Summer/Master Project/Correlation")
load("Correlation.RData")
library(tidyverse)
#------------------------------------------------------------------
# READ DATA
## Data importing function
library(tidyverse)
read_longitudinal = function(file_name){
  # Specify the data type
  col_type0 = "ffnddfnn"
  col_typeA = paste(strrep("_",8), strrep('d', 55), sep = "")
  col_typeB = paste(strrep("_",63), strrep('d', 55), sep = "")
  col_typeC = paste(strrep("_",118), strrep('d', 55), sep = "")
  col_typeD = paste(strrep("_",173), strrep('d', 55), sep = "")
  col_typeE = paste(strrep("_",228), strrep('d', 55), sep = "")
  col_typeF = paste(strrep("_",283), strrep('d', 55), sep = "")
  
  # Read the data of different timepoints
  data0 = read_csv(file_name, col_types = col_type0, skip = 1)
  dataA = read_csv(file_name, col_types = col_typeA, skip = 1)
  dataB = read_csv(file_name, col_types = col_typeB, skip = 1)
  dataC = read_csv(file_name, col_types = col_typeC, skip = 1)
  dataD = read_csv(file_name, col_types = col_typeD, skip = 1)
  dataE = read_csv(file_name, col_types = col_typeE, skip = 1)
  dataF = read_csv(file_name, col_types = col_typeF, skip = 1)
  
  ## Make the colnames consistent
  colnames(dataB) = colnames(dataA)
  colnames(dataC) = colnames(dataA)
  colnames(dataD) = colnames(dataA)
  colnames(dataE) = colnames(dataA)
  colnames(dataF) = colnames(dataA)
  
  # Combine the data.frame of different timepoints
  data_A = cbind(data0, dataA)
  data_A$Timepoint = 'A'
  
  data_B = cbind(data0, dataB)
  data_B$Timepoint = 'B'
  
  data_C = cbind(data0, dataC)
  data_C$Timepoint = 'C'
  
  data_D = cbind(data0, dataD)
  data_D$Timepoint = 'D'
  
  data_E = cbind(data0, dataE)
  data_E$Timepoint = 'E'
  
  data_F = cbind(data0, dataF)
  data_F$Timepoint = 'F'
  
  org_data = rbind(data_A, data_B, data_C, data_D, data_E, data_F)
  org_data = org_data[!is.na(org_data[, 1]), ]
  # Male 0, female 1
  org_data$Sex = ifelse(org_data$Sex == 'M', 0, 1)
  
  # Yes 1, no 0
  if("Survival (y/n)" %in%  colnames(org_data)){
    org_data$`Survival (y/n)` = ifelse(org_data$`Survival (y/n)` == 'N', 0, 1)
  }
  return(org_data)
}

# Read the data
septic = read_longitudinal("septic.csv")
septic_MFI = read_longitudinal("septic_MFI.csv")
non_septic = read_longitudinal("non_septic.csv")
non_septic_MFI = read_longitudinal("non_septic_MFI.csv")
healthy = read_csv("healthy.csv", 
                   col_types = paste("ffnd", strrep('d', 54), sep = ""),
                   skip = 1)
healthy_MFI = read_csv("healthy_MFI.csv", 
                       col_types = paste("ffnd", strrep('d', 54), sep = ""),
                       skip = 1)
healthy$Sex = ifelse(healthy$Sex == 'M', 0, 1)
healthy_MFI$Sex = ifelse(healthy_MFI$Sex == 'M', 0, 1)

#-----------------------------------------------------------------------------
# Correlation function
correlation = function(var, data, TP){
  var1 = var[1]
  var2 = var[2]
  if("Timepoint" %in%  colnames(data)){
  data = data %>% filter(Timepoint == TP)
  }
  x1 = get(var1, data)
  x2 = get(var2, data)
  if(length(na.omit(x1)) < 4 | length(na.omit(x2)) < 4){
    p_value = "No enough observations"
    cor_coef = "No enough observations"
  }else{
    cor_test = cor.test(x1, x2, na.action = na.omit)
    #cor_coef = cor(x1, x2, use = 'complete.obs')
    p_value = cor_test$p.value
    cor_coef = cor_test$estimate
  }
  return (list(p.value = p_value, cor_coef = cor_coef))
}


correlation_pv = function(var, data, TP, cname){
  var1 = var[[1]]
  var2 = cname[var == 1]
  var1 = rep(var1, length(var2))
  cor_result = apply(cbind(var1, var2), 1, correlation, 
                     data = data, TP = TP)
  p_value = sapply(cor_result, `[[`, 1)
  res = rep(NA, length(var))
  res[var == 1] = p_value
  res[1] = var[[1]]
  return (res)
}

correlation_coef = function(var, data, TP, cname){
  var1 = var[[1]]
  var2 = cname[var == 1]
  var1 = rep(var1, length(var2))
  cor_result = apply(cbind(var1, var2), 1, correlation, 
                     data = data, TP = TP)
  coef = sapply(cor_result, `[[`, 2)
  res = rep(NA, length(var))
  res[var == 1] = coef
  res[1] = var[[1]]
  return (res)
}


#------------------------------------------------------------------------
# Correlation analysis (Table 1)
cor_var1 = read_csv("correlation variables1.csv")

## septic
cor_result = apply(cor_var1, 1, correlation, 
                 data = septic, TP = 'A')
p_value = sapply(cor_result, `[[`, 1)
cor_coef = sapply(cor_result, `[[`, 2)
septic_result1 = cbind(cor_var1[, 1:2], p_value, cor_coef)

## septic_MFI
cor_result = apply(cor_var1, 1, correlation, 
                   data = septic_MFI, TP = 'A')
p_value = sapply(cor_result, `[[`, 1)
cor_coef = sapply(cor_result, `[[`, 2)
septic_MFI_result1 = cbind(cor_var1[, 1:2], p_value, cor_coef)

## non_septic
cor_result = apply(cor_var1, 1, correlation, 
                   data = non_septic, TP = 'A')
p_value = sapply(cor_result, `[[`, 1)
cor_coef = sapply(cor_result, `[[`, 2)
non_septic_result1 = cbind(cor_var1[, 1:2], p_value, cor_coef)

## non_septic_MFI
cor_result = apply(cor_var1, 1, correlation, 
                   data = non_septic_MFI, TP = 'A')
p_value = sapply(cor_result, `[[`, 1)
cor_coef = sapply(cor_result, `[[`, 2)
non_septic_MFI_result1 = cbind(cor_var1[, 1:2], p_value, cor_coef)

## healthy
cor_result = apply(cor_var1[-c(1:4, 7, 11:19), ], 1, correlation, 
                   data = healthy, TP = 'A')
p_value = sapply(cor_result, `[[`, 1)
cor_coef = sapply(cor_result, `[[`, 2)
healthy_result1 = cbind(cor_var1[-c(1:4, 7, 11:19), 1:2], p_value, cor_coef)

## healthy_MFI
cor_result = apply(cor_var1[-c(1:4, 7, 11:19), ], 1, correlation, 
                   data = healthy_MFI, TP = 'A')
p_value = sapply(cor_result, `[[`, 1)
cor_coef = sapply(cor_result, `[[`, 2)
healthy_MFI_result1 = cbind(cor_var1[-c(1:4, 7, 11:19), 1:2 ], p_value, cor_coef)


#------------------------------------------------------------------------
# Correlation analysis (Table 2)
cor_var2 = read_csv("correlation variables2.csv")

# Timepoint A
## septic
septic_A_pv = apply(cor_var2, 1, correlation_pv, 
                   data = septic, TP = 'A', cname = colnames(cor_var2))

septic_A_coef = apply(cor_var2, 1, correlation_coef, 
                    data = septic, TP = 'A', cname = colnames(cor_var2))

septic_A_result = rbind(t(septic_A_pv), rep(' ',40), t(septic_A_coef))
colnames(septic_A_result) = colnames(cor_var2)

## septic_MFI
septic_MFI_A_pv = apply(cor_var2, 1, correlation_pv, 
                    data = septic_MFI, TP = 'A', cname = colnames(cor_var2))

septic_MFI_A_coef = apply(cor_var2, 1, correlation_coef, 
                      data = septic_MFI, TP = 'A', cname = colnames(cor_var2))
septic_MFI_A_result = rbind(t(septic_MFI_A_pv), rep(' ',40), t(septic_MFI_A_coef))
colnames(septic_MFI_A_result) = colnames(cor_var2)

## non_septic
non_septic_A_pv = apply(cor_var2, 1, correlation_pv, 
                    data = non_septic, TP = 'A', cname = colnames(cor_var2))

non_septic_A_coef = apply(cor_var2, 1, correlation_coef, 
                      data = non_septic, TP = 'A', cname = colnames(cor_var2))
non_septic_A_result = rbind(t(non_septic_A_pv), rep(' ',40), t(non_septic_A_coef))
colnames(non_septic_A_result) = colnames(cor_var2)

## non_septic_MFI
non_septic_MFI_A_pv = apply(cor_var2, 1, correlation_pv, 
                    data = non_septic_MFI, TP = 'A', cname = colnames(cor_var2))

non_septic_MFI_A_coef = apply(cor_var2, 1, correlation_coef, 
                      data = non_septic_MFI, TP = 'A', cname = colnames(cor_var2))
non_septic_MFI_A_result = rbind(t(non_septic_MFI_A_pv), rep(' ',40), t(non_septic_MFI_A_coef))
colnames(non_septic_MFI_A_result) = colnames(cor_var2)

## healthy
healthy_pv = apply(cor_var2[-6, ], 1, correlation_pv, 
                    data = healthy, TP = 'A', cname = colnames(cor_var2))

healthy_coef = apply(cor_var2[-6, ], 1, correlation_coef, 
                      data = healthy, TP = 'A', cname = colnames(cor_var2))
healthy_result = rbind(t(healthy_pv), rep(' ',40), t(healthy_coef))
colnames(healthy_result) = colnames(cor_var2)

## healthy_MFI
healthy_MFI_pv = apply(cor_var2[-6, ], 1, correlation_pv, 
                    data = healthy_MFI, TP = 'A', cname = colnames(cor_var2))

healthy_MFI_coef = apply(cor_var2[-6, ], 1, correlation_coef, 
                      data = healthy_MFI, TP = 'A', cname = colnames(cor_var2))

healthy_MFI_result = rbind(t(healthy_MFI_pv), rep(' ',40), t(healthy_MFI_coef))
colnames(healthy_MFI_result) = colnames(cor_var2)


# Timepoint B
## septic
septic_B_pv = apply(cor_var2, 1, correlation_pv, 
                    data = septic, TP = 'B', cname = colnames(cor_var2))

septic_B_coef = apply(cor_var2, 1, correlation_coef, 
                      data = septic, TP = 'B', cname = colnames(cor_var2))

septic_B_result = rbind(t(septic_B_pv), rep(' ',40), t(septic_B_coef))
colnames(septic_B_result) = colnames(cor_var2)

## septic_MFI
septic_MFI_B_pv = apply(cor_var2, 1, correlation_pv, 
                        data = septic_MFI, TP = 'B', cname = colnames(cor_var2))

septic_MFI_B_coef = apply(cor_var2, 1, correlation_coef, 
                          data = septic_MFI, TP = 'B', cname = colnames(cor_var2))
septic_MFI_B_result = rbind(t(septic_MFI_B_pv), rep(' ',40), t(septic_MFI_B_coef))
colnames(septic_MFI_B_result) = colnames(cor_var2)

## non_septic
non_septic_B_pv = apply(cor_var2, 1, correlation_pv, 
                        data = non_septic, TP = 'B', cname = colnames(cor_var2))

non_septic_B_coef = apply(cor_var2, 1, correlation_coef, 
                          data = non_septic, TP = 'B', cname = colnames(cor_var2))
non_septic_B_result = rbind(t(non_septic_B_pv), rep(' ',40), t(non_septic_B_coef))
colnames(non_septic_B_result) = colnames(cor_var2)

## non_septic_MFI
non_septic_MFI_B_pv = apply(cor_var2, 1, correlation_pv, 
                            data = non_septic_MFI, TP = 'B', cname = colnames(cor_var2))

non_septic_MFI_B_coef = apply(cor_var2, 1, correlation_coef, 
                              data = non_septic_MFI, TP = 'B', cname = colnames(cor_var2))
non_septic_MFI_B_result = rbind(t(non_septic_MFI_B_pv), rep(' ',40), t(non_septic_MFI_B_coef))
colnames(non_septic_MFI_B_result) = colnames(cor_var2)


# Timepoint C
## septic
septic_C_pv = apply(cor_var2, 1, correlation_pv, 
                    data = septic, TP = 'C', cname = colnames(cor_var2))

septic_C_coef = apply(cor_var2, 1, correlation_coef, 
                      data = septic, TP = 'C', cname = colnames(cor_var2))

septic_C_result = rbind(t(septic_C_pv), rep(' ',40), t(septic_C_coef))
colnames(septic_C_result) = colnames(cor_var2)

## septic_MFI
septic_MFI_C_pv = apply(cor_var2, 1, correlation_pv, 
                        data = septic_MFI, TP = 'C', cname = colnames(cor_var2))

septic_MFI_C_coef = apply(cor_var2, 1, correlation_coef, 
                          data = septic_MFI, TP = 'C', cname = colnames(cor_var2))
septic_MFI_C_result = rbind(t(septic_MFI_C_pv), rep(' ',40), t(septic_MFI_C_coef))
colnames(septic_MFI_C_result) = colnames(cor_var2)

## non_septic
non_septic_C_pv = apply(cor_var2, 1, correlation_pv, 
                        data = non_septic, TP = 'C', cname = colnames(cor_var2))

non_septic_C_coef = apply(cor_var2, 1, correlation_coef, 
                          data = non_septic, TP = 'C', cname = colnames(cor_var2))
non_septic_C_result = rbind(t(non_septic_C_pv), rep(' ',40), t(non_septic_C_coef))
colnames(non_septic_C_result) = colnames(cor_var2)

## non_septic_MFI
non_septic_MFI_C_pv = apply(cor_var2, 1, correlation_pv, 
                            data = non_septic_MFI, TP = 'C', cname = colnames(cor_var2))

non_septic_MFI_C_coef = apply(cor_var2, 1, correlation_coef, 
                              data = non_septic_MFI, TP = 'C', cname = colnames(cor_var2))
non_septic_MFI_C_result = rbind(t(non_septic_MFI_C_pv), rep(' ',40), t(non_septic_MFI_C_coef))
colnames(non_septic_MFI_C_result) = colnames(cor_var2)


# Timepoint D
## septic
septic_D_pv = apply(cor_var2, 1, correlation_pv, 
                    data = septic, TP = 'D', cname = colnames(cor_var2))

septic_D_coef = apply(cor_var2, 1, correlation_coef, 
                      data = septic, TP = 'D', cname = colnames(cor_var2))

septic_D_result = rbind(t(septic_D_pv), rep(' ',40), t(septic_D_coef))
colnames(septic_D_result) = colnames(cor_var2)

## septic_MFI
septic_MFI_D_pv = apply(cor_var2, 1, correlation_pv, 
                        data = septic_MFI, TP = 'D', cname = colnames(cor_var2))

septic_MFI_D_coef = apply(cor_var2, 1, correlation_coef, 
                          data = septic_MFI, TP = 'D', cname = colnames(cor_var2))
septic_MFI_D_result = rbind(t(septic_MFI_D_pv), rep(' ',40), t(septic_MFI_D_coef))
colnames(septic_MFI_D_result) = colnames(cor_var2)

## non_septic
non_septic_D_pv = apply(cor_var2, 1, correlation_pv, 
                        data = non_septic, TP = 'D', cname = colnames(cor_var2))

non_septic_D_coef = apply(cor_var2, 1, correlation_coef, 
                          data = non_septic, TP = 'D', cname = colnames(cor_var2))
non_septic_D_result = rbind(t(non_septic_D_pv), rep(' ',40), t(non_septic_D_coef))
colnames(non_septic_D_result) = colnames(cor_var2)

## non_septic_MFI
non_septic_MFI_D_pv = apply(cor_var2, 1, correlation_pv, 
                            data = non_septic_MFI, TP = 'D', cname = colnames(cor_var2))

non_septic_MFI_D_coef = apply(cor_var2, 1, correlation_coef, 
                              data = non_septic_MFI, TP = 'D', cname = colnames(cor_var2))
non_septic_MFI_D_result = rbind(t(non_septic_MFI_D_pv), rep(' ',40), t(non_septic_MFI_D_coef))
colnames(non_septic_MFI_D_result) = colnames(cor_var2)


# Timepoint E
## septic
septic_E_pv = apply(cor_var2, 1, correlation_pv, 
                    data = septic, TP = 'E', cname = colnames(cor_var2))

septic_E_coef = apply(cor_var2, 1, correlation_coef, 
                      data = septic, TP = 'E', cname = colnames(cor_var2))

septic_E_result = rbind(t(septic_E_pv), rep(' ',40), t(septic_E_coef))
colnames(septic_E_result) = colnames(cor_var2)

## septic_MFI
septic_MFI_E_pv = apply(cor_var2, 1, correlation_pv, 
                        data = septic_MFI, TP = 'E', cname = colnames(cor_var2))

septic_MFI_E_coef = apply(cor_var2, 1, correlation_coef, 
                          data = septic_MFI, TP = 'E', cname = colnames(cor_var2))
septic_MFI_E_result = rbind(t(septic_MFI_E_pv), rep(' ',40), t(septic_MFI_E_coef))
colnames(septic_MFI_E_result) = colnames(cor_var2)

## non_septic
non_septic_E_pv = apply(cor_var2, 1, correlation_pv, 
                        data = non_septic, TP = 'E', cname = colnames(cor_var2))

non_septic_E_coef = apply(cor_var2, 1, correlation_coef, 
                          data = non_septic, TP = 'E', cname = colnames(cor_var2))
non_septic_E_result = rbind(t(non_septic_E_pv), rep(' ',40), t(non_septic_E_coef))
colnames(non_septic_E_result) = colnames(cor_var2)

## non_septic_MFI
non_septic_MFI_E_pv = apply(cor_var2, 1, correlation_pv, 
                            data = non_septic_MFI, TP = 'E', cname = colnames(cor_var2))

non_septic_MFI_E_coef = apply(cor_var2, 1, correlation_coef, 
                              data = non_septic_MFI, TP = 'E', cname = colnames(cor_var2))
non_septic_MFI_E_result = rbind(t(non_septic_MFI_E_pv), rep(' ',40), t(non_septic_MFI_E_coef))
colnames(non_septic_MFI_E_result) = colnames(cor_var2)


# Timepoint F
## septic
septic_F_pv = apply(cor_var2, 1, correlation_pv, 
                    data = septic, TP = 'F', cname = colnames(cor_var2))

septic_F_coef = apply(cor_var2, 1, correlation_coef, 
                      data = septic, TP = 'F', cname = colnames(cor_var2))

septic_F_result = rbind(t(septic_F_pv), rep(' ',40), t(septic_F_coef))
colnames(septic_F_result) = colnames(cor_var2)

## septic_MFI
septic_MFI_F_pv = apply(cor_var2, 1, correlation_pv, 
                        data = septic_MFI, TP = 'F', cname = colnames(cor_var2))

septic_MFI_F_coef = apply(cor_var2, 1, correlation_coef, 
                          data = septic_MFI, TP = 'F', cname = colnames(cor_var2))
septic_MFI_F_result = rbind(t(septic_MFI_F_pv), rep(' ',40), t(septic_MFI_F_coef))
colnames(septic_MFI_F_result) = colnames(cor_var2)

## non_septic
non_septic_F_pv = apply(cor_var2, 1, correlation_pv, 
                        data = non_septic, TP = 'F', cname = colnames(cor_var2))

non_septic_F_coef = apply(cor_var2, 1, correlation_coef, 
                          data = non_septic, TP = 'F', cname = colnames(cor_var2))
non_septic_F_result = rbind(t(non_septic_F_pv), rep(' ',40), t(non_septic_F_coef))
colnames(non_septic_F_result) = colnames(cor_var2)

## non_septic_MFI
non_septic_MFI_F_pv = apply(cor_var2, 1, correlation_pv, 
                            data = non_septic_MFI, TP = 'F', cname = colnames(cor_var2))

non_septic_MFI_F_coef = apply(cor_var2, 1, correlation_coef, 
                              data = non_septic_MFI, TP = 'F', cname = colnames(cor_var2))
non_septic_MFI_F_result = rbind(t(non_septic_MFI_F_pv), rep(' ',40), t(non_septic_MFI_F_coef))
colnames(non_septic_MFI_F_result) = colnames(cor_var2)

#------------------------------------------------------------------------
# Multivariate analysis
## The function for logistic regression model
multi_sur = function(data, TP, cname){
  if("Timepoint" %in%  colnames(data)){
    data = data %>% filter(Timepoint == TP)
  }
  aa = which(colnames(data) %in% c(cname, 'Survival (y/n)'))
  data = na.omit(data[ ,aa])
  if(nrow(data) < 3) return ('Not enough observations')
  logit_model = glm(`Survival (y/n)` ~ ., data = data, family = 'binomial')
  null_model = glm(`Survival (y/n)` ~ 1, data = data, family = 'binomial')
  bb = anova(logit_model, null_model, test = 'Chisq')
  return (bb$`Pr(>Chi)`[2])
}

cname = colnames(cor_var2)[-(1:6)]
Timepoints = c('A', 'B', 'C', 'D', 'E', 'F')

## Septic
septic_logit = sapply(Timepoints, multi_sur, data = septic, cname = cname)

## Septic_MFI
septic_MFI_logit = sapply(Timepoints, multi_sur, data = septic_MFI, cname = cname)

## Non_septic
non_septic_logit = sapply(Timepoints, multi_sur, data = non_septic, cname = cname)

## Non_septic_MFI
non_septic_MFI_logit = sapply(Timepoints, multi_sur, data = non_septic_MFI, cname = cname)

multi_res = data.frame(septic = septic_logit, septic_MFI = septic_MFI_logit,
                       non_septic = non_septic_logit, non_septic_MFI = non_septic_MFI_logit)

library(xlsx)
write.xlsx(multi_res, file="Multivariate.xlsx", row.names=T)

#------------------------------------------------------------------------
# Export results
## Table 1
library(xlsx)
write.xlsx(septic_result1, file="Correltion1.xlsx", sheetName="septic(%)", row.names=FALSE)
write.xlsx(septic_MFI_result1, file="Correltion1.xlsx", sheetName="septic(MFI)", append=TRUE, row.names=FALSE)
write.xlsx(non_septic_result1, file="Correltion1.xlsx", sheetName="non_septic(%)", append=TRUE, row.names=FALSE)
write.xlsx(non_septic_MFI_result1, file="Correltion1.xlsx", sheetName="non_septic(MFI)", append=TRUE, row.names=FALSE)
write.xlsx(healthy_result1, file="Correltion1.xlsx", sheetName="healthy(%)", append=TRUE, row.names=FALSE)
write.xlsx(healthy_MFI_result1, file="Correltion1.xlsx", sheetName="healthy(MFI)", append=TRUE, row.names=FALSE)

## Table 2 (A)
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "septic(%)")
addWorksheet(wb, "septic(MFI)")
addWorksheet(wb, "non_septic(%)")
addWorksheet(wb, "non_septic(MFI)")
addWorksheet(wb, "healthy(%)")
addWorksheet(wb, "healthy(MFI)")

negStyle <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
## rule applies to all each cell in range
aa <- septic_A_result[,1]
septic_A_result = septic_A_result[, -1]
septic_A_result[is.na(septic_A_result)] = 10000
septic_A_result = apply(septic_A_result, 2, as.numeric)
rownames(septic_A_result) = aa

aa <- septic_MFI_A_result[,1]
septic_MFI_A_result = septic_MFI_A_result[, -1]
septic_MFI_A_result[is.na(septic_MFI_A_result)] = 10000
septic_MFI_A_result = apply(septic_MFI_A_result, 2, as.numeric)
rownames(septic_MFI_A_result) = aa

aa <- non_septic_A_result[,1]
non_septic_A_result = non_septic_A_result[, -1]
non_septic_A_result[is.na(non_septic_A_result)] = 10000
non_septic_A_result = apply(non_septic_A_result, 2, as.numeric)
rownames(non_septic_A_result) = aa

aa <- non_septic_MFI_A_result[,1]
non_septic_MFI_A_result = non_septic_MFI_A_result[, -1]
non_septic_MFI_A_result[is.na(non_septic_MFI_A_result)] = 10000
non_septic_MFI_A_result = apply(non_septic_MFI_A_result, 2, as.numeric)
rownames(non_septic_MFI_A_result) = aa

aa <- healthy_result[,1]
healthy_result = healthy_result[, -1]
healthy_result[is.na(healthy_result)] = 10000
healthy_result = apply(healthy_result, 2, as.numeric)
rownames(healthy_result) = aa

aa <- healthy_MFI_result[,1]
healthy_MFI_result = healthy_MFI_result[, -1]
healthy_MFI_result[is.na(healthy_MFI_result)] = 10000
healthy_MFI_result = apply(healthy_MFI_result, 2, as.numeric)
rownames(healthy_MFI_result) = aa

writeData(wb, "septic(%)", septic_A_result, rowNames = T)
writeData(wb, "septic(MFI)", septic_MFI_A_result, rowNames = T)
writeData(wb, "non_septic(%)", non_septic_A_result, rowNames = T)
writeData(wb, "non_septic(MFI)", non_septic_MFI_A_result, rowNames = T)
writeData(wb, "healthy(%)", healthy_result, rowNames = T)
writeData(wb, "healthy(MFI)", healthy_MFI_result, rowNames = T)

conditionalFormatting(wb, "septic(%)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)
conditionalFormatting(wb, "septic(MFI)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)
conditionalFormatting(wb, "non_septic(%)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)
conditionalFormatting(wb, "non_septic(MFI)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)
conditionalFormatting(wb, "healthy(%)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)
conditionalFormatting(wb, "healthy(MFI)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)

saveWorkbook(wb, "Timepoint A.xlsx", TRUE)

## Table 2 (B)
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "septic(%)")
addWorksheet(wb, "septic(MFI)")
addWorksheet(wb, "non_septic(%)")
addWorksheet(wb, "non_septic(MFI)")

negStyle <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
## rule applies to all each cell in range
aa <- septic_B_result[,1]
septic_B_result = septic_B_result[, -1]
septic_B_result[is.na(septic_B_result)] = 10000
septic_B_result = apply(septic_B_result, 2, as.numeric)
rownames(septic_B_result) = aa

aa <- septic_MFI_B_result[,1]
septic_MFI_B_result = septic_MFI_B_result[, -1]
septic_MFI_B_result[is.na(septic_MFI_B_result)] = 10000
septic_MFI_B_result = apply(septic_MFI_B_result, 2, as.numeric)
rownames(septic_MFI_B_result) = aa

aa <- non_septic_B_result[,1]
non_septic_B_result = non_septic_B_result[, -1]
non_septic_B_result[is.na(non_septic_B_result)] = 10000
non_septic_B_result = apply(non_septic_B_result, 2, as.numeric)
rownames(non_septic_B_result) = aa

aa <- non_septic_MFI_B_result[,1]
non_septic_MFI_B_result = non_septic_MFI_B_result[, -1]
non_septic_MFI_B_result[is.na(non_septic_MFI_B_result)] = 10000
non_septic_MFI_B_result = apply(non_septic_MFI_B_result, 2, as.numeric)
rownames(non_septic_MFI_B_result) = aa


writeData(wb, "septic(%)", septic_B_result, rowNames = T)
writeData(wb, "septic(MFI)", septic_MFI_B_result, rowNames = T)
writeData(wb, "non_septic(%)", non_septic_B_result, rowNames = T)
writeData(wb, "non_septic(MFI)", non_septic_MFI_B_result, rowNames = T)

conditionalFormatting(wb, "septic(%)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)
conditionalFormatting(wb, "septic(MFI)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)
conditionalFormatting(wb, "non_septic(%)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)
conditionalFormatting(wb, "non_septic(MFI)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)

saveWorkbook(wb, "Timepoint B.xlsx", TRUE)


## Table 2 (C)
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "septic(%)")
addWorksheet(wb, "septic(MFI)")
addWorksheet(wb, "non_septic(%)")
addWorksheet(wb, "non_septic(MFI)")

negStyle <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
## rule applies to all each cell in range
aa <- septic_C_result[,1]
septic_C_result = septic_C_result[, -1]
septic_C_result[is.na(septic_C_result)] = 10000
septic_C_result = apply(septic_C_result, 2, as.numeric)
rownames(septic_C_result) = aa

aa <- septic_MFI_C_result[,1]
septic_MFI_C_result = septic_MFI_C_result[, -1]
septic_MFI_C_result[is.na(septic_MFI_C_result)] = 10000
septic_MFI_C_result = apply(septic_MFI_C_result, 2, as.numeric)
rownames(septic_MFI_C_result) = aa

aa <- non_septic_C_result[,1]
non_septic_C_result = non_septic_C_result[, -1]
non_septic_C_result[is.na(non_septic_C_result)] = 10000
non_septic_C_result = apply(non_septic_C_result, 2, as.numeric)
rownames(non_septic_C_result) = aa

aa <- non_septic_MFI_C_result[,1]
non_septic_MFI_C_result = non_septic_MFI_C_result[, -1]
non_septic_MFI_C_result[is.na(non_septic_MFI_C_result)] = 10000
non_septic_MFI_C_result = apply(non_septic_MFI_C_result, 2, as.numeric)
rownames(non_septic_MFI_C_result) = aa


writeData(wb, "septic(%)", septic_C_result, rowNames = T)
writeData(wb, "septic(MFI)", septic_MFI_C_result, rowNames = T)
writeData(wb, "non_septic(%)", non_septic_C_result, rowNames = T)
writeData(wb, "non_septic(MFI)", non_septic_MFI_C_result, rowNames = T)

conditionalFormatting(wb, "septic(%)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)
conditionalFormatting(wb, "septic(MFI)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)
conditionalFormatting(wb, "non_septic(%)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)
conditionalFormatting(wb, "non_septic(MFI)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)

saveWorkbook(wb, "Timepoint C.xlsx", TRUE)


## Table 2 (D)
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "septic(%)")
addWorksheet(wb, "septic(MFI)")
addWorksheet(wb, "non_septic(%)")
addWorksheet(wb, "non_septic(MFI)")

negStyle <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
## rule applies to all each cell in range
aa <- septic_D_result[,1]
septic_D_result = septic_D_result[, -1]
septic_D_result[is.na(septic_D_result)] = 10000
septic_D_result = apply(septic_D_result, 2, as.numeric)
rownames(septic_D_result) = aa

aa <- septic_MFI_D_result[,1]
septic_MFI_D_result = septic_MFI_D_result[, -1]
septic_MFI_D_result[is.na(septic_MFI_D_result)] = 10000
septic_MFI_D_result = apply(septic_MFI_D_result, 2, as.numeric)
rownames(septic_MFI_D_result) = aa

aa <- non_septic_D_result[,1]
non_septic_D_result = non_septic_D_result[, -1]
non_septic_D_result[is.na(non_septic_D_result)] = 10000
non_septic_D_result = apply(non_septic_D_result, 2, as.numeric)
rownames(non_septic_D_result) = aa

aa <- non_septic_MFI_D_result[,1]
non_septic_MFI_D_result = non_septic_MFI_D_result[, -1]
non_septic_MFI_D_result[is.na(non_septic_MFI_D_result)] = 10000
non_septic_MFI_D_result = apply(non_septic_MFI_D_result, 2, as.numeric)
rownames(non_septic_MFI_D_result) = aa


writeData(wb, "septic(%)", septic_D_result, rowNames = T)
writeData(wb, "septic(MFI)", septic_MFI_D_result, rowNames = T)
writeData(wb, "non_septic(%)", non_septic_D_result, rowNames = T)
writeData(wb, "non_septic(MFI)", non_septic_MFI_D_result, rowNames = T)

conditionalFormatting(wb, "septic(%)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)
conditionalFormatting(wb, "septic(MFI)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)
conditionalFormatting(wb, "non_septic(%)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)
conditionalFormatting(wb, "non_septic(MFI)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)

saveWorkbook(wb, "Timepoint D.xlsx", TRUE)


## Table 2 (E)
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "septic(%)")
addWorksheet(wb, "septic(MFI)")
addWorksheet(wb, "non_septic(%)")
addWorksheet(wb, "non_septic(MFI)")

negStyle <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
## rule applies to all each cell in range
aa <- septic_E_result[,1]
septic_E_result = septic_E_result[, -1]
septic_E_result[is.na(septic_E_result)] = 10000
septic_E_result = apply(septic_E_result, 2, as.numeric)
rownames(septic_E_result) = aa

aa <- septic_MFI_E_result[,1]
septic_MFI_E_result = septic_MFI_E_result[, -1]
septic_MFI_E_result[is.na(septic_MFI_E_result)] = 10000
septic_MFI_E_result = apply(septic_MFI_E_result, 2, as.numeric)
rownames(septic_MFI_E_result) = aa

aa <- non_septic_E_result[,1]
non_septic_E_result = non_septic_E_result[, -1]
non_septic_E_result[is.na(non_septic_E_result)] = 10000
non_septic_E_result = apply(non_septic_E_result, 2, as.numeric)
rownames(non_septic_E_result) = aa

aa <- non_septic_MFI_E_result[,1]
non_septic_MFI_E_result = non_septic_MFI_E_result[, -1]
non_septic_MFI_E_result[is.na(non_septic_MFI_E_result)] = 10000
non_septic_MFI_E_result = apply(non_septic_MFI_E_result, 2, as.numeric)
rownames(non_septic_MFI_E_result) = aa


writeData(wb, "septic(%)", septic_E_result, rowNames = T)
writeData(wb, "septic(MFI)", septic_MFI_E_result, rowNames = T)
writeData(wb, "non_septic(%)", non_septic_E_result, rowNames = T)
writeData(wb, "non_septic(MFI)", non_septic_MFI_E_result, rowNames = T)

conditionalFormatting(wb, "septic(%)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)
conditionalFormatting(wb, "septic(MFI)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)
conditionalFormatting(wb, "non_septic(%)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)
conditionalFormatting(wb, "non_septic(MFI)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)

saveWorkbook(wb, "Timepoint E.xlsx", TRUE)

## Table 2 (F)
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "septic(%)")
addWorksheet(wb, "septic(MFI)")
addWorksheet(wb, "non_septic(%)")
addWorksheet(wb, "non_septic(MFI)")

negStyle <- createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
## rule applies to all each cell in range
aa <- septic_F_result[,1]
septic_F_result = septic_F_result[, -1]
septic_F_result[is.na(septic_F_result)] = 10000
septic_F_result = apply(septic_F_result, 2, as.numeric)
rownames(septic_F_result) = aa

aa <- septic_MFI_F_result[,1]
septic_MFI_F_result = septic_MFI_F_result[, -1]
septic_MFI_F_result[is.na(septic_MFI_F_result)] = 10000
septic_MFI_F_result = apply(septic_MFI_F_result, 2, as.numeric)
rownames(septic_MFI_F_result) = aa

aa <- non_septic_F_result[,1]
non_septic_F_result = non_septic_F_result[, -1]
non_septic_F_result[is.na(non_septic_F_result)] = 10000
non_septic_F_result = apply(non_septic_F_result, 2, as.numeric)
rownames(non_septic_F_result) = aa

aa <- non_septic_MFI_F_result[,1]
non_septic_MFI_F_result = non_septic_MFI_F_result[, -1]
non_septic_MFI_F_result[is.na(non_septic_MFI_F_result)] = 10000
non_septic_MFI_F_result = apply(non_septic_MFI_F_result, 2, as.numeric)
rownames(non_septic_MFI_F_result) = aa


writeData(wb, "septic(%)", septic_F_result, rowNames = T)
writeData(wb, "septic(MFI)", septic_MFI_F_result, rowNames = T)
writeData(wb, "non_septic(%)", non_septic_F_result, rowNames = T)
writeData(wb, "non_septic(MFI)", non_septic_MFI_F_result, rowNames = T)

conditionalFormatting(wb, "septic(%)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)
conditionalFormatting(wb, "septic(MFI)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)
conditionalFormatting(wb, "non_septic(%)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)
conditionalFormatting(wb, "non_septic(MFI)", cols=2:41, rows=2:9, rule="<0.05", style = negStyle)

saveWorkbook(wb, "Timepoint F.xlsx", TRUE)
#----------------------------------------------------------------------
# Data Visualization
## First group of correlation plots
library(ggpubr)
Mcorp1 = ggplot(septic[septic$Timepoint == 'A', ], aes(x = `M-cd69`, y = `M-gzB pae`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

Mcorp2 = ggplot(septic[septic$Timepoint == 'A', ], aes(x = `M-cd69`, y = `M-IFNy pae`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

Mcorp3 = ggplot(septic[septic$Timepoint == 'A', ], aes(x = `M-cd69`, y = `M-IFNy sau`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

Mcorp4 = ggplot(septic[septic$Timepoint == 'A', ], aes(x = `M-cd69`, y = `M-IFNy sty`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

Mcorp5 = ggplot(septic[septic$Timepoint == 'A', ], aes(x = `M-cd69`, y = `M-IL-2 sau`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

Mcorp6 = ggplot(septic[septic$Timepoint == 'A', ], aes(x = `M-cd69`, y = `M-IL-2 sty`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

library(cowplot)
plot_grid(Mcorp1, Mcorp2, Mcorp3, Mcorp4, Mcorp5, Mcorp6, labels = "AUTO")

## Second group of correlation plots
CDcorp1 = ggplot(septic_MFI[septic_MFI$Timepoint == 'B', ], aes(x = `CD14 HLADR`, y = `M-gzB media`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

CDcorp2 = ggplot(septic_MFI[septic_MFI$Timepoint == 'B', ], aes(x = `CD14 HLADR`, y = `M-gzB 12+18`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

CDcorp3 = ggplot(septic_MFI[septic_MFI$Timepoint == 'B', ], aes(x = `CD14 HLADR`, y = `M-IFNy media`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

CDcorp4 = ggplot(septic_MFI[septic_MFI$Timepoint == 'B', ], aes(x = `CD14 HLADR`, y = `M-IFNy 12+cd28+ligand`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

CDcorp5 = ggplot(septic_MFI[septic_MFI$Timepoint == 'B', ], aes(x = `CD14 HLADR`, y = `M-IL-10 Ecoli`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

CDcorp6 = ggplot(septic_MFI[septic_MFI$Timepoint == 'B', ], aes(x = `CD14 HLADR`, y = `M-IL-10 kpn`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

CDcorp7 = ggplot(septic_MFI[septic_MFI$Timepoint == 'B', ], aes(x = `CD14 HLADR`, y = `M-IL-2 12+cd28+ligand`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

plot_grid(CDcorp1, CDcorp2, CDcorp3, CDcorp4, CDcorp5, CDcorp6, CDcorp7)


## Second plot without outliers
septic_MFI_B <- septic_MFI[septic_MFI$Timepoint == 'B', ]
CDcorp4_out = ggplot(septic_MFI_B[-16, ], aes(x = `CD14 HLADR`, y = `M-IFNy 12+cd28+ligand`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

CDcorp7_out = ggplot(septic_MFI_B[-16, ], aes(x = `CD14 HLADR`, y = `M-IL-2 12+cd28+ligand`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

plot_grid(CDcorp4_out, CDcorp7_out)



## Third group of correlation plots
CDcorpF1 = ggplot(septic_MFI[septic_MFI$Timepoint == 'F', ], aes(x = `CD14 HLADR`, y = `M-lag-3`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

CDcorpF2 = ggplot(septic_MFI[septic_MFI$Timepoint == 'F', ], aes(x = `CD14 HLADR`, y = `M-gzB Ecoli`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

CDcorpF3 = ggplot(septic_MFI[septic_MFI$Timepoint == 'F', ], aes(x = `CD14 HLADR`, y = `M-IFNy 12+cd28+ligand`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

CDcorpF4 = ggplot(septic_MFI[septic_MFI$Timepoint == 'F', ], aes(x = `CD14 HLADR`, y = `M-IL-2 12+cd28+ligand`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

CDcorpF5 = ggplot(septic_MFI[septic_MFI$Timepoint == 'F', ], aes(x = `CD14 HLADR`, y = `M-TNFa kpn`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

CDcorpF6 = ggplot(septic_MFI[septic_MFI$Timepoint == 'F', ], aes(x = `CD14 HLADR`, y = `M-TNFa pae`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

CDcorpF7 = ggplot(septic_MFI[septic_MFI$Timepoint == 'F', ], aes(x = `CD14 HLADR`, y = `M-TNFa  12+cd28+ligand`)) + 
  geom_point() + geom_smooth(method='lm', se = F) +
  stat_cor(method = "pearson")

plot_grid(CDcorpF1, CDcorpF2, CDcorpF3, CDcorpF4, CDcorpF5, CDcorpF6, CDcorpF7)

