setwd("C:/Users/ags10001111/OneDrive/2020 Summer/Epilepsy study")
library(tidyverse)
load('Epilepsy study_2013.RData')

#--------------------------------------
# Function for comparison
compare_fun <- function(data){
  # data <- as.data.frame(demo_data)
  data <- as.data.frame(data)
  SUDEP_data <- data[data$`Cause of death:` == 'SUDEP', ]
  Non_SUDEP_data <- data[data$`Cause of death:` == 'Non-SUDEP', ]
  test_result <- compare_test(SUDEP_data, Non_SUDEP_data)
  return (test_result)
}

compare_test <- function(data1, data2){
  # data1 <- SUDEP_data
  # data2 <- Non_SUDEP_data
  data <- rbind(data1, data2)
  if (ncol(data1) != ncol(data2)){
    stop ('The number of variables in two dataframes are different')
  }
  n <- ncol(data1)
  com_summary <- data.frame(Variables = NA, `Non_SUDEP` = NA, SUDEP = NA, `p-value` = NA)
  # i = 5
  for (i in 1:(n - 1)){
    # Use different tests for different types of variables
    if (class(data1[ , i]) == 'character'){
      tbl <- table(data[, i], droplevels(as.factor(data$`Cause of death:`)))
      if(nrow(tbl) == 1){ 
        p.value <- 10000
        }else{
      # Use Fisher exact test because the sample size is less than 1000
      fisher_result <- fisher.test(tbl, workspace = 2e8)
      p.value <- fisher_result$p.value
      var_name0 <- paste(colnames(data1)[i], ' ,', ' n (%)', sep = '')
      var_name <- c(var_name0, NA , NA, round(p.value, 3))
      com_summary <- rbind(com_summary, var_name)
      
      # Fisher exact test for each value
      require(caret)
      data_i <- as.data.frame(data[,i])
      colnames(data_i) <- colnames(data)[i]
      dmy <- dummyVars(" ~ .", data = data_i)
      trsf <- data.frame(predict(dmy, newdata = data_i))
      num_level <- length(levels(as.factor(data[,i])))
      p.value_ev <- rep(NA, num_level)
      for (j in 1:num_level){
        tbl_ev <- table(trsf[, j], data$`Cause of death:`) 
        fisher_result_ev <- fisher.test(tbl_ev, workspace = 2e8)
        p.value_ev[j] <- fisher_result_ev$p.value
      }
        }
      prop_tbl <- prop.table(tbl, 2)
      prop_tbl <- round(100*prop_tbl, 1)
      com_tbl <- paste(tbl, ' (', prop_tbl, ')', sep = '')
      com_tbl <- matrix(com_tbl, ncol = 2)
      tbl1 <- cbind(rownames(tbl), com_tbl, round(p.value_ev, 3))
      colnames(tbl1) <- colnames(com_summary)
      com_summary <- rbind(com_summary, tbl1)
      com_summary <- rbind(com_summary, NA)
    }
    if (class(data1[, i]) == 'numeric'){
      # Compute p-value using normal approximation with continuity correction
      MW_result <- wilcox.test(data1[, i], data2[, i])
      p.value <- MW_result$p.value
      var_name0 <- paste(colnames(data1)[i], ' ,', 'mean', sep = '')
      var_name <- c(var_name0, NA , NA, round(p.value, 3))
      com_summary <- rbind(com_summary, var_name)
      vec_mean <- c(NA, round(mean(data1[, i], na.rm = T), 2), round(mean(data2[, i], na.rm = T), 2), NA)
      vec_mean <- as.data.frame(t(vec_mean))
      colnames(vec_mean) <- colnames(com_summary)
      com_summary <- rbind(com_summary, vec_mean)
      com_summary <- rbind(com_summary, NA)
    }
  }
  return (com_summary)
}



demo_summary <- compare_fun(demo_data)
#demo_summary <- compare_fun(demo_data[demo_data$`Sex:` == 'Male', c(9,8)])
#write_csv(demo_summary, 'male_age1.csv')
clinical.info_summary <- compare_fun(clinical.info_data)
cod_summary <- compare_fun(cod_data)
neuro_summary <- compare_fun(neuro_data)
AP_summary <- compare_fun(AP_data)

# Export the result
demo_summary$p.value <- as.numeric(demo_summary$p.value)
clinical.info_summary$p.value <- as.numeric(clinical.info_summary$p.value)
cod_summary$p.value <- as.numeric(cod_summary$p.value)
neuro_summary$p.value <- as.numeric(neuro_summary$p.value)
AP_summary$p.value <- as.numeric(AP_summary$p.value)

library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "Patient demographics")
addWorksheet(wb, "Clinical information")
addWorksheet(wb, "Circumstances of death")
addWorksheet(wb, "Neuropathology findings")
addWorksheet(wb, "AP findings")

negStyle <- createStyle(fontColour = "#9C0006")
flagged_style <- createStyle(fgFill = '#c4e68a')

## rule applies to all each cell in range
writeData(wb, "Patient demographics", demo_summary[-1, ], rowNames = F)
writeData(wb, "Clinical information", clinical.info_summary[-1, ], rowNames = F)
writeData(wb, "Circumstances of death", cod_summary[-1, ], rowNames = F)
writeData(wb, "Neuropathology findings", neuro_summary[-1, ], rowNames = F)
writeData(wb, "AP findings", AP_summary[-1, ], rowNames = F)


conditionalFormatting(wb, "Patient demographics", cols=4, rows=1:nrow(demo_summary), rule="<0.05", style = negStyle)
conditionalFormatting(wb, "Clinical information", cols=4, rows=1:nrow(clinical.info_summary), rule="<0.05", style = negStyle)
conditionalFormatting(wb, "Circumstances of death", cols=4, rows=1:nrow(cod_summary), rule="<0.05", style = negStyle)
conditionalFormatting(wb, "Neuropathology findings", cols=4, rows=1:nrow(neuro_summary), rule="<0.05", style = negStyle)
conditionalFormatting(wb, "AP findings", cols=4, rows=1:nrow(AP_summary), rule="<0.05", style = negStyle)

addStyle(wb, "Patient demographics", style = flagged_style, cols=1:4, rows= which(is.na(demo_summary[, 4])) + 1, gridExpand = T, stack = T)
addStyle(wb, "Clinical information", style = flagged_style, cols=1:4, rows= which(is.na(clinical.info_summary[, 4])) + 1, gridExpand = T, stack = T)
addStyle(wb, "Circumstances of death", style = flagged_style, cols=1:4, rows= which(is.na(cod_summary[, 4]))+ 1, gridExpand = T, stack = T)
addStyle(wb, "Neuropathology findings", style = flagged_style, cols=1:4, rows= which(is.na(neuro_summary[, 4])) + 1, gridExpand = T, stack = T)
addStyle(wb, "AP findings", style = flagged_style, cols=1:4, rows= which(is.na(AP_summary[, 4])) + 1, gridExpand = T, stack = T)

saveWorkbook(wb, "Summary_2013.xlsx", TRUE)

