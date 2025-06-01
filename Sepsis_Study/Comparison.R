setwd("C:/Users/ags10001111/OneDrive/2020 Summer/Master Project/Comparison")
load("Comparison.RData")
library(tidyverse)

#-----------------------------------------------------------------------------
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

#--------------------------------------------------------------------------
# EXTRACT DATA
# Extract the data according to the timepoint
## Timepoint A
septic_A = septic %>% filter(Timepoint == "A")
non_septic_A = non_septic %>% filter(Timepoint == "A")
septic_MFI_A = septic_MFI %>% filter(Timepoint == "A")
non_septic_MFI_A = non_septic_MFI %>% filter(Timepoint == "A")

## Timepoint B
septic_B = septic %>% filter(Timepoint == "B")
non_septic_B = non_septic %>% filter(Timepoint == "B")
septic_MFI_B = septic_MFI %>% filter(Timepoint == "B")
non_septic_MFI_B = non_septic_MFI %>% filter(Timepoint == "B")

## Timepoint C
septic_C = septic %>% filter(Timepoint == "C")
non_septic_C = non_septic %>% filter(Timepoint == "C")
septic_MFI_C = septic_MFI %>% filter(Timepoint == "C")
non_septic_MFI_C = non_septic_MFI %>% filter(Timepoint == "C")

## Timepoint D
septic_D = septic %>% filter(Timepoint == "D")
non_septic_D = non_septic %>% filter(Timepoint == "D")
septic_MFI_D = septic_MFI %>% filter(Timepoint == "D")
non_septic_MFI_D = non_septic_MFI %>% filter(Timepoint == "D")

## Timepoint E
septic_E = septic %>% filter(Timepoint == "E")
non_septic_E = non_septic %>% filter(Timepoint == "E")
septic_MFI_E = septic_MFI %>% filter(Timepoint == "E")
non_septic_MFI_E = non_septic_MFI %>% filter(Timepoint == "E")

## Timepoint F
septic_F = septic %>% filter(Timepoint == "F")
non_septic_F = non_septic %>% filter(Timepoint == "F")
septic_MFI_F = septic_MFI %>% filter(Timepoint == "F")
non_septic_MFI_F = non_septic_MFI %>% filter(Timepoint == "F")

#-----------------------------------------------------------------------
# COMPARISON FUNCTION
## Box-cox transformation function
powerTransform = function(y, lam1) {
    if (lam1 == 0) {
      return (log(y))
    } else {
      return ((y^lam1 - 1) / lam1)
    }
}

## t-test function
library(rcompanion)
compare_t = function(var_tc, data1, data2){
    x1 = get(var_tc, data1)
    x2 = get(var_tc, data2)
    p_value = NA
    i0 = NA
    if (length(na.omit(x1)) >= 3 && length(na.omit(x2)) >= 3){
      if(length(unique(na.omit(x1))) >= 2 && length(unique(na.omit(x2))) >= 2){
         if(shapiro.test(x1)$p.value >0.05 && shapiro.test(x2)$p.value > 0.05){
             p_value = t.test(x1, x2)$p.value
         }else{
            for (i in seq(-5, 5, by = 0.25)){
                 x1_norm = powerTransform(x1 + 1, i)
                 x2_norm = powerTransform(x2 + 1, i)
                 if(length(na.omit(x1)) == length(na.omit(x1_norm)) && length(na.omit(x2)) == length(na.omit(x2_norm))){
                    if(length(unique(na.omit(x1_norm))) > 1 && length(unique(na.omit(x2_norm))) > 1){
                       if(shapiro.test(x1_norm)$p.value >0.05 && shapiro.test(x2_norm)$p.value > 0.05){
                          i0 = i
                          break
                       }
                    }
                  }
              }
              p_value = t.test(x1_norm, x2_norm)$p.value
      }
      }else {
      p_value = "Not enough distinct observations"
    }
      }else{
      p_value = "Not enough observations"
    }
    return (list(p.value = p_value, lambda = i0, difference = mean(x1, na.rm = T) - mean(x2, na.rm = T)))
}

## Mann-Whitney U test
compare_MW = function(var_tc, data1, data2){
  x1 = get(var_tc, data1)
  x2 = get(var_tc, data2)
  if (length(na.omit(x1)) >= 3 && length(na.omit(x2)) >= 3){
  p_value = wilcox.test(x1, x2, exact = F)$p.value
  } else {
    p_value = "Not enough observations"
    }
  return (p_value)
}

variable_to_compare = read_csv("Comparison_variable.csv", col_names = F)

#-------------------------------------------------------------------------------
# SEPTIC VS NON_SEPTIC(%)
## Timepoint A
### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                        data1 = septic_A, data2 = non_septic_A)
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

sep_non_A_t = cbind(variable_to_compare, p_value, lambda, difference)

### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = septic_A, data2 =non_septic_A)

sep_non_A_MW = cbind(variable_to_compare, p_value)

## Timepoint B
### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = septic_B, data2 = non_septic_B)
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

sep_non_B_t = cbind(variable_to_compare, p_value, lambda, difference)

### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = septic_B, data2 =non_septic_B)

sep_non_B_MW = cbind(variable_to_compare, p_value)

## Timepoint C
### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = septic_C, data2 = non_septic_C)
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

sep_non_C_t = cbind(variable_to_compare, p_value, lambda, difference)

### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = septic_C, data2 =non_septic_C)

sep_non_C_MW = cbind(variable_to_compare, p_value)


## Timepoint D
### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = septic_D, data2 = non_septic_D)
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

sep_non_D_t = cbind(variable_to_compare, p_value, lambda, difference)

### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = septic_D, data2 =non_septic_D)

sep_non_D_MW = cbind(variable_to_compare, p_value)


## Timepoint E
### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = septic_E, data2 = non_septic_E)
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

sep_non_E_t = cbind(variable_to_compare, p_value, lambda, difference)

### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = septic_E, data2 =non_septic_E)

sep_non_E_MW = cbind(variable_to_compare, p_value)


## Timepoint F
### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = septic_F, data2 = non_septic_F)
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

sep_non_F_t = cbind(variable_to_compare, p_value, lambda, difference)

### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = septic_F, data2 =non_septic_F)

sep_non_F_MW = cbind(variable_to_compare, p_value)

#-----------------------------------------------------------------------
# SEPTIC VS NON_SEPTIC(MFI)
## Timepoint A
### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = septic_MFI_A, data2 = non_septic_MFI_A)
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

sep_non_MFI_A_t = cbind(variable_to_compare, p_value, lambda, difference)

### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = septic_MFI_A, data2 =non_septic_MFI_A)

sep_non_MFI_A_MW = cbind(variable_to_compare, p_value)

## Timepoint B
### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = septic_MFI_B, data2 = non_septic_MFI_B)
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

sep_non_MFI_B_t = cbind(variable_to_compare, p_value, lambda, difference)

### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = septic_MFI_B, data2 =non_septic_MFI_B)

sep_non_MFI_B_MW = cbind(variable_to_compare, p_value)

## Timepoint C
### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = septic_MFI_C, data2 = non_septic_MFI_C)
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

sep_non_MFI_C_t = cbind(variable_to_compare, p_value, lambda, difference)

### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = septic_MFI_C, data2 =non_septic_MFI_C)

sep_non_MFI_C_MW = cbind(variable_to_compare, p_value)

## Timepoint D
### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = septic_MFI_D, data2 = non_septic_MFI_D)
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

sep_non_MFI_D_t = cbind(variable_to_compare, p_value, lambda, difference)

### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = septic_MFI_D, data2 =non_septic_MFI_D)

sep_non_MFI_D_MW = cbind(variable_to_compare, p_value)

## Timepoint E
### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = septic_MFI_E, data2 = non_septic_MFI_E)
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

sep_non_MFI_E_t = cbind(variable_to_compare, p_value, lambda, difference)

### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = septic_MFI_E, data2 =non_septic_MFI_E)

sep_non_MFI_E_MW = cbind(variable_to_compare, p_value)

## Timepoint F
### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = septic_MFI_F, data2 = non_septic_MFI_F)
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

sep_non_MFI_F_t = cbind(variable_to_compare, p_value, lambda, difference)

### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = septic_MFI_F, data2 =non_septic_MFI_F)

sep_non_MFI_F_MW = cbind(variable_to_compare, p_value)

#-------------------------------------------------------------------
# SEPTIC VS HEALTHY(%)
### t-test
t_result = apply(variable_to_compare[-1, ], 1, compare_t, 
                 data1 = septic_A, data2 = healthy)
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

sep_hty_t = cbind(variable_to_compare[-1, ], p_value, lambda, difference)

### Mann-Whitney U test
p_value = apply(variable_to_compare[-1, ], 1, compare_MW, 
                data1 = septic_A, data2 = healthy)

sep_hty_MW = cbind(variable_to_compare[-1, ], p_value)

#-----------------------------------------------------------------------
# SEPTIC VS HEALTHY(MFI)
### t-test
t_result = apply(variable_to_compare[-1, ], 1, compare_t, 
                 data1 = septic_MFI_A, data2 = healthy_MFI)
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

sep_hty_MFI_t = cbind(variable_to_compare[-1, ], p_value, lambda, difference)

### Mann-Whitney U test
p_value = apply(variable_to_compare[-1, ], 1, compare_MW, 
                data1 = septic_MFI_A, data2 = healthy_MFI)

sep_hty_MFI_MW = cbind(variable_to_compare[-1, ], p_value)

#---------------------------------------------------------------------
# NON_SEPTIC VS HEALTHY(%)
### t-test
t_result = apply(variable_to_compare[-1, ], 1, compare_t, 
                 data1 = non_septic_A, data2 = healthy)
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

non_hty_t = cbind(variable_to_compare[-1, ], p_value, lambda, difference)

### Mann-Whitney U test
p_value = apply(variable_to_compare[-1, ], 1, compare_MW, 
                data1 = non_septic_A, data2 = healthy)

non_hty_MW = cbind(variable_to_compare[-1, ], p_value)

#-----------------------------------------------------------------------
# NON_SEPTIC VS HEALTHY(MFI)
### t-test
t_result = apply(variable_to_compare[-1, ], 1, compare_t, 
                 data1 = non_septic_MFI_A, data2 = healthy_MFI)
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

non_hty_MFI_t = cbind(variable_to_compare[-1, ], p_value, lambda, difference)

### Mann-Whitney U test
p_value = apply(variable_to_compare[-1, ], 1, compare_MW, 
                data1 = non_septic_MFI_A, data2 = healthy_MFI)

non_hty_MFI_MW = cbind(variable_to_compare[-1, ], p_value)

#---------------------------------------------------------------------
# Compare MALE with FEMALE
# The function to split male and female
gender_split = function(data){
  data_male = data %>% filter(Sex == 'M')
  data_female = data %>% filter(Sex == 'F')
  return (list(male = data_male, female = data_female))
}
#-------------------------------------------------------
## septic(%)
### Timepoint A
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(septic_A)[[1]], data2 = gender_split(septic_A)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_sep_A_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(septic_A)[[1]], data2 = gender_split(septic_A)[[2]])

gender_sep_A_MW = cbind(variable_to_compare, p_value)

### Timepoint B
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(septic_B)[[1]], data2 = gender_split(septic_B)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_sep_B_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(septic_B)[[1]], data2 = gender_split(septic_B)[[2]])

gender_sep_B_MW = cbind(variable_to_compare, p_value)

### Timepoint C
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(septic_C)[[1]], data2 = gender_split(septic_C)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_sep_C_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(septic_C)[[1]], data2 = gender_split(septic_C)[[2]])

gender_sep_C_MW = cbind(variable_to_compare, p_value)

### Timepoint D
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(septic_D)[[1]], data2 = gender_split(septic_D)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_sep_D_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(septic_D)[[1]], data2 = gender_split(septic_D)[[2]])

gender_sep_D_MW = cbind(variable_to_compare, p_value)

### Timepoint E
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(septic_E)[[1]], data2 = gender_split(septic_E)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_sep_E_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(septic_E)[[1]], data2 = gender_split(septic_E)[[2]])

gender_sep_E_MW = cbind(variable_to_compare, p_value)

### Timepoint F
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(septic_F)[[1]], data2 = gender_split(septic_F)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_sep_F_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(septic_F)[[1]], data2 = gender_split(septic_F)[[2]])

gender_sep_F_MW = cbind(variable_to_compare, p_value)
#--------------------------------------------------------------------------------------
## septic(MFI)
### Timepoint A
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(septic_MFI_A)[[1]], data2 = gender_split(septic_MFI_A)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_sep_MFI_A_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(septic_MFI_A)[[1]], data2 = gender_split(septic_MFI_A)[[2]])

gender_sep_MFI_A_MW = cbind(variable_to_compare, p_value)

### Timepoint B
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(septic_MFI_B)[[1]], data2 = gender_split(septic_MFI_B)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_sep_MFI_B_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(septic_MFI_B)[[1]], data2 = gender_split(septic_MFI_B)[[2]])

gender_sep_MFI_B_MW = cbind(variable_to_compare, p_value)

### Timepoint C
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(septic_MFI_C)[[1]], data2 = gender_split(septic_MFI_C)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_sep_MFI_C_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(septic_MFI_C)[[1]], data2 = gender_split(septic_MFI_C)[[2]])

gender_sep_MFI_C_MW = cbind(variable_to_compare, p_value)

### Timepoint D
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(septic_MFI_D)[[1]], data2 = gender_split(septic_MFI_D)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_sep_MFI_D_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(septic_MFI_D)[[1]], data2 = gender_split(septic_MFI_D)[[2]])

gender_sep_MFI_D_MW = cbind(variable_to_compare, p_value)

### Timepoint E
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(septic_MFI_E)[[1]], data2 = gender_split(septic_MFI_E)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_sep_MFI_E_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(septic_MFI_E)[[1]], data2 = gender_split(septic_MFI_E)[[2]])

gender_sep_MFI_E_MW = cbind(variable_to_compare, p_value)

### Timepoint F
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(septic_MFI_F)[[1]], data2 = gender_split(septic_MFI_F)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_sep_MFI_F_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(septic_MFI_F)[[1]], data2 = gender_split(septic_MFI_F)[[2]])

gender_sep_MFI_F_MW = cbind(variable_to_compare, p_value)

#------------------------------------------------------------------------------------
## non_septic(%)
### Timepoint A
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(non_septic_A)[[1]], data2 = gender_split(non_septic_A)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_non_A_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(non_septic_A)[[1]], data2 = gender_split(non_septic_A)[[2]])

gender_non_A_MW = cbind(variable_to_compare, p_value)

### Timepoint B
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(non_septic_B)[[1]], data2 = gender_split(non_septic_B)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_non_B_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(non_septic_B)[[1]], data2 = gender_split(non_septic_B)[[2]])

gender_non_B_MW = cbind(variable_to_compare, p_value)

### Timepoint C
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(non_septic_C)[[1]], data2 = gender_split(non_septic_C)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_non_C_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(non_septic_C)[[1]], data2 = gender_split(non_septic_C)[[2]])

gender_non_C_MW = cbind(variable_to_compare, p_value)

### Timepoint D
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(non_septic_D)[[1]], data2 = gender_split(non_septic_D)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_non_D_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(non_septic_D)[[1]], data2 = gender_split(non_septic_D)[[2]])

gender_non_D_MW = cbind(variable_to_compare, p_value)

### Timepoint E
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(non_septic_E)[[1]], data2 = gender_split(non_septic_E)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_non_E_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(non_septic_E)[[1]], data2 = gender_split(non_septic_E)[[2]])

gender_non_E_MW = cbind(variable_to_compare, p_value)

### Timepoint F
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(non_septic_F)[[1]], data2 = gender_split(non_septic_F)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_non_F_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(non_septic_F)[[1]], data2 = gender_split(non_septic_F)[[2]])

gender_non_F_MW = cbind(variable_to_compare, p_value)
#--------------------------------------------------------------------------------------
## non_septic(MFI)
### Timepoint A
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(non_septic_MFI_A)[[1]], data2 = gender_split(non_septic_MFI_A)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_non_MFI_A_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(non_septic_MFI_A)[[1]], data2 = gender_split(non_septic_MFI_A)[[2]])

gender_non_MFI_A_MW = cbind(variable_to_compare, p_value)

### Timepoint B
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(non_septic_MFI_B)[[1]], data2 = gender_split(non_septic_MFI_B)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_non_MFI_B_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(non_septic_MFI_B)[[1]], data2 = gender_split(non_septic_MFI_B)[[2]])

gender_non_MFI_B_MW = cbind(variable_to_compare, p_value)

### Timepoint C
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(non_septic_MFI_C)[[1]], data2 = gender_split(non_septic_MFI_C)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_non_MFI_C_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(non_septic_MFI_C)[[1]], data2 = gender_split(non_septic_MFI_C)[[2]])

gender_non_MFI_C_MW = cbind(variable_to_compare, p_value)

### Timepoint D
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(non_septic_MFI_D)[[1]], data2 = gender_split(non_septic_MFI_D)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_non_MFI_D_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(non_septic_MFI_D)[[1]], data2 = gender_split(non_septic_MFI_D)[[2]])

gender_non_MFI_D_MW = cbind(variable_to_compare, p_value)

### Timepoint E
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(non_septic_MFI_E)[[1]], data2 = gender_split(non_septic_MFI_E)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_non_MFI_E_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(non_septic_MFI_E)[[1]], data2 = gender_split(non_septic_MFI_E)[[2]])

gender_non_MFI_E_MW = cbind(variable_to_compare, p_value)

### Timepoint F
#### t-test
t_result = apply(variable_to_compare, 1, compare_t, 
                 data1 = gender_split(non_septic_MFI_F)[[1]], data2 = gender_split(non_septic_MFI_F)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_non_MFI_F_t = cbind(variable_to_compare, p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare, 1, compare_MW, 
                data1 = gender_split(non_septic_MFI_F)[[1]], data2 = gender_split(non_septic_MFI_F)[[2]])

gender_non_MFI_F_MW = cbind(variable_to_compare, p_value)
#----------------------------------------------------------------------------------------------------------
## Healthy (%)
#### t-test
t_result = apply(variable_to_compare[-1, ], 1, compare_t, 
                 data1 = gender_split(healthy)[[1]], data2 = gender_split(healthy)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_hty_t = cbind(variable_to_compare[-1, ], p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare[-1, ], 1, compare_MW, 
                data1 = gender_split(healthy)[[1]], data2 = gender_split(healthy)[[2]])

gender_hty_MW = cbind(variable_to_compare[-1, ], p_value)

## Healthy (MFI)
#### t-test
t_result = apply(variable_to_compare[-1, ], 1, compare_t, 
                 data1 = gender_split(healthy_MFI)[[1]], data2 = gender_split(healthy_MFI)[[2]])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

gender_hty_MFI_t = cbind(variable_to_compare[-1, ], p_value, lambda, difference)

#### Mann-Whitney U test
p_value = apply(variable_to_compare[-1, ], 1, compare_MW, 
                data1 = gender_split(healthy_MFI)[[1]], data2 = gender_split(healthy_MFI)[[2]])

gender_hty_MFI_MW = cbind(variable_to_compare[-1, ], p_value)
#----------------------------------------------------------------------
# Compare stim condtions for sepsis survivors and nonsurvivors
variable_stimulated = read_csv("Simulated_variable.csv", col_names = F)
variable_stimulated = as.data.frame(as.vector(t(variable_stimulated)))

# Septic(%)
## t-test
t_result = apply(variable_stimulated, 1, compare_t, 
                 data1 = septic_A[septic_A$`Survival (y/n)` == 'Y', ], data2 = septic_A[septic_A$`Survival (y/n)` == 'N', ])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

stimulated_t = cbind(variable_stimulated, p_value, lambda, difference)


#### Mann-Whitney U test
p_value = apply(variable_stimulated, 1, compare_MW, 
                data1 = septic_A[septic_A$`Survival (y/n)` == 'Y', ], data2 = septic_A[septic_A$`Survival (y/n)` == 'N', ])
stimulated_MW = cbind(variable_stimulated, p_value)

# Septic (MFI)
## t-test
t_result = apply(variable_stimulated, 1, compare_t, 
                 data1 = septic_MFI_A[septic_MFI_A$`Survival (y/n)` == 'Y', ], data2 = septic_MFI_A[septic_MFI_A$`Survival (y/n)` == 'N', ])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

stimulated_MFI_t = cbind(variable_stimulated, p_value, lambda, difference)


#### Mann-Whitney U test
p_value = apply(variable_stimulated, 1, compare_MW, 
                data1 = septic_MFI_A[septic_MFI_A$`Survival (y/n)` == 'Y', ], data2 = septic_MFI_A[septic_MFI_A$`Survival (y/n)` == 'N', ])
stimulated_MFI_MW = cbind(variable_stimulated, p_value)

stimulated_result = cbind(stimulated_t, stimulated_MW)
stimulated_MFI_result = cbind(stimulated_MFI_t, stimulated_MFI_MW)

library(xlsx)
write.xlsx(stimulated_result, file="Additional1.xlsx", sheetName="septic_A(%)", row.names=FALSE)
write.xlsx(stimulated_MFI_result, file="Additional1.xlsx", sheetName="septic_A(MFI)", append=TRUE, row.names=FALSE)
#-------------------------------------------------------
# Compare stim condtions by sau culture
sau_culture = read_csv("sau_culture.csv", col_names = T)

# Timepoint A
# Septic(%)
## t-test
t_result = apply(variable_stimulated, 1, compare_t, 
                 data1 = septic_A[sau_culture[, 2] == 'pos', ], data2 = septic_A[sau_culture[, 2] == 'neg', ])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

stimulated_sau_A_t = cbind(variable_stimulated, p_value, lambda, difference)


#### Mann-Whitney U test
p_value = apply(variable_stimulated, 1, compare_MW, 
                data1 = septic_A[sau_culture[, 2] == 'pos', ], data2 = septic_A[sau_culture[, 2] == 'neg', ])
stimulated_sau_A_MW = cbind(variable_stimulated, p_value)

# Septic (MFI)
## t-test
t_result = apply(variable_stimulated, 1, compare_t, 
                 data1 = septic_MFI_A[sau_culture[, 2] == 'pos', ], data2 = septic_MFI_A[sau_culture[, 2] == 'neg', ])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

stimulated_MFI_sau_A_t = cbind(variable_stimulated, p_value, lambda, difference)


#### Mann-Whitney U test
p_value = apply(variable_stimulated, 1, compare_MW, 
                data1 = septic_MFI_A[sau_culture[, 2] == 'pos', ], data2 = septic_MFI_A[sau_culture[, 2] == 'neg', ])
stimulated_MFI_sau_A_MW = cbind(variable_stimulated, p_value)

stimulated_sau_A = cbind(stimulated_sau_A_t, stimulated_sau_A_MW)
stimulated_MFI_sau_A = cbind(stimulated_MFI_sau_A_t, stimulated_MFI_sau_A_MW)

# Timepoint B
# Septic(%)
## t-test
t_result = apply(variable_stimulated, 1, compare_t, 
                 data1 = septic_B[sau_culture[, 2] == 'pos', ], data2 = septic_B[sau_culture[, 2] == 'neg', ])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

stimulated_sau_B_t = cbind(variable_stimulated, p_value, lambda, difference)


#### Mann-Whitney U test
p_value = apply(variable_stimulated, 1, compare_MW, 
                data1 = septic_B[sau_culture[, 2] == 'pos', ], data2 = septic_B[sau_culture[, 2] == 'neg', ])
stimulated_sau_B_MW = cbind(variable_stimulated, p_value)

# Septic (MFI)
## t-test
t_result = apply(variable_stimulated, 1, compare_t, 
                 data1 = septic_MFI_B[sau_culture[, 2] == 'pos', ], data2 = septic_MFI_B[sau_culture[, 2] == 'neg', ])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

stimulated_MFI_sau_B_t = cbind(variable_stimulated, p_value, lambda, difference)


#### Mann-Whitney U test
p_value = apply(variable_stimulated, 1, compare_MW, 
                data1 = septic_MFI_B[sau_culture[, 2] == 'pos', ], data2 = septic_MFI_B[sau_culture[, 2] == 'neg', ])
stimulated_MFI_sau_B_MW = cbind(variable_stimulated, p_value)

stimulated_sau_B = cbind(stimulated_sau_B_t, stimulated_sau_B_MW)
stimulated_MFI_sau_B = cbind(stimulated_MFI_sau_B_t, stimulated_MFI_sau_B_MW)


# Timepoint C
# Septic(%)
## t-test
t_result = apply(variable_stimulated, 1, compare_t, 
                 data1 = septic_C[sau_culture[, 2] == 'pos', ], data2 = septic_C[sau_culture[, 2] == 'neg', ])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

stimulated_sau_C_t = cbind(variable_stimulated, p_value, lambda, difference)


#### Mann-Whitney U test
p_value = apply(variable_stimulated, 1, compare_MW, 
                data1 = septic_C[sau_culture[, 2] == 'pos', ], data2 = septic_C[sau_culture[, 2] == 'neg', ])
stimulated_sau_C_MW = cbind(variable_stimulated, p_value)

# Septic (MFI)
## t-test
t_result = apply(variable_stimulated, 1, compare_t, 
                 data1 = septic_MFI_C[sau_culture[, 2] == 'pos', ], data2 = septic_MFI_C[sau_culture[, 2] == 'neg', ])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

stimulated_MFI_sau_C_t = cbind(variable_stimulated, p_value, lambda, difference)


#### Mann-Whitney U test
p_value = apply(variable_stimulated, 1, compare_MW, 
                data1 = septic_MFI_C[sau_culture[, 2] == 'pos', ], data2 = septic_MFI_C[sau_culture[, 2] == 'neg', ])
stimulated_MFI_sau_C_MW = cbind(variable_stimulated, p_value)

stimulated_sau_C = cbind(stimulated_sau_C_t, stimulated_sau_C_MW)
stimulated_MFI_sau_C = cbind(stimulated_MFI_sau_C_t, stimulated_MFI_sau_C_MW)


# Timepoint D
# Septic(%)
## t-test
t_result = apply(variable_stimulated, 1, compare_t, 
                 data1 = septic_D[sau_culture[, 2] == 'pos', ], data2 = septic_D[sau_culture[, 2] == 'neg', ])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

stimulated_sau_D_t = cbind(variable_stimulated, p_value, lambda, difference)


#### Mann-Whitney U test
p_value = apply(variable_stimulated, 1, compare_MW, 
                data1 = septic_D[sau_culture[, 2] == 'pos', ], data2 = septic_D[sau_culture[, 2] == 'neg', ])
stimulated_sau_D_MW = cbind(variable_stimulated, p_value)

# Septic (MFI)
## t-test
t_result = apply(variable_stimulated, 1, compare_t, 
                 data1 = septic_MFI_D[sau_culture[, 2] == 'pos', ], data2 = septic_MFI_D[sau_culture[, 2] == 'neg', ])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

stimulated_MFI_sau_D_t = cbind(variable_stimulated, p_value, lambda, difference)


#### Mann-Whitney U test
p_value = apply(variable_stimulated, 1, compare_MW, 
                data1 = septic_MFI_D[sau_culture[, 2] == 'pos', ], data2 = septic_MFI_D[sau_culture[, 2] == 'neg', ])
stimulated_MFI_sau_D_MW = cbind(variable_stimulated, p_value)

stimulated_sau_D = cbind(stimulated_sau_D_t, stimulated_sau_D_MW)
stimulated_MFI_sau_D = cbind(stimulated_MFI_sau_D_t, stimulated_MFI_sau_D_MW)


# Timepoint E
# Septic(%)
## t-test
t_result = apply(variable_stimulated, 1, compare_t, 
                 data1 = septic_E[sau_culture[, 2] == 'pos', ], data2 = septic_E[sau_culture[, 2] == 'neg', ])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

stimulated_sau_E_t = cbind(variable_stimulated, p_value, lambda, difference)


#### Mann-Whitney U test
p_value = apply(variable_stimulated, 1, compare_MW, 
                data1 = septic_E[sau_culture[, 2] == 'pos', ], data2 = septic_E[sau_culture[, 2] == 'neg', ])
stimulated_sau_E_MW = cbind(variable_stimulated, p_value)

# Septic (MFI)
## t-test
t_result = apply(variable_stimulated, 1, compare_t, 
                 data1 = septic_MFI_E[sau_culture[, 2] == 'pos', ], data2 = septic_MFI_E[sau_culture[, 2] == 'neg', ])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

stimulated_MFI_sau_E_t = cbind(variable_stimulated, p_value, lambda, difference)


#### Mann-Whitney U test
p_value = apply(variable_stimulated, 1, compare_MW, 
                data1 = septic_MFI_E[sau_culture[, 2] == 'pos', ], data2 = septic_MFI_E[sau_culture[, 2] == 'neg', ])
stimulated_MFI_sau_E_MW = cbind(variable_stimulated, p_value)

stimulated_sau_E = cbind(stimulated_sau_E_t, stimulated_sau_E_MW)
stimulated_MFI_sau_E = cbind(stimulated_MFI_sau_E_t, stimulated_MFI_sau_E_MW)


# Timepoint F
# Septic(%)
## t-test
t_result = apply(variable_stimulated, 1, compare_t, 
                 data1 = septic_F[sau_culture[, 2] == 'pos', ], data2 = septic_F[sau_culture[, 2] == 'neg', ])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

stimulated_sau_F_t = cbind(variable_stimulated, p_value, lambda, difference)


#### Mann-Whitney U test
p_value = apply(variable_stimulated, 1, compare_MW, 
                data1 = septic_F[sau_culture[, 2] == 'pos', ], data2 = septic_F[sau_culture[, 2] == 'neg', ])
stimulated_sau_F_MW = cbind(variable_stimulated, p_value)

# Septic (MFI)
## t-test
t_result = apply(variable_stimulated, 1, compare_t, 
                 data1 = septic_MFI_F[sau_culture[, 2] == 'pos', ], data2 = septic_MFI_F[sau_culture[, 2] == 'neg', ])
p_value = sapply(t_result, `[[`, 1)
lambda = sapply(t_result, `[[`, 2)
difference = sapply(t_result, `[[`, 3)

stimulated_MFI_sau_F_t = cbind(variable_stimulated, p_value, lambda, difference)


#### Mann-Whitney U test
p_value = apply(variable_stimulated, 1, compare_MW, 
                data1 = septic_MFI_F[sau_culture[, 2] == 'pos', ], data2 = septic_MFI_F[sau_culture[, 2] == 'neg', ])
stimulated_MFI_sau_F_MW = cbind(variable_stimulated, p_value)

stimulated_sau_F = cbind(stimulated_sau_F_t, stimulated_sau_F_MW)
stimulated_MFI_sau_F = cbind(stimulated_MFI_sau_F_t, stimulated_MFI_sau_F_MW)

#-------------------------------
stimulated_sau_A$Timepoint = 'A'
stimulated_MFI_sau_A$Timepoint = 'A'

stimulated_sau_B$Timepoint = 'B'
stimulated_MFI_sau_B$Timepoint = 'B'

stimulated_sau_C$Timepoint = 'C'
stimulated_MFI_sau_C$Timepoint = 'C'

stimulated_sau_D$Timepoint = 'D'
stimulated_MFI_sau_D$Timepoint = 'D'

stimulated_sau_E$Timepoint = 'E'
stimulated_MFI_sau_E$Timepoint = 'E'

stimulated_sau_F$Timepoint = 'F'
stimulated_MFI_sau_F$Timepoint = 'F'

stimulated_sau = rbind(stimulated_sau_A, stimulated_sau_B, stimulated_sau_C, stimulated_sau_D, stimulated_sau_E, stimulated_sau_F)
stimulated_sau_MFI = rbind(stimulated_MFI_sau_A, stimulated_MFI_sau_B, stimulated_MFI_sau_C, stimulated_MFI_sau_D, stimulated_MFI_sau_E, stimulated_MFI_sau_F)


library(xlsx)
write.xlsx(stimulated_sau, file="Additional2.xlsx", sheetName="septic(%)", row.names=FALSE)
write.xlsx(stimulated_sau_MFI, file="Additional2.xlsx", sheetName="septic(MFI)", append=TRUE, row.names=FALSE)

#----------------------------------------------------
# t-test VS Mann-Whitney U test
t_MW = function(res_t, res_MW){
  eno_obs_t = which(res_t[, 2] != 'Not enough observations')
  eno_obs_MW = which(res_MW[, 2] != 'Not enough observations')
  eno_obs = intersect(eno_obs_t, eno_obs_MW)
  n_total = length(eno_obs)
  if(n_total == 0){
    return (rep(0, 4))
  } else{

## t-test has smaller p-value
small_t = sum(as.numeric(as.character(res_t[eno_obs, 2])) < as.numeric(as.character(res_MW[eno_obs, 2])), na.rm = T)

## t-test significant but M-W test not
t_MW = sum(as.numeric(as.character(res_t[eno_obs, 2])) < 0.05 & as.numeric(as.character(res_MW[eno_obs, 2])) > 0.05, na.rm = T)

## MW test significant but t-test not
MW_t = sum(as.numeric(as.character(res_t[eno_obs, 2])) > 0.05 & as.numeric(as.character(res_MW[eno_obs, 2])) < 0.05, na.rm = T)

return (c(n_total, small_t, t_MW, MW_t))
}
}

tps = c('A', 'B', 'C', 'D', 'E', 'F')
com_grps = c('sep_non_', 'sep_non_MFI_', 'gender_sep_', 'gender_sep_MFI_', 'gender_non_', 'gender_non_MFI_')
n_total = small_t = t_MW0 = MW_t0 = 0
for(com_grp in com_grps){
  for (tp in tps){
    res_t = get(paste(com_grp, tp, '_t', sep = ''))
    res_MW = get(paste(com_grp, tp, '_MW', sep = ''))
    counts = t_MW(res_t, res_MW)
    n_total = n_total + counts[1]
    small_t = small_t + counts[2]
    t_MW0 = t_MW0 + counts[3]
    MW_t0 = MW_t0 + counts[4]
  }
}

t_MW_res = c(n_total, small_t, t_MW0, MW_t0)
t_MW_res = t_MW_res + t_MW(sep_hty_t, sep_hty_MW)
t_MW_res = t_MW_res + t_MW(sep_hty_MFI_t, sep_hty_MFI_MW)
t_MW_res = t_MW_res + t_MW(non_hty_t, non_hty_MW)
t_MW_res = t_MW_res + t_MW(non_hty_MFI_t, non_hty_MFI_MW)
t_MW_res = t_MW_res + t_MW(gender_hty_t, gender_hty_MW)
t_MW_res = t_MW_res + t_MW(gender_hty_MFI_t, gender_hty_MFI_MW)
write.csv(t_MW_res, 't_MW_res.csv')

#------------------------------------------------------------
# Data visulization
## T-tim-3 at Timepoint A
T_tim_3_1 = data.frame(`T-tim-3` = septic_A$`T-tim-3`)
T_tim_3_1$cohort = 'septic'
T_tim_3_2 = data.frame(`T-tim-3` =non_septic_A$`T-tim-3`)
T_tim_3_2$cohort = 'non_septic'
T_tim_3A = rbind(T_tim_3_1, T_tim_3_2)

ggplot(T_tim_3A, aes(x = cohort, y = T.tim.3)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5, stackratio = 1.5) +
  ylab('T-tim-3') + ggtitle('T-tim-3 at Timepoint A') +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_summary(fun.y=mean, geom="point", shape=18,
               size=3, color="red")

## WBC
WBC_1 = data.frame(WBC = septic$WBC, Timepoint = septic$Timepoint)
WBC_1$cohort = 'septic'
WBC_2 = data.frame(WBC =non_septic$WBC, Timepoint = non_septic$Timepoint)
WBC_2$cohort = 'non_septic'
WBC = rbind(WBC_1, WBC_2)

WBC_mean = WBC %>%
  group_by(Timepoint, cohort) %>% 
  summarise_at(vars("WBC"), mean, na.rm = T)      

ggplot(WBC_mean, aes(x=Timepoint, y=WBC, group = cohort, color=cohort)) + 
  geom_line() + geom_point() + ggtitle('WBC at each timepoint') +
  ylab('Population means of WBC') +
  theme(plot.title = element_text(hjust = 0.5))

## Septic WBC by gender
WBC_sep_A = septic_A[, which(names(septic_A) %in% c('WBC', 'Sex', 'T-tim-3'))]
WBC_sep_A$Sex = factor(WBC_sep_A$Sex, levels = c('M', 'F'))
WBC_non_A = non_septic_A[, which(names(non_septic_A) %in% c('WBC', 'Sex', 'T-tim-3'))]
WBC_hty = healthy[, which(names(he) %in% c('T-tim-3', 'Sex'))]

WBC_sep_plot = ggplot(WBC_sep_A[!is.na(WBC_sep_A$Sex),], aes(x = Sex, y = WBC)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5, stackratio = 1.5) +
  ylab('WBC') + ggtitle('WBC of septic patients at Timepoint A') +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_summary(fun.y=mean, geom="point", shape=18,
               size=5, color="red")

WBC_non_plot = ggplot(WBC_non_A[!is.na(WBC_non_A$Sex),], aes(x = Sex, y = WBC)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5, stackratio = 1.5) +
  ylab('WBC') + ggtitle('WBC of non-septic patients at Timepoint A') +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_summary(fun.y=mean, geom="point", shape=18,
               size=5, color="red")

library(cowplot)
plot_grid(WBC_sep_plot, WBC_non_plot, labels = "AUTO")

TT3_sep_plot = ggplot(WBC_sep_A[!is.na(WBC_sep_A$Sex),], aes(x = Sex, y = `T-tim-3`)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5, stackratio = 1.5) +
  ylab('WBC') + ggtitle('WBC of septic patients at Timepoint A') +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_summary(fun.y=mean, geom="point", shape=18,
               size=5, color="red")

TT3_non_plot = ggplot(WBC_non_A[!is.na(WBC_non_A$Sex),], aes(x = Sex, y = `T-tim-3`)) + 
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5, stackratio = 1.5) +
  ylab('WBC') + ggtitle('WBC of non-septic patients at Timepoint A') +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_summary(fun.y=mean, geom="point", shape=18,
               size=5, color="red")

plot_grid(TT3_sep_plot, TT3_non_plot, labels = "AUTO")

#--------------------------------------------
## MAIT % of total at Timepoint A
MAIT_1 = data.frame(`MAIT % of total` = septic_A$`MAIT % of total`)
MAIT_1$cohort = 'septic'
MAIT_2 = data.frame(`MAIT % of total` =non_septic_A$`MAIT % of total`)
MAIT_2$cohort = 'non_septic'
MAIT_0 = data.frame(`MAIT % of total` = healthy$`MAIT % of total`)
MAIT_0$cohort = 'healthy'
MAIT_A = rbind(MAIT_0, MAIT_1, MAIT_2)

MAIT_sur = data.frame(`MAIT % of total` = septic_A[septic_A$`Survival (y/n)` == 'Y', ]$`MAIT % of total`)
MAIT_death = data.frame(`MAIT % of total` = septic_A[septic_A$`Survival (y/n)` == 'N', ]$`MAIT % of total`)
MAIT_sur$cohort = 'septic_survivors'
MAIT_death$cohort = 'septic_deaths'
MAIT_A1 = rbind(MAIT_0, MAIT_sur, MAIT_death)
MAIT_A1$cohort = factor(MAIT_A1$cohort, levels = c('healthy', 'septic_survivors', 'septic_deaths'))

ggplot(MAIT_A, aes(x = cohort, y = log10(MAIT...of.total))) + 
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5, stackratio = 1.5) +
  ylab('log_10(MAIT % of total)') + ggtitle('The common logarithm of `MAIT % of total` at Timepoint A') +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20)) +
  stat_summary(fun.y=mean, geom="point", shape=18,
               size=3, color="red")

ggplot(MAIT_A1, aes(x = cohort, y = log10(MAIT...of.total))) + 
  geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 1.5, stackratio = 1.5) +
  ylab('log_10(MAIT % of total)') + ggtitle('The common logarithm of `MAIT % of total` at Timepoint A') +
  theme(plot.title = element_text(hjust = 0.5), text = element_text(size=20)) +
  stat_summary(fun.y=mean, geom="point", shape=18,
               size=3, color="red")
