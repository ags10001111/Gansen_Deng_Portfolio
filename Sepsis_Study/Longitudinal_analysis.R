setwd("C:/Users/ags10001111/OneDrive/2020 Summer/Master Project/Longitudinal analysis")
load("Longitudinal.RData")
library(tidyverse)

#---------------------------------------------------------------------
# Prepare the data
colnames(septic)[1] = "id"
colnames(septic_MFI)[1] = "id"
colnames(non_septic)[1] = "id"
colnames(non_septic_MFI)[1] = "id"

#---------------------------------------------------------------------
# Prepare functions
## The function for pair comparison
### Box-cox transformation function
powerTransform = function(y, lam1) {
  if (lam1 == 0) {
    return (log(y))
  } else {
    return ((y^lam1 - 1) / lam1)
  }
}

### t-test function
library(rcompanion)
compare_t1 = function(var_tc, data, timepoint){
  data_A = data %>% filter(Timepoint == 'A')
  data_T = data %>% filter(Timepoint == timepoint)
  x1 = get(var_tc, data_A)
  x2 = get(var_tc, data_T)
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


## The function for diagnosing time effect
library(geepack)
compare_gee = function(var_tc, data){
  # Construct the formula using the target varaible
  formula0 = as.formula(paste('`',var_tc, '`', "~ Timepoint", sep = ""))
  
  # Fit the gee model
  gee_model = geeglm(formula0, id = id, data = data, corstr="ar1")
  
  # Get the coefficients and p-values
  coef = gee_model$coefficients
  p_value = summary(gee_model)$coefficients[, 4]
  std = summary(gee_model)$coefficients[, 2]
  
  if(length(coef) < 6) {
    coef = c(coef, rep(NA, 6 - length(coef)))
    names(coef) = c("(Intercept)", "TimepointB",  "TimepointC",  "TimepointD", 
                    "TimepointE", "TimepointF")}
  if(length(p_value) < 6) {
    p_value = c(p_value, rep(NA, 6 - length(p_value)))
    names(p_value) = c("(Intercept)", "TimepointB",  "TimepointC",  "TimepointD", 
                    "TimepointE", "TimepointF")}
  
  if(length(std) < 6) {
    std = c(std, rep(NA, 6 - length(std)))
    names(std) = c("(Intercept)", "TimepointB",  "TimepointC",  "TimepointD", 
                       "TimepointE", "TimepointF")}
  
  return (list(p.value = p_value, coefficient = coef, standard_error = std))
}

#--------------------------------------------------------------------------------------
# Unstimulated PBMC samples
## Read the variables to compare
variable_unsimulated = read_csv("Unsimulated_variable.csv", col_names = F)


## Time effect
# Septic(%)
gee_result = apply(variable_unsimulated, 1, compare_gee, 
                 data = septic)
p_value_septic = sapply(gee_result, `[[`, 1)
coefficient_septic = sapply(gee_result, `[[`, 2)
std_septic = sapply(gee_result, `[[`, 3)

sep_unsti_pv = cbind(variable_unsimulated, t(p_value_septic))
colnames(sep_unsti_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
sep_unsti_coef = cbind(variable_unsimulated, t(coefficient_septic))
colnames(sep_unsti_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
sep_unsti_std = cbind(variable_unsimulated, t(std_septic))
colnames(sep_unsti_std) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')

write.csv(sep_unsti_coef, '222.csv')
sep_unsti = cbind(sep_unsti_pv, rep(' ', nrow(sep_unsti_coef)), sep_unsti_coef)

# Septic(MFI)
gee_result = apply(variable_unsimulated, 1, compare_gee, 
                   data = septic_MFI)
p_value_septic_MFI = sapply(gee_result, `[[`, 1)
coefficient_septic_MFI = sapply(gee_result, `[[`, 2)

sep_MFI_unsti_pv = cbind(variable_unsimulated, t(p_value_septic_MFI))
colnames(sep_MFI_unsti_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
sep_MFI_unsti_coef = cbind(variable_unsimulated, t(coefficient_septic_MFI))
colnames(sep_MFI_unsti_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')

sep_MFI_unsti = cbind(sep_MFI_unsti_pv, rep(' ', nrow(sep_MFI_unsti_coef)), sep_MFI_unsti_coef)


# non_septic(%)
gee_result = apply(variable_unsimulated, 1, compare_gee, 
                   data = non_septic)
p_value_non = sapply(gee_result, `[[`, 1)
coefficient_non = sapply(gee_result, `[[`, 2)

non_unsti_pv = cbind(variable_unsimulated, t(p_value_non))
colnames(non_unsti_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E')
non_unsti_coef = cbind(variable_unsimulated, t(coefficient_non))
colnames(non_unsti_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')

non_unsti = cbind(non_unsti_pv, rep(' ', nrow(non_unsti_coef)), non_unsti_coef)


# non_septic(MFI)
gee_result = apply(variable_unsimulated, 1, compare_gee, 
                   data = non_septic_MFI)
p_value_non_MFI = sapply(gee_result, `[[`, 1)
coefficient_non_MFI = sapply(gee_result, `[[`, 2)

non_MFI_unsti_pv = cbind(variable_unsimulated, t(p_value_non_MFI))
colnames(non_MFI_unsti_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
non_MFI_unsti_coef = cbind(variable_unsimulated, t(coefficient_non_MFI))
colnames(non_MFI_unsti_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')

non_MFI_unsti = cbind(non_MFI_unsti_pv, rep(' ', nrow(non_MFI_unsti_coef)), non_MFI_unsti_coef)


#------------------------------------------------------------------------
# Simulated PBMC sample
variable_simulated = read_csv("Simulated_variable.csv", col_names = F)
Ecoli = variable_simulated[, 1]
kpn = variable_simulated[, 2]
media = variable_simulated[, 3]
pae = variable_simulated[, 4]
sau = variable_simulated[, 5]
sty = variable_simulated[, 6]
ligand = variable_simulated[, 7]
add_1218 = variable_simulated[, 8]

## Ecoli
# Septic(%)
gee_result = apply(Ecoli, 1, compare_gee, 
                   data = septic)
p_value_Ecoli_septic = sapply(gee_result, `[[`, 1)
coefficient_Ecoli_septic = sapply(gee_result, `[[`, 2)

septic_Ecoli_pv = cbind(Ecoli, t(p_value_Ecoli_septic))
colnames(septic_Ecoli_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
septic_Ecoli_coef = cbind(Ecoli, t(coefficient_Ecoli_septic))
colnames(septic_Ecoli_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')



# Septic(MFI)
gee_result = apply(Ecoli, 1, compare_gee, 
                   data = septic_MFI)
p_value_Ecoli_septic_MFI = sapply(gee_result, `[[`, 1)
coefficient_Ecoli_septic_MFI = sapply(gee_result, `[[`, 2)

septic_MFI_Ecoli_pv = cbind(Ecoli, t(p_value_Ecoli_septic_MFI))
colnames(septic_MFI_Ecoli_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
septic_MFI_Ecoli_coef = cbind(Ecoli, t(coefficient_Ecoli_septic_MFI))
colnames(septic_MFI_Ecoli_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')



# non_septic(%)
gee_result = apply(Ecoli, 1, compare_gee, 
                   data = non_septic)
p_value_Ecoli_non = sapply(gee_result, `[[`, 1)
coefficient_Ecoli_non = sapply(gee_result, `[[`, 2)

non_Ecoli_pv = cbind(Ecoli, t(p_value_Ecoli_non))
colnames(non_Ecoli_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
non_Ecoli_coef = cbind(Ecoli, t(coefficient_Ecoli_non))
colnames(non_Ecoli_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# non_septic(MFI)
gee_result = apply(Ecoli, 1, compare_gee, 
                   data = non_septic_MFI)
p_value_Ecoli_non_MFI = sapply(gee_result, `[[`, 1)
coefficient_Ecoli_non_MFI = sapply(gee_result, `[[`, 2)

non_MFI_Ecoli_pv = cbind(Ecoli, t(p_value_Ecoli_non_MFI))
colnames(non_MFI_Ecoli_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
non_MFI_Ecoli_coef = cbind(Ecoli, t(coefficient_Ecoli_non_MFI))
colnames(non_MFI_Ecoli_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')



## kpn
# Septic(%)
gee_result = apply(kpn, 1, compare_gee, 
                   data = septic)
p_value_kpn_septic = sapply(gee_result, `[[`, 1)
coefficient_kpn_septic = sapply(gee_result, `[[`, 2)

septic_kpn_pv = cbind(kpn, t(p_value_kpn_septic))
colnames(septic_kpn_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
septic_kpn_coef = cbind(kpn, t(coefficient_kpn_septic))
colnames(septic_kpn_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# Septic(MFI)
gee_result = apply(kpn, 1, compare_gee, 
                   data = septic_MFI)
p_value_kpn_septic_MFI = sapply(gee_result, `[[`, 1)
coefficient_kpn_septic_MFI = sapply(gee_result, `[[`, 2)

septic_MFI_kpn_pv = cbind(kpn, t(p_value_kpn_septic_MFI))
colnames(septic_MFI_kpn_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
septic_MFI_kpn_coef = cbind(kpn, t(coefficient_kpn_septic_MFI))
colnames(septic_MFI_kpn_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# non_septic(%)
gee_result = apply(kpn, 1, compare_gee, 
                   data = non_septic)
p_value_kpn_non = sapply(gee_result, `[[`, 1)
coefficient_kpn_non = sapply(gee_result, `[[`, 2)

non_kpn_pv = cbind(kpn, t(p_value_kpn_non))
colnames(non_kpn_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
non_kpn_coef = cbind(kpn, t(coefficient_kpn_non))
colnames(non_kpn_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# non_septic(MFI)
gee_result = apply(kpn, 1, compare_gee, 
                   data = non_septic_MFI)
p_value_kpn_non_MFI = sapply(gee_result, `[[`, 1)
coefficient_kpn_non_MFI = sapply(gee_result, `[[`, 2)

non_MFI_kpn_pv = cbind(kpn, t(p_value_kpn_non_MFI))
colnames(non_MFI_kpn_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
non_MFI_kpn_coef = cbind(kpn, t(coefficient_kpn_non_MFI))
colnames(non_MFI_kpn_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')

## media
# Septic(%)
gee_result = apply(media, 1, compare_gee, 
                   data = septic)
p_value_media_septic = sapply(gee_result, `[[`, 1)
coefficient_media_septic = sapply(gee_result, `[[`, 2)

septic_media_pv = cbind(media, t(p_value_media_septic))
colnames(septic_media_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
septic_media_coef = cbind(media, t(coefficient_media_septic))
colnames(septic_media_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# Septic(MFI)
gee_result = apply(media, 1, compare_gee, 
                   data = septic_MFI)
p_value_media_septic_MFI = sapply(gee_result, `[[`, 1)
coefficient_media_septic_MFI = sapply(gee_result, `[[`, 2)

septic_MFI_media_pv = cbind(media, t(p_value_media_septic_MFI))
colnames(septic_MFI_media_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
septic_MFI_media_coef = cbind(media, t(coefficient_media_septic_MFI))
colnames(septic_MFI_media_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# non_septic(%)
gee_result = apply(media, 1, compare_gee, 
                   data = non_septic)
p_value_media_non = sapply(gee_result, `[[`, 1)
coefficient_media_non = sapply(gee_result, `[[`, 2)

non_media_pv = cbind(media, t(p_value_media_non))
colnames(non_media_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
non_media_coef = cbind(media, t(coefficient_media_non))
colnames(non_media_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# non_septic(MFI)
gee_result = apply(media, 1, compare_gee, 
                   data = non_septic_MFI)
p_value_media_non_MFI = sapply(gee_result, `[[`, 1)
coefficient_media_non_MFI = sapply(gee_result, `[[`, 2)

non_MFI_media_pv = cbind(media, t(p_value_media_non_MFI))
colnames(non_MFI_media_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
non_MFI_media_coef = cbind(media, t(coefficient_media_non_MFI))
colnames(non_MFI_media_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')

## pae
# Septic(%)
gee_result = apply(pae, 1, compare_gee, 
                   data = septic)
p_value_pae_septic = sapply(gee_result, `[[`, 1)
coefficient_pae_septic = sapply(gee_result, `[[`, 2)

septic_pae_pv = cbind(pae, t(p_value_pae_septic))
colnames(septic_pae_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
septic_pae_coef = cbind(pae, t(coefficient_pae_septic))
colnames(septic_pae_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# Septic(MFI)
gee_result = apply(pae, 1, compare_gee, 
                   data = septic_MFI)
p_value_pae_septic_MFI = sapply(gee_result, `[[`, 1)
coefficient_pae_septic_MFI = sapply(gee_result, `[[`, 2)

septic_MFI_pae_pv = cbind(pae, t(p_value_pae_septic_MFI))
colnames(septic_MFI_pae_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
septic_MFI_pae_coef = cbind(pae, t(coefficient_pae_septic_MFI))
colnames(septic_MFI_pae_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# non_septic(%)
gee_result = apply(pae, 1, compare_gee, 
                   data = non_septic)
p_value_pae_non = sapply(gee_result, `[[`, 1)
coefficient_pae_non = sapply(gee_result, `[[`, 2)

non_pae_pv = cbind(pae, t(p_value_pae_non))
colnames(non_pae_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
non_pae_coef = cbind(pae, t(coefficient_pae_non))
colnames(non_pae_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# non_septic(MFI)
gee_result = apply(pae, 1, compare_gee, 
                   data = non_septic_MFI)
p_value_pae_non_MFI = sapply(gee_result, `[[`, 1)
coefficient_pae_non_MFI = sapply(gee_result, `[[`, 2)

non_MFI_pae_pv = cbind(pae, t(p_value_pae_non_MFI))
colnames(non_MFI_pae_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
non_MFI_pae_coef = cbind(pae, t(coefficient_pae_non_MFI))
colnames(non_MFI_pae_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')

## sau
# Septic(%)
gee_result = apply(sau, 1, compare_gee, 
                   data = septic)
p_value_sau_septic = sapply(gee_result, `[[`, 1)
coefficient_sau_septic = sapply(gee_result, `[[`, 2)

septic_sau_pv = cbind(sau, t(p_value_sau_septic))
colnames(septic_sau_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
septic_sau_coef = cbind(sau, t(coefficient_sau_septic))
colnames(septic_sau_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# Septic(MFI)
gee_result = apply(sau, 1, compare_gee, 
                   data = septic_MFI)
p_value_sau_septic_MFI = sapply(gee_result, `[[`, 1)
coefficient_sau_septic_MFI = sapply(gee_result, `[[`, 2)

septic_MFI_sau_pv = cbind(sau, t(p_value_sau_septic_MFI))
colnames(septic_MFI_sau_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
septic_MFI_sau_coef = cbind(sau, t(coefficient_sau_septic_MFI))
colnames(septic_MFI_sau_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# non_septic(%)
gee_result = apply(sau, 1, compare_gee, 
                   data = non_septic)
p_value_sau_non = sapply(gee_result, `[[`, 1)
coefficient_sau_non = sapply(gee_result, `[[`, 2)

non_sau_pv = cbind(sau, t(p_value_sau_non))
colnames(non_sau_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
non_sau_coef = cbind(sau, t(coefficient_sau_non))
colnames(non_sau_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# non_septic(MFI)
gee_result = apply(sau, 1, compare_gee, 
                   data = non_septic_MFI)
p_value_sau_non_MFI = sapply(gee_result, `[[`, 1)
coefficient_sau_non_MFI = sapply(gee_result, `[[`, 2)

non_MFI_sau_pv = cbind(sau, t(p_value_sau_non_MFI))
colnames(non_MFI_sau_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
non_MFI_sau_coef = cbind(sau, t(coefficient_sau_non_MFI))
colnames(non_MFI_sau_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')

## sty
# Septic(%)
gee_result = apply(sty, 1, compare_gee, 
                   data = septic)
p_value_sty_septic = sapply(gee_result, `[[`, 1)
coefficient_sty_septic = sapply(gee_result, `[[`, 2)

septic_sty_pv = cbind(sty, t(p_value_sty_septic))
colnames(septic_sty_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
septic_sty_coef = cbind(sty, t(coefficient_sty_septic))
colnames(septic_sty_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# Septic(MFI)
gee_result = apply(sty, 1, compare_gee, 
                   data = septic_MFI)
p_value_sty_septic_MFI = sapply(gee_result, `[[`, 1)
coefficient_sty_septic_MFI = sapply(gee_result, `[[`, 2)

septic_MFI_sty_pv = cbind(sty, t(p_value_sty_septic_MFI))
colnames(septic_MFI_sty_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
septic_MFI_sty_coef = cbind(sty, t(coefficient_sty_septic_MFI))
colnames(septic_MFI_sty_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# non_septic(%)
gee_result = apply(sty, 1, compare_gee, 
                   data = non_septic)
p_value_sty_non = sapply(gee_result, `[[`, 1)
coefficient_sty_non = sapply(gee_result, `[[`, 2)

non_sty_pv = cbind(sty, t(p_value_sty_non))
colnames(non_sty_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
non_sty_coef = cbind(sty, t(coefficient_sty_non))
colnames(non_sty_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# non_septic(MFI)
gee_result = apply(sty, 1, compare_gee, 
                   data = non_septic_MFI)
p_value_sty_non_MFI = sapply(gee_result, `[[`, 1)
coefficient_sty_non_MFI = sapply(gee_result, `[[`, 2)

non_MFI_sty_pv = cbind(sty, t(p_value_sty_non_MFI))
colnames(non_MFI_sty_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
non_MFI_sty_coef = cbind(sty, t(coefficient_sty_non_MFI))
colnames(non_MFI_sty_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')

## ligand
# Septic(%)
gee_result = apply(ligand, 1, compare_gee, 
                   data = septic)
p_value_ligand_septic = sapply(gee_result, `[[`, 1)
coefficient_ligand_septic = sapply(gee_result, `[[`, 2)

septic_ligand_pv = cbind(ligand, t(p_value_ligand_septic))
colnames(septic_ligand_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
septic_ligand_coef = cbind(ligand, t(coefficient_ligand_septic))
colnames(septic_ligand_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# Septic(MFI)
gee_result = apply(ligand, 1, compare_gee, 
                   data = septic_MFI)
p_value_ligand_septic_MFI = sapply(gee_result, `[[`, 1)
coefficient_ligand_septic_MFI = sapply(gee_result, `[[`, 2)

septic_MFI_ligand_pv = cbind(ligand, t(p_value_ligand_septic_MFI))
colnames(septic_MFI_ligand_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
septic_MFI_ligand_coef = cbind(ligand, t(coefficient_ligand_septic_MFI))
colnames(septic_MFI_ligand_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# non_septic(%)
gee_result = apply(ligand, 1, compare_gee, 
                   data = non_septic)
p_value_ligand_non = sapply(gee_result, `[[`, 1)
coefficient_ligand_non = sapply(gee_result, `[[`, 2)

non_ligand_pv = cbind(ligand, t(p_value_ligand_non))
colnames(non_ligand_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
non_ligand_coef = cbind(ligand, t(coefficient_ligand_non))
colnames(non_ligand_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# non_septic(MFI)
gee_result = apply(ligand, 1, compare_gee, 
                   data = non_septic_MFI)
p_value_ligand_non_MFI = sapply(gee_result, `[[`, 1)
coefficient_ligand_non_MFI = sapply(gee_result, `[[`, 2)

non_MFI_ligand_pv = cbind(ligand, t(p_value_ligand_non_MFI))
colnames(non_MFI_ligand_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
non_MFI_ligand_coef = cbind(ligand, t(coefficient_ligand_non_MFI))
colnames(non_MFI_ligand_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')

## add_1218
# Septic(%)
gee_result = apply(add_1218, 1, compare_gee, 
                   data = septic)
p_value_add_1218_septic = sapply(gee_result, `[[`, 1)
coefficient_add_1218_septic = sapply(gee_result, `[[`, 2)

septic_add_1218_pv = cbind(add_1218, t(p_value_add_1218_septic))
colnames(septic_add_1218_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
septic_add_1218_coef = cbind(add_1218, t(coefficient_add_1218_septic))
colnames(septic_add_1218_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# Septic(MFI)
gee_result = apply(add_1218, 1, compare_gee, 
                   data = septic_MFI)
p_value_add_1218_septic_MFI = sapply(gee_result, `[[`, 1)
coefficient_add_1218_septic_MFI = sapply(gee_result, `[[`, 2)

septic_MFI_add_1218_pv = cbind(add_1218, t(p_value_add_1218_septic_MFI))
colnames(septic_MFI_add_1218_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
septic_MFI_add_1218_coef = cbind(add_1218, t(coefficient_add_1218_septic_MFI))
colnames(septic_MFI_add_1218_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# non_septic(%)
gee_result = apply(add_1218, 1, compare_gee, 
                   data = non_septic)
p_value_add_1218_non = sapply(gee_result, `[[`, 1)
coefficient_add_1218_non = sapply(gee_result, `[[`, 2)

non_add_1218_pv = cbind(add_1218, t(p_value_add_1218_non))
colnames(non_add_1218_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
non_add_1218_coef = cbind(add_1218, t(coefficient_add_1218_non))
colnames(non_add_1218_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')


# non_septic(MFI)
gee_result = apply(add_1218, 1, compare_gee, 
                   data = non_septic_MFI)
p_value_add_1218_non_MFI = sapply(gee_result, `[[`, 1)
coefficient_add_1218_non_MFI = sapply(gee_result, `[[`, 2)

non_MFI_add_1218_pv = cbind(add_1218, t(p_value_add_1218_non_MFI))
colnames(non_MFI_add_1218_pv) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')
non_MFI_add_1218_coef = cbind(add_1218, t(coefficient_add_1218_non_MFI))
colnames(non_MFI_add_1218_coef) = c('Variables', 'A', 'B', 'C', 'D', 'E', 'F')

#---------------------------------------------------------------------------------
sep_Ecoli = cbind(septic_Ecoli_pv, rep(' ', nrow(septic_Ecoli_coef)), septic_Ecoli_coef)
sep_MFI_Ecoli = cbind(septic_MFI_Ecoli_pv, rep(' ', nrow(septic_MFI_Ecoli_coef)), septic_MFI_Ecoli_coef)
non_Ecoli = cbind(non_Ecoli_pv, rep(' ', nrow(non_Ecoli_coef)), non_Ecoli_coef)
non_MFI_Ecoli = cbind(non_MFI_Ecoli_pv, rep(' ', nrow(non_MFI_Ecoli_coef)), non_MFI_Ecoli_coef)

sep_kpn = cbind(septic_kpn_pv, rep(' ', nrow(septic_kpn_coef)), septic_kpn_coef)
sep_MFI_kpn = cbind(septic_MFI_kpn_pv, rep(' ', nrow(septic_MFI_kpn_coef)), septic_MFI_kpn_coef)
non_kpn = cbind(non_kpn_pv, rep(' ', nrow(non_kpn_coef)), non_kpn_coef)
non_MFI_kpn = cbind(non_MFI_kpn_pv, rep(' ', nrow(non_MFI_kpn_coef)), non_MFI_kpn_coef)

sep_media = cbind(septic_media_pv, rep(' ', nrow(septic_media_coef)), septic_media_coef)
sep_MFI_media = cbind(septic_MFI_media_pv, rep(' ', nrow(septic_MFI_media_coef)), septic_MFI_media_coef)
non_media = cbind(non_media_pv, rep(' ', nrow(non_media_coef)), non_media_coef)
non_MFI_media = cbind(non_MFI_media_pv, rep(' ', nrow(non_MFI_media_coef)), non_MFI_media_coef)

sep_pae = cbind(septic_pae_pv, rep(' ', nrow(septic_pae_coef)), septic_pae_coef)
sep_MFI_pae = cbind(septic_MFI_pae_pv, rep(' ', nrow(septic_MFI_pae_coef)), septic_MFI_pae_coef)
non_pae = cbind(non_pae_pv, rep(' ', nrow(non_pae_coef)), non_pae_coef)
non_MFI_pae = cbind(non_MFI_pae_pv, rep(' ', nrow(non_MFI_pae_coef)), non_MFI_pae_coef)

sep_sau = cbind(septic_sau_pv, rep(' ', nrow(septic_sau_coef)), septic_sau_coef)
sep_MFI_sau = cbind(septic_MFI_sau_pv, rep(' ', nrow(septic_MFI_sau_coef)), septic_MFI_sau_coef)
non_sau = cbind(non_sau_pv, rep(' ', nrow(non_sau_coef)), non_sau_coef)
non_MFI_sau = cbind(non_MFI_sau_pv, rep(' ', nrow(non_MFI_sau_coef)), non_MFI_sau_coef)

sep_sty = cbind(septic_sty_pv, rep(' ', nrow(septic_sty_coef)), septic_sty_coef)
sep_MFI_sty = cbind(septic_MFI_sty_pv, rep(' ', nrow(septic_MFI_sty_coef)), septic_MFI_sty_coef)
non_sty = cbind(non_sty_pv, rep(' ', nrow(non_sty_coef)), non_sty_coef)
non_MFI_sty = cbind(non_MFI_sty_pv, rep(' ', nrow(non_MFI_sty_coef)), non_MFI_sty_coef)

sep_ligand = cbind(septic_ligand_pv, rep(' ', nrow(septic_ligand_coef)), septic_ligand_coef)
sep_MFI_ligand = cbind(septic_MFI_ligand_pv, rep(' ', nrow(septic_MFI_ligand_coef)), septic_MFI_ligand_coef)
non_ligand = cbind(non_ligand_pv, rep(' ', nrow(non_ligand_coef)), non_ligand_coef)
non_MFI_ligand = cbind(non_MFI_ligand_pv, rep(' ', nrow(non_MFI_ligand_coef)), non_MFI_ligand_coef)

sep_add_1218 = cbind(septic_add_1218_pv, rep(' ', nrow(septic_add_1218_coef)), septic_add_1218_coef)
sep_MFI_add_1218 = cbind(septic_MFI_add_1218_pv, rep(' ', nrow(septic_MFI_add_1218_coef)), septic_MFI_add_1218_coef)
non_add_1218 = cbind(non_add_1218_pv, rep(' ', nrow(non_add_1218_coef)), non_add_1218_coef)
non_MFI_add_1218 = cbind(non_MFI_add_1218_pv, rep(' ', nrow(non_MFI_add_1218_coef)), non_MFI_add_1218_coef)


#-------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------
# Export the data
library(xlsx)
write.xlsx(sep_unsti, file="longitudinal.xlsx", sheetName="sep_unsti(%)", row.names=FALSE)
write.xlsx(sep_MFI_unsti, file="longitudinal.xlsx", sheetName="sep_unsti(MFI)", append=TRUE, row.names=FALSE)
write.xlsx(non_unsti, file="longitudinal.xlsx", sheetName="non_unsti(%)", append=TRUE, row.names=FALSE)
write.xlsx(non_MFI_unsti, file="longitudinal.xlsx", sheetName="non_unsti(MFI)", append=TRUE, row.names=FALSE)

write.xlsx(sep_Ecoli, file="Ecoli.xlsx", sheetName="sep_Ecoli(%)", row.names=FALSE)
write.xlsx(sep_MFI_Ecoli, file="Ecoli.xlsx", sheetName="sep_Ecoli(MFI)", append=TRUE, row.names=FALSE)
write.xlsx(non_Ecoli, file="Ecoli.xlsx", sheetName="non_Ecoli(%)", append=TRUE, row.names=FALSE)
write.xlsx(non_MFI_Ecoli, file="Ecoli.xlsx", sheetName="non_Ecoli(MFI)", append=TRUE, row.names=FALSE)

write.xlsx(sep_kpn, file="kpn.xlsx", sheetName="sep_kpn(%)", row.names=FALSE)
write.xlsx(sep_MFI_kpn, file="kpn.xlsx", sheetName="sep_kpn(MFI)", append=TRUE, row.names=FALSE)
write.xlsx(non_kpn, file="kpn.xlsx", sheetName="non_kpn(%)", append=TRUE, row.names=FALSE)
write.xlsx(non_MFI_kpn, file="kpn.xlsx", sheetName="non_kpn(MFI)", append=TRUE, row.names=FALSE)

write.xlsx(sep_media, file="media.xlsx", sheetName="sep_media(%)", row.names=FALSE)
write.xlsx(sep_MFI_media, file="media.xlsx", sheetName="sep_media(MFI)", append=TRUE, row.names=FALSE)
write.xlsx(non_media, file="media.xlsx", sheetName="non_media(%)", append=TRUE, row.names=FALSE)
write.xlsx(non_MFI_media, file="media.xlsx", sheetName="non_media(MFI)", append=TRUE, row.names=FALSE)

write.xlsx(sep_pae, file="pae.xlsx", sheetName="sep_pae(%)", row.names=FALSE)
write.xlsx(sep_MFI_pae, file="pae.xlsx", sheetName="sep_pae(MFI)", append=TRUE, row.names=FALSE)
write.xlsx(non_pae, file="pae.xlsx", sheetName="non_pae(%)", append=TRUE, row.names=FALSE)
write.xlsx(non_MFI_pae, file="pae.xlsx", sheetName="non_pae(MFI)", append=TRUE, row.names=FALSE)

write.xlsx(sep_sau, file="sau.xlsx", sheetName="sep_sau(%)", row.names=FALSE)
write.xlsx(sep_MFI_sau, file="sau.xlsx", sheetName="sep_sau(MFI)", append=TRUE, row.names=FALSE)
write.xlsx(non_sau, file="sau.xlsx", sheetName="non_sau(%)", append=TRUE, row.names=FALSE)
write.xlsx(non_MFI_sau, file="sau.xlsx", sheetName="non_sau(MFI)", append=TRUE, row.names=FALSE)

write.xlsx(sep_sty, file="sty.xlsx", sheetName="sep_sty(%)", row.names=FALSE)
write.xlsx(sep_MFI_sty, file="sty.xlsx", sheetName="sep_sty(MFI)", append=TRUE, row.names=FALSE)
write.xlsx(non_sty, file="sty.xlsx", sheetName="non_sty(%)", append=TRUE, row.names=FALSE)
write.xlsx(non_MFI_sty, file="sty.xlsx", sheetName="non_sty(MFI)", append=TRUE, row.names=FALSE)

write.xlsx(sep_ligand, file="ligand.xlsx", sheetName="sep_ligand(%)", row.names=FALSE)
write.xlsx(sep_MFI_ligand, file="ligand.xlsx", sheetName="sep_ligand(MFI)", append=TRUE, row.names=FALSE)
write.xlsx(non_ligand, file="ligand.xlsx", sheetName="non_ligand(%)", append=TRUE, row.names=FALSE)
write.xlsx(non_MFI_ligand, file="ligand.xlsx", sheetName="non_ligand(MFI)", append=TRUE, row.names=FALSE)

write.xlsx(sep_add_1218, file="add_1218.xlsx", sheetName="sep_add_1218(%)", row.names=FALSE)
write.xlsx(sep_MFI_add_1218, file="add_1218.xlsx", sheetName="sep_add_1218(MFI)", append=TRUE, row.names=FALSE)
write.xlsx(non_add_1218, file="add_1218.xlsx", sheetName="non_add_1218(%)", append=TRUE, row.names=FALSE)
write.xlsx(non_MFI_add_1218, file="add_1218.xlsx", sheetName="non_add_1218(MFI)", append=TRUE, row.names=FALSE)

#--------------------------------------------------------------------
# Data Visualization
## Significant trends of septic patients
long1 = septic %>%
  group_by(Timepoint) %>% 
  summarise_at(vars('T-cd69', 'M-pd-1', 'M-lag-3'), mean, na.rm = T)      

ggplot(long1, aes(x=Timepoint)) + 
  geom_line(aes(y = `T-cd69`, group = 1, colour = "T-cd69")) + geom_point(aes(y = `T-cd69`)) +
  geom_line(aes(y = `M-pd-1`, group = 1, colour = "M-pd-1")) + geom_point(aes(y = `M-pd-1`)) +
  geom_line(aes(y = `M-lag-3`, group = 1, colour = "M-lag-3")) + geom_point(aes(y = `M-lag-3`)) +
  ggtitle('T-cd69, M-pd-1 and M-lag-3 at each timepoint (Septic)') +
  ylab('Population mean') +
  theme(plot.title = element_text(hjust = 0.5))

## M-cd69 and CD14 HLADR
MCD_1 = data.frame(`M-cd69` = septic$`M-cd69`, `CD14 HLADR` = septic$`CD14 HLADR`, Timepoint = septic$Timepoint)
MCD_1$cohort = 'septic'
MCD_2 = data.frame(`M-cd69` = non_septic$`M-cd69`, `CD14 HLADR` = non_septic$`CD14 HLADR`, Timepoint = non_septic$Timepoint)
MCD_2$cohort = 'non_septic'
MCD = rbind(MCD_1, MCD_2)

MCD_mean = MCD %>%
  group_by(Timepoint, cohort) %>% 
  summarise_at(vars("M.cd69", "CD14.HLADR"), mean, na.rm = T)      

MCD_plot1 = ggplot(MCD_mean, aes(x=Timepoint, y=M.cd69, group = cohort, color=cohort)) + 
  geom_line() + geom_point() + ggtitle('M-cd69 at each timepoint') +
  ylab('Mean of M-cd69') +
  theme(plot.title = element_text(hjust = 0.5))

MCD_plot2 = ggplot(MCD_mean, aes(x=Timepoint, y=CD14.HLADR, group = cohort, color=cohort)) + 
  geom_line() + geom_point() + ggtitle('CD14 HLADR at each timepoint') +
  ylab('Mean of CD14 HLADR') +
  theme(plot.title = element_text(hjust = 0.5))

library(cowplot)
plot_grid(MCD_plot1, MCD_plot2, labels = "AUTO")

#---------------------------------------------------
MCD_plot1 <- ggplot(septic, aes(x = Timepoint, y = `M-cd69`, group = septic, colour = septic)) + 
  geom_line(linetype = 'dashed') +
  stat_summary(fun.y=mean, aes(group=1), geom="line", colour="blue", size = 1.5) +
  stat_summary(fun.y=mean, aes(group=1), geom="point", colour="blue", size=3, shape=4) + 
  ggtitle('M-cd69 at each timepoint (septic)') +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

MCD_plot2 <-  ggplot(septic, aes(x = Timepoint, y = `CD14 HLADR`, group = septic, colour = septic)) + 
  geom_line(linetype = 'dashed') +
  stat_summary(fun.y=mean, aes(group=1), geom="line", colour="blue", size = 1.5) +
  stat_summary(fun.y=mean, aes(group=1), geom="point", colour="blue", size=3, shape=4) + 
  ggtitle('CD14 HLADR at each timepoint (septic)') +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

library(cowplot)
plot_grid(MCD_plot1, MCD_plot2, labels = "AUTO")

#-------------------------------------------
sss_plot1 <- ggplot(septic, aes(x = Timepoint, y = `T-cd69`, group = septic, colour = septic)) + 
  geom_line(linetype = 'dashed') +
  stat_summary(fun.y=mean, aes(group=1), geom="line", colour="blue", size = 1.5) +
  stat_summary(fun.y=mean, aes(group=1), geom="point", colour="blue", size=3, shape=4) + 
  ggtitle('T-cd69 at each timepoint (septic)') +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

sss_plot2 <-  ggplot(septic, aes(x = Timepoint, y = `M-pd-1`, group = septic, colour = septic)) + 
  geom_line(linetype = 'dashed') +
  stat_summary(fun.y=mean, aes(group=1), geom="line", colour="blue", size = 1.5) +
  stat_summary(fun.y=mean, aes(group=1), geom="point", colour="blue", size=3, shape=4) + 
  ggtitle('M-pd-1 at each timepoint (septic)') +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

sss_plot3 <-  ggplot(septic, aes(x = Timepoint, y = `M-lag-3`, group = septic, colour = septic)) + 
  geom_line(linetype = 'dashed') +
  stat_summary(fun.y=mean, aes(group=1), geom="line", colour="blue", size = 1.5) +
  stat_summary(fun.y=mean, aes(group=1), geom="point", colour="blue", size=3, shape=4) + 
  ggtitle('M-lag-3 at each timepoint (septic)') +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

library(cowplot)
plot_grid(sss_plot1, sss_plot2, sss_plot3, labels = "AUTO")



