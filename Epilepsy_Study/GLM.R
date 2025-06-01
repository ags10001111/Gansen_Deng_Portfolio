setwd("C:/Users/ags10001111/OneDrive/2020 Summer/Epilepsy study")
library(tidyverse)
load('Epilepsy study_2013.RData')
#-----------------------------------------------------------
# Combining the data
GLM_data <- cbind(demo_data1[, -ncol(demo_data1)], 
                  clinical.info_data1[, -ncol(clinical.info_data1)],
                  cod_data1[, -ncol(cod_data1)],
                  neuro_data1[, -ncol(neuro_data1)],
                  AP_data1[, -ncol(AP_data1)],
                  demo_data1[, ncol(demo_data1)])

colnames(GLM_data)[ncol(GLM_data)] <- 'Cause of death:'

GLM_data$`Cause of death:` <- ifelse(GLM_data$`Cause of death:` == 'SUDEP', 1, 0)
GLM_data <- GLM_data[, !duplicated(colnames(GLM_data))]

# Remove columns with more than 40% missing values
na_prop <- function(vec){
  return (sum(is.na(vec))/length(vec))
}

na_pro <- apply(GLM_data, 2, na_prop)
plot(na_pro)
GLM_data_clean <- GLM_data[ ,which(na_pro<0.1)]

# Remove columns with only one level
num_lev <- function(vec){
  return(length(levels((as.factor(vec)))))
}

nol <- apply(GLM_data_clean[complete.cases(GLM_data_clean),], 2, num_lev)
GLM_data_clean <- GLM_data_clean[, which(nol > 1)]

# Convert character to factor
GLM_data_clean <- mutate_if(GLM_data_clean, is.character, as.factor)


sum(complete.cases(GLM_data_clean)) # Only 91 complete observations out of 146 (How to deal with missing values)
#------------------------------------------------------------
# Fit logistic regression models
## Complete case analysis
full_model <- glm(`Cause of death:` ~ . , data = GLM_data_clean[complete.cases(GLM_data_clean), ], family = 'binomial')

### Step wise regression
small_model <- step(full_model, k = log(sum(complete.cases(GLM_data_clean))))
ca_sw_result <- small_model$coefficients
summary(small_model)

### LASSO
library(glmnet)
x_train <- model.matrix(~ ., GLM_data_clean[complete.cases(GLM_data_clean), which(colnames(GLM_data_clean) != 'Cause of death:')])
lasso_model <- glmnet(x_train,
    GLM_data_clean[complete.cases(GLM_data_clean), which(colnames(GLM_data_clean) == 'Cause of death:')],
    family = 'binomial', intercept = T)

lasso.beta = lasso_model$beta
obj = -deviance(lasso_model)
k = lasso_model$df
n = lasso_model$nobs

BIC_lasso = log(n)*k - obj

lambda.lasso = which.min(BIC_lasso)
lasso.beta = lasso.beta[, lambda.lasso]
ca_lasso_result <- lasso.beta[lasso.beta!= 0]

### Random forest
library(randomForest)
GLM_data_clean1 <- GLM_data_clean
colnames(GLM_data_clean1) <- make.names(colnames(GLM_data_clean1))
ca_rf <- randomForest(Cause.of.death. ~ .,  data= GLM_data_clean1, 
                     na.action = na.omit, n.tree = 500, importance = T)
varImpPlot(ca_rf, type = 1)

#------------------------------------------------------
# Deal with missing values
addq <- function(x) paste0("`", x, "`")
my_formula = as.formula(paste("`Cause of death:` ~", paste(addq(colnames(GLM_data_clean)[!colnames(GLM_data_clean) %in% "Cause of death:"]), collapse = " + ")))

library(Hmisc)
data_imputed = aregImpute(my_formula, 
                          data = GLM_data_clean, n.impute = 5, type = 'regression')
