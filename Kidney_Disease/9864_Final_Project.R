library(tibble)
library(xtable)
library(tidyverse)
library(e1071)
library(randomForest)
library(varImp)
library(caret)
library(pROC)

# Read the data
kidney_data = read_csv("Original_data.csv", na = c("", "NA", "?"))
kidney_data = read_csv("chronic_kidney_disease.csv", na = c("", "NA", "?"))
as.tibble(kidney_data)
# make colnames meaningful
colnames(kidney_data) = gsub("'", "", colnames(kidney_data))
#----------------------------------------------------------------------
# Change variables' class
factor_no = c(3, 4, 5, 6, 7, 8, 9, 19:25) + 1
factor_name = colnames(kidney_data)[factor_no]
numeric_name = colnames(kidney_data)[-c(1, factor_no)]
kidney_data = kidney_data %>% mutate_each_(funs(factor), factor_name) %>% mutate_each_(funs(as.numeric), numeric_name) 


# Remove column ID
data1 = kidney_data[, -1]

# Listwise deletion
data2 = drop_na(data1)
#-------------------------------------------------------------------------------
# Packages that deal with missing data
library(mice)
library(VIM)

# Get the distribution of NA
mice_plot <- aggr(data1, col=c('navyblue','yellow'),
                    numbers=TRUE, sortVars=TRUE,
                    labels=names(data1), cex.axis= 0.3,
                    gap=3, ylab=c("Missing data","Pattern"))
# Get the count of NA in each variable
NA_count = arrange(mice_plot$missings, Count)

# Plot the distribution of NA by variable
plot_NA_count = ggplot(NA_count, aes(x = reorder(Variable, -Count), y = Count)) + 
  geom_bar(position="dodge", stat="identity") + xlab("Variables") +
  ylab("NA count") + ggtitle("NA distribution by variables")

# Plot the distribution of NA in dataframe
NA_dist = as.data.frame(ifelse(is.na(data1), 1, 0))
colnames(NA_dist) = colnames(data1)
NA_dist$row = 1:nrow(NA_dist)
NA_dist1 = NA_dist[1:250, ]
NA_dist2 = NA_dist[251:400, ]

NA_dist = gather(NA_dist, col, value, -row) %>% 
  mutate(col = factor(col, levels=colnames(NA_dist)))
NA_dist1 = gather(NA_dist1, col, value, -row) %>% 
  mutate(col = factor(col, levels=colnames(NA_dist1)))
NA_dist2 = gather(NA_dist2, col, value, -row) %>% 
  mutate(col = factor(col, levels=colnames(NA_dist2)))

plot_NA_dist1 = ggplot(NA_dist1, aes(col, row, fill=factor(value))) +
  geom_tile(colour="grey50") +
  scale_fill_manual(values=c("1"="green", "0"="white")) +
  scale_y_reverse(breaks=1:50, expand=c(0,0)) +
  scale_x_discrete(position="top") +
  labs(fill="NA") +
  theme_classic() + ggtitle("NA distribution of class = 'ckd'")

plot_NA_dist2 = ggplot(NA_dist2, aes(col, row, fill=factor(value))) +
  geom_tile(colour="grey50") +
  scale_fill_manual(values=c("1"="green", "0"="white")) +
  scale_y_reverse(breaks=1:50, expand=c(0,0)) +
  scale_x_discrete(position="top") +
  labs(fill="NA") +
  theme_classic() + ggtitle("NA distribution of class = 'notckd'")

# Drop rbc, rbcc and wbcc
drop.cols <- c('rbc', 'rbcc', 'wbcc')
data3 = data1 %>% select(-one_of(drop.cols))


# Diagnose the missing data mechanism
## Check MCAR
library(BaylorEdPsych)
little_test = LittleMCAR(as.data.frame(data3))

# Multiple Imputation
data4 = data3
al5 = (data4$al == 5)
al5[is.na(al5)] <- FALSE
data4[al5,]$al = 4
data4$al = factor(data4$al)

su5 = (data4$su == 5)
su5[is.na(su5)] <- FALSE
data4[su5,]$su = 4
data4$su = factor(data4$su)

library(Hmisc)
data_imputed = aregImpute(class ~ age + bp + al + su + pc + sg + pcc + ba + bgr + bu + 
             sc + sod + pot + hemo + pcv + htn + dm + cad + appet + pe + 
             ane, 
                             data = data4, n.impute = 5, type = 'pmm')
fill_data <- function(impute = data_imputed, data = data4, im = 5) {
  cbind.data.frame(impute.transcan(x = impute, 
                                   imputation = im, 
                                   data = data, 
                                   list.out = TRUE, 
                                   pr = FALSE))
}
kidney_imputed <- fill_data()
levels(kidney_imputed$al) = c(levels(kidney_imputed$al), 5)
levels(kidney_imputed$su) = c(levels(kidney_imputed$su), 5)
kidney_imputed[al5,]$al = 5
kidney_imputed[su5,]$su = 5

# The kidney imputed data is the cleaned data
#-------------------------------------------------------------------
# EDA
# Correlation table
load("9864_FP.RData")
library(corrplot)
kidney_imputed_numeric = dplyr::select_if(kidney_imputed, is.numeric) 
res = cor(kidney_imputed_numeric)
corrplot(res, method = "pie", type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)

# PCA
pca = prcomp(kidney_imputed_numeric, center = TRUE, scale = TRUE)
summary(pca)

# MCA
library(FactoMineR)
library(factoextra)
kidney_imputed_factor = dplyr::select_if(kidney_imputed, is.factor)
res.mca = MCA(kidney_imputed_factor) # multiple correspondence analysis
mca_plot = fviz_mca_var(res.mca, repel = TRUE, # Avoid text overlapping (slow)
                        ggtheme = theme_minimal(), col.var = "coral")

# Random Forest
set.seed(200)
rf = randomForest(class ~ ., data = kidney_imputed, importance = TRUE)
varimp_plot = varImpPlot(rf,type=2)
#--------------------------------------------------------------------------------
# Modelling
# Logistic regression
## Split the data into training set and testing set
kidney_imputed$class = ifelse(kidney_imputed$class == "ckd", 1, 0)
set.seed(222)
index = sample(nrow(kidney_imputed), 0.7*nrow(kidney_imputed))
kidney_train = kidney_imputed[index, ]
kidney_test = kidney_imputed[-index, ]

## Fit a logistic regression model including all variables
model_logistic = glm(class ~ ., data = kidney_train, family = 'binomial')
## Use AIC to do the stepwise regression
AIC_logisitc = step(model_logistic)
summary(AIC_logisitc)
pchisq(AIC_logisitc$deviance, df=AIC_logisitc$df.residual, lower.tail=FALSE)

## Draw the ROC curve
prob_logit0=predict(AIC_logisitc,kidney_test,type='response')
library(pROC)
roc_logit0=roc(response=kidney_test$class, predictor=prob_logit0)
plot(roc_logit0, main = "ROC of logistic")
auc(roc_logit0)

## Fit another logistic regression model with hemo only
model_logistic1 = glm(class ~ hemo, data = kidney_train, family = 'binomial')
summary(model_logistic1)

## Fit another logistic regression model with pcv only
model_logistic2 = glm(class ~ pcv, data = kidney_train, family = 'binomial')
summary(model_logistic2)

## Make prediction with hemo model
prob_logit1=predict(model_logistic1,kidney_test,type='response')
## Draw the ROC curve
library(pROC)
roc_logit1=roc(response=kidney_test$class, predictor=prob_logit1)
plot(roc_logit1, main="ROC curve with Hemo only")
auc(roc_logit1)

## Make prediction with pcv model
prob_logit2=predict(model_logistic2,kidney_test,type='response')
## Draw the ROC curve
library(pROC)
roc_logit2=roc(response=kidney_test$class, predictor=prob_logit2)
plot(roc_logit2, main="ROC curve with PCV only")
auc(roc_logit2)


# SVM
set.seed(222)
train = sample(1:400, 280, replace = FALSE)
kidney_imputed_train = kidney_imputed[train, ]
kidney_imputed_test = kidney_imputed[-train, ]
fit = svm(class ~ ., data = kidney_imputed_train, scale = FALSE, kernel = "radial", cost = 5)
predict_svm = predict(fit, newdata = kidney_imputed_test)
conf_mat_svm = table(predicted = predict_svm, actual = kidney_imputed_test$class)
svm_matrix = confusionMatrix(conf_mat_svm)

svm_data = kidney_imputed %>% mutate_each_(funs = as.numeric, 1)
svm_data$class = svm_data$class - 1
fit_2 = svm(class ~ ., data = svm_data[train, ], scale = FALSE, kernel = "radial", cost = 5)
predict_svm_2 = predict(fit_2, newdata = svm_data[-train, ])
roc_svm = roc(svm_data[-train, ]$class, predict_svm_2, plot = TRUE, col = "orange")
auc_svm = auc(roc_svm)

#-----------------------------------------------------------------
# Prediction
# Use 5-fold CV to find the best cutoff that maximize accuracy
## Prepare the data for cross validation
set.seed(850)
rand_index = sample(nrow(kidney_train))
kidney_cv = kidney_train[rand_index,]
folds = cut(1:nrow(kidney_train),breaks=5,labels=FALSE)

## A function that compute average accuracy
logit_cut = function(c){
  acc_logit = numeric(5)
  for(i in 1:5){
    test_index = which(folds==i)
    cv_test = kidney_cv[test_index,]
    cv_train = kidney_cv[-test_index,]
    fit_logit = glm(class ~ age + sg + bgr + hemo + pcv + dm + appet,
                    data = cv_train, family = 'binomial')
    prob_logit = predict(fit_logit,newdata=cv_test,type="response")
    pred_logit = ifelse(prob_logit > c, 1, 0)
    conf_mat_logit = table(predicted = pred_logit, actual = cv_test$class)
    cof_table = confusionMatrix(conf_mat_logit, positive = "1")
    acc_logit[i] = cof_table$overall[1]
  }
  return (mean(acc_logit))
}
## Find the optimum cutoff
cut_max=optimize(logit_cut,interval=c(0,1), maximum = T, tol = 0.00001)

## Use that cutoff to make prediction
pred_logit0=ifelse(prob_logit0>0.976, 1, 0)
## The confusion matrix
library(caret)
conf_mat_logit0 = table(predicted = pred_logit0, actual = kidney_test$class)
confusionMatrix(conf_mat_logit0, positive = "1")

#------------------------------------------------------------------------------
# Decision tree
library(rpart)   #Decision tree
library(rattle)  #Fancy tree plot
## Build a tree model and plot it
kidney_tree0=rpart(class~.,data=kidney_imputed,method='class')
fancyRpartPlot(kidney_tree0, sub = "Decision tree")

