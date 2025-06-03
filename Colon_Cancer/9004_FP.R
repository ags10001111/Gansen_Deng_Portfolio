library(survival)
library(tidyverse)
library(survminer)
library(corrplot)
library(ranger)
library(ggplot2)
library(ggfortify)

load("9004_FP.RData")

# Data preprocessing
# Proportion of complete cases
mean(complete.cases(colon))

# Remove first two columns
colon_data = colon[ ,c(-1, -2)]

# Remove observations with NA
colon_data = colon_data[complete.cases(colon), ]

# Split the data by event type
colon1 = colon_data %>% filter(etype == 1)
colon2 = colon_data %>% filter(etype == 2)

# Remove etype
colon1 = colon1[, -14]
colon2 = colon2[, -14]
#--------------------------------------------------------------
# EDA
## Correlation table
colon_numeric_rec = colon1[, c("age", "nodes", "time")]
colon_numeric_death = colon2[, c("age", "nodes", "time")]
res_rec = cor(colon_numeric_rec)
res_death = cor(colon_numeric_death)
corrplot(res_rec, type = "upper", order = "hclust", tl.col = "black", 
         tl.srt = 0)
corrplot(res_death, type = "upper", order = "hclust", tl.col = "black", 
         tl.srt = 0)

## Recurrence VS Death
colon_data = colon_data %>% mutate(etype = replace(etype, etype == "recurrence", "Recurrence")) %>%
  mutate(etype = replace(etype, etype == 2, "Death"))
rec_death_fit = survfit(Surv(time, status) ~ etype, data = colon_data)
rec_death_plot = autoplot(rec_death_fit, main = "Survival probability(Recurrence time VS Death time)")
rec_death_plot

## The effect of different treatments for recurrence and death
par(mfrow = c(1, 2))
rec_treat_fit = survfit(Surv(time, status) ~ rx, data = colon1)
rec_treat_plot = autoplot(rec_treat_fit, main = "S(t) under different treatments(Recurrence time)")
rec_treat_plot

death_treat_fit = survfit(Surv(time, status) ~ rx, data = colon2)
death_treat_plot = autoplot(death_treat_fit, main = "S(t) under different treatments(Death time)")
death_treat_plot

## Covariates change
aa_fit_rec = aareg(Surv(time, status) ~. , data = colon1)
rec_covariates_plot = autoplot(aa_fit_rec, main = "Change of covariates(Recurrence time)")
rec_covariates_plot

aa_fit_death = aareg(Surv(time, status) ~. , data = colon2)
death_covariates_plot = autoplot(aa_fit_death, main = "Change of covariates(Death time)")
death_covariates_plot

# The K-M estimate
KM1 = survfit(Surv(time, status) ~ 1, data = colon1, ctype = 1)
KM2 = survfit(Surv(time, status) ~ 1, data = colon2)
#--------------------------------------------------------------
# Model selection
plot(log(-log(KM1$surv)) ~ log(KM1$time))

#--------------------------------------------------------------
# Cox PH model(1)
## Check the nonlinearity for continuous variable nodes
ggcoxfunctional(Surv(time, status) ~ age + log(age) + sqrt(age), data = colon1)
ggcoxfunctional(Surv(time, status) ~ nodes + log(nodes+1) + sqrt(nodes), data = colon1)

## Fit a Cox PH model with all covariates
colon_cox = coxph(Surv(time, status) ~ rx+sex+log(age)+obstruct+perfor+
                    adhere+nodes+differ+extent+surg+node4 ,
                  data = colon1)
summary(colon_cox)

## Use AIC to do variable selection
colon_AIC = step(colon_cox)
summary(colon_AIC)

## Using sex, perfor and adhere as stratum
colon_cox2 = coxph(Surv(time, status) ~ rx+strata(sex)+obstruct+strata(perfor)+
                    strata(adhere)+nodes+differ+extent+surg+node4 ,
                  data = colon1)
summary(colon_cox2)

## Test for PH
cox.zph(colon_cox2)
ggcoxzph(cox.zph(colon_cox2))

var_name = c("rxLev" , "rxLev+5FU", "obstruct" ,  
             "nodes" , "differ" , "extent" , "surg" , "node4")
for(i in var_name){
  plot(cox.zph(colon_cox2), col = c('blue', 'red'), var = i)
  abline(a = 0, b = 0, col = 'green')
}

## Improve the model
colon_cox3 = coxph(Surv(time, status) ~ rx+strata(sex)+strata(obstruct)+strata(perfor)+
                     strata(adhere)+nodes+strata(differ)+extent+surg+strata(node4) ,
                   data = colon1)
summary(colon_cox3)

## Test the PH again
cox.zph(colon_cox3)
ggcoxzph(cox.zph(colon_cox3))

## Check for exponential link
ggcoxdiagnostics(colon_cox3, type = "deviance",
                linear.predictions = T, ggtheme = theme_bw())

## Obtain nonparametric baseline hazard estimate
colon_baseh = survfit(colon_cox3, type="aalen")
colon_baseh =survfit(colon_cox3, type="breslow")
plot(colon_baseh)

## Cumulative hazard function
plot(colon_baseh$time, -log(colon_baseh$surv), main="Breslow Estimate for CBH", type="s")

## Influential observation detection: estimated coefficient change divided by it's standard error
ggcoxdiagnostics(colon_cox3, type = "dfbeta", linear.predictions = TRUE)

#---------------------------------------------------------------
# Cox PH model(2)
## Check the nonlinearity for continuous variable nodes
ggcoxfunctional(Surv(time, status) ~ age + log(age) + sqrt(age), data = colon2)
ggcoxfunctional(Surv(time, status) ~ nodes + log(nodes+1) + sqrt(nodes), data = colon2)

## Fit a Cox PH model with all covariates
colon_cox_b = coxph(Surv(time, status) ~ rx+sex+log(age)+obstruct+perfor+
                    adhere+nodes+differ+extent+surg+node4 ,
                  data = colon2)
summary(colon_cox_b)

## Use AIC to do variable selection
colon_AIC_b = step(colon_cox_b)
summary(colon_AIC_b)

## Using sex, perfor and adhere as stratum
colon_cox_b2 = coxph(Surv(time, status) ~ rx+strata(sex)+ log(age) + obstruct+strata(perfor)+
                     strata(adhere)+nodes+differ+extent+surg+node4 ,
                   data = colon2)
summary(colon_cox_b2)

## Test for PH
cox.zph(colon_cox_b2)

## Improve the model
colon_cox_b3 = coxph(Surv(time, status) ~ rx+strata(sex)+log(age)+strata(obstruct)+strata(perfor)+
                     strata(adhere)+log(nodes+1)+strata(differ)+extent+surg+strata(node4) ,
                   data = colon2)
summary(colon_cox_b3)

## Test the PH again
cox.zph(colon_cox_b3)
ggcoxzph(cox.zph(colon_cox_b3))

var_name_b = c("rxLev" , "rxLev+5FU", "age", "nodes" , 
                 "extent" , "surg" )
for(i in var_name_b){
  plot(cox.zph(colon_cox_b2), col = c('blue', 'red'), var = i)
  abline(a = 0, b = 0, col = 'green')
}

## Check for exponential link
ggcoxdiagnostics(colon_cox_b3, type = "deviance",
                 linear.predictions = T, ggtheme = theme_bw())

## Inference
### Solve for the mean survival time
MST = function(x){
  X_pre = data.frame(rx = "Lev", sex = 1, age = 60, obstruct = 1,
                     perfor = 0, adhere = 0, nodes = 10, differ = 2,
                     extent = 2, surg = 1, node4 = 1, status = 1, time = x)
  S = exp(- predict(colon_cox_b3, X_pre, type = "expected"))
  return (S - 0.5)
}
confint(predict(colon_cox_b3, X_pre))

uniroot(MST, interval = c(0, 3329))


#---------------------------------------------------------------------------
# AFT model
## Test the exponential and weibull distribution
fit1<-survreg(Surv(time, status)~.,data=colon2)
cs.res1<-exp(-fit1$linear.predictor/fit1$scale)* (Surv(colon2$time, colon2$status)[,1])^(1/fit1$scale)
cs.fit1<-survfit(Surv(cs.res1,colon2$status)~1, type="fleming-harrington")
par(mfrow = c(1,2))
plot(cs.fit1$time, -log(cs.fit1$surv), type="s", 
     main = "-log(S(t)) vs. t for Exponential")
plot(log(cs.fit1$time), log(-log(cs.fit1$surv)), type="s", 
     main = "log(-log(S(t))) vs. log(t) for Weibull")

## Test the log-logistic distribution
fit <- survreg(Surv(time,status)~. , data = colon2 ,dist = "loglogistic")
cs.res <-  exp((log(Surv(colon2$time, colon2$status)[,1])-fit$linear.predictor)/fit$scale)
cs.fit <- survfit(Surv(cs.res,colon2$status)~1, type="fh2") 
plot(log(cs.fit$time), log(exp(-log(cs.fit$surv))-1), type="s", 
     main = "log(exp(-log(S(t)))) vs. log(t) for Loglogistic")

## Variable selection
fit2 <- survreg(Surv(time,status)~ rx+age+obstruct+perfor+adhere+nodes+differ+extent+surg+node4, data = colon2 , dist = "loglogistic")
fit3 <- survreg(Surv(time,status)~ rx+age+obstruct+adhere+nodes+differ+extent+surg+node4, data = colon2 , dist = "loglogistic")
fit4 <- survreg(Surv(time,status)~ rx+age+obstruct+nodes+differ+extent+surg+node4, data = colon2 , dist = "loglogistic")
fit5 <- survreg(Surv(time,status)~ rx+age+obstruct+nodes+extent+surg+node4, data = colon2 , dist = "loglogistic")

## Inference on T|X
coefficient <- as.matrix(fit5$coefficients, nrow=1 )
var_coef <- as.matrix(fit5$var)
var_coef[,10] <- var_coef[,10]*fit5$scale

## Inference on median death time
quan <- qlogis(0.5)
coef <- as.matrix(c(fit5$coefficients, fit5$scale), nrow = 1)
condition <- matrix(c(1,0,1,60,1,10,2,1,1), nrow = 1)
condition1 <- matrix(c(1,0,1,60,1,10,2,1,1,quan), nrow = 1)
y_median <- condition1 %*% coef
t_median <- exp(y_median)
var_y <- condition1 %*% var_coef %*% t(condition1)
CI_y <- c( y_median-1.96*sqrt(var_y),y_median+1.96*sqrt(var_y))
CI_t <- exp(CI_y)
var_coef[10,] <- var_coef[10,]*fit5$scale

## Get the point estimation and CI for S(t|X)
epsilon_hat <- (log(1500) - condition %*% coefficient)/fit5$scale
condition2 <- matrix(-c(1,0,1,60,1,10,2,1,1,epsilon_hat)/fit5$scale, nrow=1)
var_epsilon <- condition2 %*% var_coef %*% t(condition2)
CI_epsilon <- c(epsilon_hat-1.96*sqrt(var_epsilon),epsilon_hat+1.96*sqrt(var_epsilon))
CI_S <- 1/(1+exp(-CI_epsilon))
expect_S <- 1/(1+exp(-epsilon_hat))

## Model assessment
de.res<-residuals(fit, type="deviance")
par(mfrow = c(1,2))
plot(fit$linear.predictor, de.res)
plot(Surv(colon2$time, colon2$status)[,1], de.res)

## Analysis on recurrence data
### Check the error distribution
fit6<-survreg(Surv(time, status)~.,data=colon1)
cs.res2<-exp(-fit6$linear.predictor/fit6$scale)* (Surv(colon1$time, colon1$status)[,1])^(1/fit6$scale)
cs.fit2<-survfit(Surv(cs.res2,colon1$status)~1, type="fleming-harrington")
par(mfrow = c(1,2))
plot(cs.fit2$time, -log(cs.fit2$surv), type="s", 
     main = "-log(S(t)) vs. t for Exponential")
plot(log(cs.fit2$time), log(-log(cs.fit2$surv)), type="s", 
     main = "log(-log(S(t))) vs. log(t) for Weibull")

fit7 <- survreg(Surv(time,status)~. , data = colon1 ,dist = "loglogistic")
cs.res3 <-  exp((log(Surv(colon1$time, colon1$status)[,1])-fit7$linear.predictor)/fit7$scale)
cs.fit3 <- survfit(Surv(cs.res3,colon1$status)~1, type="fh2") 
plot(log(cs.fit3$time), log(exp(-log(cs.fit3$surv))-1), type="s", 
     main = "log(exp(-log(S(t)))) vs. log(t) for Loglogistic")

### Variable selection
summary(fit7)
fit8 <- survreg(Surv(time,status)~ rx+sex+obstruct+perfor+adhere+nodes+differ+extent+surg+node4, data = colon1 , dist = "loglogistic")
summary(fit8)
fit9 <- survreg(Surv(time,status)~ rx+sex+obstruct+adhere+nodes+differ+extent+surg+node4, data = colon1 , dist = "loglogistic")
summary(fit9)
fit10 <- survreg(Surv(time,status)~ rx+sex+obstruct+nodes+differ+extent+surg+node4, data = colon1 , dist = "loglogistic")
summary(fit10)
fit11 <- survreg(Surv(time,status)~ rx+obstruct+nodes+differ+extent+surg+node4, data = colon1 , dist = "loglogistic")
summary(fit11)
fit12 <- survreg(Surv(time,status)~ rx+obstruct+nodes+extent+surg+node4, data = colon1 , dist = "loglogistic")
summary(fit12)
