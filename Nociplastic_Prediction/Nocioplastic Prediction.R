library(readxl)
library(tidyverse)
library(pROC)
library(caret)
library(randomForest)
library(glmnet)
library(hrbrthemes)
library(e1071)
library(MASS)

# Data Preprocessing--------------------------------------------------------------------------
df0 <- read.csv('SPARC_PPT.csv')

## Change the column type to numeric
df <- df0 %>% mutate_if(is.character, as.numeric)

## Only include the complete cases 
#!!! Note: Missing Values can be imputed but not sure how PPT should be imputed !!!#
df_clean <- df[complete.cases(df),]

## Encode the PMD to Nociplastic (Binary)
df_clean <- df_clean %>% 
  mutate(Nociplastic = case_when(
    `PMD` %in% c(3,5,6,7) ~ 1,
    `PMD` %in% c(1,2,4) ~ 0
  ))

df_clean$Nociplastic <- as.factor(df_clean$Nociplastic)

## The Normal PPT Values
normal_PPT_M <- data.frame(t(c('Trapezius Left' = 47.62, 'Trapezius Right' = 51.96,
                             'Levator Left' = 56.86, 'Levator Right' = 53.92,
                             'Latissimus Left' = 62.75, 'Latissimus Right' = 61.76,
                             'Multifidus Left' = 76.47, 'Multifidus Right' = 84.31,
                             'Gluteus medius Left' = 68.63, 'Gluteus medius Right' = 65.69,
                             'Piriformis Left' = 68.63, 'Piriformis Right' = 65.69)))
normal_PPT_M <- do.call("rbind", replicate( 
  sum(df_clean$Sex == 1), normal_PPT_M, simplify = FALSE)) 

normal_PPT_F <- data.frame(t(c('Trapezius Left' = 36.27, 'Trapezius Right' = 35.29,
                             'Levator Left' = 44.12, 'Levator Right' = 44.12,
                             'Latissimus Left' = 62.75, 'Latissimus Right' = 61.76,
                             'Multifidus Left' = 55.88, 'Multifidus Right' = 58.82,
                             'Gluteus medius Left' = 62.75, 'Gluteus medius Right' = 65.69,
                             'Piriformis Left' = 62.75, 'Piriformis Right' = 65.69)))
normal_PPT_F <- do.call("rbind", replicate( 
  sum(df_clean$Sex == 2), normal_PPT_F, simplify = FALSE)) 

df_PPT <- df_clean[,colnames(normal_PPT_M)]
df_PPT_diff <- df_PPT
df_PPT_diff[df_clean$Sex == 1,] <- df_PPT[df_clean$Sex == 1,] - normal_PPT_M
df_PPT_diff[df_clean$Sex == 2,] <- df_PPT[df_clean$Sex == 2,] - normal_PPT_F
df_PPT_diff[df_PPT == -1] <- 0

#df_pre <- cbind(df_PPT_diff, df_clean[,21:23])

if (T){
df_PPT_diff$APD <- rowMeans(replace(df_PPT_diff, df_PPT_diff <= 0, NA), na.rm = TRUE)
df_PPT_diff$AND <- abs(rowMeans(replace(df_PPT_diff, df_PPT_diff >= 0, NA), na.rm = TRUE))
df_PPT_diff[is.na(df_PPT_diff)] <- 0
df_pre <- cbind(df_PPT_diff[,13:14], df_clean[,21:23])}
colnames(df_pre)[3] <- 'MPD'
df_pre <- df_pre[df_pre$APD < 100,]

# Calculate the skewness statistics
skewness(df_pre$APD)
skewness(df_pre$AND)
skewness(df_pre$MPD)
skewness(df_pre$CSI)

cor(df_pre$CSI, df_pre$MPD)
cor(df_pre$MPD, df_pre$AND)
cor(df_pre$MPD, df_pre$APD)

# Box-Cox Transformation
if (F){
  b <- boxcox(lm(df_pre$APD + 1 ~ 1))
  bl <- b$x[which.max(b$y)]
  df_pre$APD <- ((df_pre$APD + 1) ^ bl - 1) / bl
  
  b <- boxcox(lm(df_pre$AND + 1 ~ 1))
  bl <- b$x[which.max(b$y)]
  df_pre$AND <- ((df_pre$AND + 1) ^ bl - 1) / bl
  
  b <- boxcox(lm(df_pre$MPD + 1 ~ 1))
  bl <- b$x[which.max(b$y)]
  df_pre$MPD <- ((df_pre$MPD + 1) ^ bl - 1) / bl
}

write.csv(df_pre, 'df_pre_orig.csv')
df_pre$APD <- log(df_pre$APD + 1)
df_pre$AND <- log(df_pre$AND + 1)
df_pre$MPD <- log(df_pre$MPD + 1)
skewness(df_pre$APD)
skewness(df_pre$AND)
skewness(df_pre$MPD)
skewness(df_pre$CSI)

# EDA --------------------------
panel.hist <- function(x, ...)
{
  usr <- par("usr")
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", border = "white")
}
pairs(
  df_pre[1:4],            # Columns to include
  col = df_pre$Nociplastic,
  diag.panel = panel.hist
      # Color by group
)

write.csv(df_pre, 'df_pre_clean.csv')
cor(df_pre$CSI, as.numeric(df_pre$Nociplastic))

# Predict MPD using CSI & PPT --------------------------------------------
n = 100
lr_acc <- numeric(n)
rf_acc <- numeric(n)
svml_acc <- numeric(n)
svmr_acc <- numeric(n)

for (i in 1:n){
  ## Split train and test
  set.seed(i)
  trn_idx <- sample(1:nrow(df_pre), 0.7*nrow(df_pre))
  train <- df_pre[trn_idx,]
  test <- df_pre[-trn_idx,]
  
  ## Logistic Regression---------------------------------------------------------------
  model_lr <- glm(Nociplastic ~ CSI, data = train, family = binomial)
  predict_prob <- predict(model_lr, test, type = "response")
  predict_class <- as.factor(ifelse(predict_prob > 0.5, 1, 0))
  lr_acc[i] <- mean(test$Nociplastic == predict_class)
  
  ## Random Forest
  control <- trainControl(
    method = "cv",      # Cross-validation
    number = 5,         # Number of folds
    verboseIter = FALSE  # Print progress
  )
  
  # Grid for tuning
  tune_grid <- expand.grid(
    mtry = c(1, 2, 3, 4)  # Number of predictors at each split
  )
  
  # Train Random Forest with caret
  rf_model <- train(
    Nociplastic ~ CSI,
    data = train,
    method = "rf",  # Random Forest
    trControl = control,
    tuneGrid = tune_grid,
    ntree = 500  # Number of trees
  )
  
  predict_prob <- predict(rf_model, test, type = "prob")[,2]
  predict_class <- as.factor(ifelse(predict_prob > 0.5, 1, 0))
  rf_acc[i] <- mean(test$Nociplastic == predict_class)
  
  ## SVM
  model_svm <- train(Nociplastic ~ CSI , data = train, 
                     method = "svmLinear")
  svml_acc[i] <- mean(predict(model_svm, test) == test$Nociplastic)
  
  # Train the model with a radial basis function (RBF) kernel
  svm_model_radial <- train(
    Nociplastic ~ CSI,
    data = train,
    method = "svmRadial",  # Radial kernel
    trControl = trainControl(method = "cv", number = 5),  # 5-fold cross-validation
    tuneLength = 5  # Number of tuning parameter combinations to try
  )
  
  # View the model results
  svmr_acc[i] <- mean(predict(svm_model_radial, test) == test$Nociplastic)
}

mean(lr_acc)
mean(rf_acc)
mean(svml_acc)
mean(svmr_acc)

# Find the best prediction model --------------------------------------------
n = 100
lr_acc <- numeric(n)
rf_acc <- numeric(n)
svml_acc <- numeric(n)
svmr_acc <- numeric(n)

for (i in 1:n){
  ## Split train and test
  set.seed(i)
  trn_idx <- sample(1:nrow(df_pre), 0.7*nrow(df_pre))
  train <- df_pre[trn_idx,]
  test <- df_pre[-trn_idx,]
  
  ## Logistic Regression---------------------------------------------------------------
  model_lr <- glm(Nociplastic ~ CSI + MPD, data = train, family = binomial)
  predict_prob <- predict(model_lr, test, type = "response")
  predict_class <- as.factor(ifelse(predict_prob > 0.5, 1, 0))
  lr_acc[i] <- mean(test$Nociplastic == predict_class)
  
  ## Random Forest
  control <- trainControl(
    method = "cv",      # Cross-validation
    number = 5,         # Number of folds
    verboseIter = FALSE  # Print progress
  )
  
  # Grid for tuning
  tune_grid <- expand.grid(
    mtry = c(1, 2, 3, 4)  # Number of predictors at each split
  )
  
  # Train Random Forest with caret
  rf_model <- train(
    Nociplastic ~ CSI + MPD,
    data = train,
    method = "rf",  # Random Forest
    trControl = control,
    tuneGrid = tune_grid,
    ntree = 500  # Number of trees
  )
  
  predict_prob <- predict(rf_model, test, type = "prob")[,2]
  predict_class <- as.factor(ifelse(predict_prob > 0.5, 1, 0))
  rf_acc[i] <- mean(test$Nociplastic == predict_class)
  
  ## SVM
  model_svm <- train(Nociplastic ~ CSI + MPD , data = train, 
                     method = "svmLinear")
  svml_acc[i] <- mean(predict(model_svm, test) == test$Nociplastic)
  
  # Train the model with a radial basis function (RBF) kernel
  svm_model_radial <- train(
    Nociplastic ~ CSI + MPD,
    data = train,
    method = "svmRadial",  # Radial kernel
    trControl = trainControl(method = "cv", number = 5),  # 5-fold cross-validation
    tuneLength = 5  # Number of tuning parameter combinations to try
  )
  
  # View the model results
  svmr_acc[i] <- mean(predict(svm_model_radial, test) == test$Nociplastic)
}

mean(lr_acc)
mean(rf_acc)
mean(svml_acc)
mean(svmr_acc)

# Selecting the best cut of CSI -------------------------------------------------------------
CSI_cut <- function(x, train_data = df_pre){
  train_data$MPD <- ifelse(train_data$MPD > x, 1, 0)
  train_control <- trainControl(method = "repeatedcv", 
                                number = 10, repeats = 10)
  model <- train(Nociplastic ~ CSI + MPD, data = train_data, 
                 method = "svmLinear",
                 trControl = train_control)
  return(model$results$Accuracy[1])
}

x_values <- log(seq(1, 30, by = 1))
y_values <- c()
for (x in x_values){
  y_values <- c(y_values, CSI_cut(x))
}
data <- data.frame(x = x_values, y = y_values)

# Plot the function using ggplot2
ggplot(data, aes(x = exp(x)-1, y = y)) +
  geom_line(color = "blue", size = 1) +
  labs(title = "CV Accuracy VS MPD Cutoff", x = "MPD Cutoff", y = "Mean Accuracy") +
  theme_minimal()

train_data <- df_pre
train_data$MPD <- ifelse(train_data$MPD > log(7), 1, 0)
train_control <- trainControl(method = "repeatedcv", 
                              number = 10, repeats = 10)
model <- train(Nociplastic ~ CSI + MPD, data = train_data, 
               method = "svmLinear",
               trControl = train_control)
model$results$AccuracySD

# Final Model-------------------------------------------------------------------
model_svm <- svm(Nociplastic ~ CSI + MPD, data = df_pre, kernel = "linear", scale = F)
plot(model_svm, df_pre[,3:5])

x_seq <- seq(min(df_pre$CSI), max(df_pre$CSI), length.out = 100)
y_seq <- seq(min(df_pre$MPD), max(df_pre$MPD), length.out = 100)
prediction_grid <- expand.grid(CSI = x_seq, MPD = y_seq)

# Predict classes for the grid
prediction_grid$predicted_class <- predict(model_svm, newdata = prediction_grid)
mean(predict(model_svm, newdata = df_pre) == df_pre$Nociplastic)

ggplot(df_pre, aes(x = CSI, y = MPD, color = Nociplastic)) +
  geom_point() +  # Plot original data points
  stat_contour(data = prediction_grid, aes(z = as.numeric(predicted_class)), color = "black") +
  labs(title = "SVM Decision Boundary") + ylab('log (MPD+1)')

# Extract support vectors and their coefficients
support_vectors <- model_svm$SV
coefficients <- model_svm$coefs
intercept <- model_svm$rho  # Note the negative sign for `rho`

# Compute the weights (w)
weights <- t(coefficients) %*% as.matrix(support_vectors)
print(weights)
print(intercept)

## --------------
# Calculate centroids for each class
centroids <- df_pre %>%
  group_by(Nociplastic) %>%
  summarize(
    mean_CSI = mean(CSI),
    mean_log_MPD = mean(log(MPD + 1))
  )

# Create the plot with centroids
ggplot(df_pre, aes(x = CSI, y = MPD, color = Nociplastic)) +
  geom_point() +  # Plot original data points
  stat_contour(data = prediction_grid, aes(x = CSI, y = MPD, z = as.numeric(predicted_class)), color = "black") +
  geom_point(
    data = centroids, aes(x = mean_CSI, y = mean_log_MPD),
    color = "purple", size = 6, shape = 8  # Red star-shaped points for centroids
  ) +
  labs(
    title = "SVM Decision Boundary with Centroids",
    x = "CSI",
    y = "log(MPD + 1)"
  ) +
  theme_minimal()

# Selecting the best cut of CSI (Using 70/30 split rather than CV) -------------------------------------------------------------
CSI_cut <- function(x, train_data = df_pre){
  svml_acc <- numeric(100)
  train_data$CSI <- ifelse(train_data$CSI > x, 1, 0)
  for (i in 1:100){
    set.seed(i)
    trn_idx <- sample(1:nrow(train_data), 0.7*nrow(train_data))
    train <- train_data[trn_idx,]
    test <- train_data[-trn_idx,]
    model_svm <- train(Nociplastic ~ CSI + MPD , data = train, 
                       method = "svmLinear")
    svml_acc[i] <- mean(predict(model_svm, test) == test$Nociplastic)
  }
  
  return(mean(svml_acc))
}

x_values <- seq(30, 60, by = 1)
y_values <- c()
for (x in x_values){
  y_values <- c(y_values, CSI_cut(x))
}
data <- data.frame(x = x_values, y = y_values)

# Plot the function using ggplot2
ggplot(data, aes(x = x, y = y)) +
  geom_line(color = "blue", size = 1) +
  labs(title = "Mean Prediction Accuracy VS CSI Cutoff", x = "CSI Cutoff", y = "Mean Accuracy") +
  theme_minimal()

#---------------------------------------------------------
library(caret)
library(ggplot2)
library(reshape2)

CSI_MPD_analysis <- function(CSI_range, MPD_range, train_data){
  
  if(!all(c("CSI", "Nociplastic") %in% names(train_data))){
    stop("Error: Required columns ('CSI', 'Nociplastic') are missing in train_data.")
  }
  
  # Create an empty matrix to store accuracies
  accuracy_matrix <- matrix(NA, nrow = length(CSI_range), ncol = length(MPD_range),
                            dimnames = list(CSI_range, MPD_range))
  
  for (x in CSI_range) {
    for (y in MPD_range) {
      
      svml_acc <- numeric(10)
      
      # Define binary MPD and CSI thresholds
      train_data$CSI_bin <- ifelse(train_data$CSI > x, 1, 0)
      train_data$MPD_bin <- ifelse(train_data$MPD > log(y+1), 1, 0)
      
      for (i in 1:10) {
        set.seed(i)
        
        trn_idx <- sample(1:nrow(train_data), 0.7 * nrow(train_data))
        train <- train_data[trn_idx, ]
        test <- train_data[-trn_idx, ]
        
        # Train SVM model with both CSI and MPD as predictors
        model_svm <- train(Nociplastic ~ CSI_bin + MPD_bin, data = train, method = "svmLinear")
        
        # Store accuracy
        svml_acc[i] <- mean(predict(model_svm, test) == test$Nociplastic)
      }
      
      # Store the mean accuracy for (CSI threshold = x, MPD threshold = y)
      accuracy_matrix[as.character(x), as.character(y)] <- mean(svml_acc)
    }
  }
  
  # Convert matrix to dataframe for visualization
  accuracy_df <- melt(accuracy_matrix, varnames = c("CSI_threshold", "MPD_threshold"),
                      value.name = "Accuracy")
  
  return(accuracy_df)
}

# Example usage:
CSI_vals <- seq(30, 60, by = 1)  # Define range of CSI thresholds
MPD_vals <- seq(1, 30, by = 1)  # Define range of MPD thresholds

accuracy_df <- CSI_MPD_analysis(CSI_vals, MPD_vals, train_data = df_pre)
write.csv(accuracy_df, 'accuracy_2d.csv')

# Plot heatmap
ggplot(accuracy_df, aes(x = as.numeric(CSI_threshold), 
                        y = as.numeric(MPD_threshold), 
                        fill = Accuracy)) +
  geom_tile(color = "white") +  # Adds grid lines for better readability
  scale_fill_viridis_c(option = "magma", direction = -1) +  # Improved color scale
  labs(title = "SVM Accuracy for Different CSI & MPD Cutoffs",
       x = "CSI Cutoff",
       y = "MPD Cutoff",
       fill = "Accuracy") +
  theme_minimal(base_size = 14) +  # Slightly larger font for readability
  theme(panel.grid = element_blank(),  # Remove background grid
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

#-------------------------------
df_pre$CSI_bin <- ifelse(df_pre$CSI > 46, 1, 0)
df_pre$MPD_bin <- ifelse(df_pre$MPD > log(20+1), 1, 0)
xtabs(~ Nociplastic + CSI_bin + MPD_bin, data = df_pre)

library(janitor)

df_pre %>%
  tabyl(CSI_bin, MPD_bin, Nociplastic) %>%
  adorn_totals()

# Sample Size Estimates ---------------------------------
library(e1071)  # SVM package
library(caret)

# Simulated data
set.seed(123)
sample_sizes <- seq(50, 1000, by = 50)  # Varying sample sizes
accuracies <- c()

for (n in sample_sizes) {
  # Generate synthetic dataset
  train_data <- twoClassSim(n = n, linearVars = 2)
  
  # Train SVM
  model <- train(Class ~ ., data = train_data, method = "svmLinear",
                 trControl = trainControl(method = "cv", number = 5))
  
  # Store accuracy
  accuracies <- c(accuracies, mean(model$results$Accuracy))
}

# Plot Learning Curve
plot(sample_sizes, accuracies, type = "b", pch = 19,
     xlab = "Sample Size", ylab = "Accuracy",
     main = "Learning Curve for SVM")

# Circos Plot ----------------------------------------------------------------------
library(circlize)

df <- read.csv("df_pre_orig.csv")

# Remove the first column if it's an index
df <- df[ , -1]

# Compute the correlation matrix
cor_matrix <- cor(df, use = "complete.obs")

# Convert correlation matrix to long format for Circos plotting
cor_long <- as.data.frame(as.table(cor_matrix))

# Remove self-correlations (diagonal)
cor_long <- cor_long[cor_long$Var1 != cor_long$Var2, ]

# Set up the Circos plot
circos.clear()
circos.par(gap.after = 10)  # Adjust spacing between sectors

# Define unique sectors (column names)
sectors <- unique(c(cor_long$Var1, cor_long$Var2))

# Initialize the Circos plot
circos.initialize(factors = sectors, xlim = c(0, 1))

# Add sectors (labels)
circos.track(ylim = c(0, 1), panel.fun = function(x, y) {
  sector_name <- get.cell.meta.data("sector.index")
  circos.text(CELL_META$xcenter, CELL_META$ycenter, sector_name, facing = "inside", niceFacing = TRUE)
})

# Define color scale for correlation values
color_scale <- colorRampPalette(c("blue", "white", "red"))(100)

# Add correlation links
for (i in 1:nrow(cor_long)) {
  from <- as.character(cor_long$Var1[i])
  to <- as.character(cor_long$Var2[i])
  value <- cor_long$Freq[i]
  
  # Assign color based on correlation value
  color_index <- round((value + 1) * 49) + 1  # Scale [-1,1] to [1,100]
  link_color <- color_scale[color_index]
  
  circos.link(from, 0.5, to, 0.5, col = link_color, border = NA, lwd = 3 * abs(value))
}

# Finalize the Circos plot
circos.clear()



