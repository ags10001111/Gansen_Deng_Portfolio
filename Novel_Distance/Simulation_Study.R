source("Distance_Definition.R")
library(MASS)
library(kamila)
library(mixdist)
#------------------------------------------------------------------
# Simulation study
## The function for calculating all distances for comparison
set.seed(916)
distmat <- function(df, num_col = NULL, cat_col = NULL, sr_col = NULL, sr_colp = NULL,w = 0.5, method = "ques"){
  require(kmed)
  if (method == "ques"){
    all_dist <- ques_dist(df, num_col, cat_col, sr_col)
  }
  else if (method == "ahmad"){
    all_dist <- sqrt(distmix(df, method = method, idnum = c(num_col, unlist(sr_col)), idcat = cat_col))
  }
  else {
    all_dist <- distmix(df, method = method, idnum = c(num_col, unlist(sr_col)), idcat = cat_col)
  }
  return (all_dist)
}

## The function for calculating ARI with simulated data
ARI_sim <- function(data, trueID, num_col = NULL, cat_col = NULL, sr_col = NULL){
  require(mclust)
  require(cluster)
  require(clValid)
  require(clusterSim)
  ARI_res <- c()
  CI_res <- c()
  SS_res <- c()
  method_vec <- c("ques", "gower", "podani", "ahmad", "harikumar")
  for (method0 in method_vec){
    print(method0)
    all_dist <- distmat(data, num_col, cat_col, sr_col, method = method0)
    hclust_avg <- hclust(as.dist(all_dist), method = 'average')
    cut_avg <- cutree(hclust_avg, k = 2)
    ARI0 <- adjustedRandIndex(trueID, cut_avg) # ARI
    CI0 <- index.C(as.dist(all_dist), cut_avg) # C-Index
    SS0 <- summary(silhouette(cut_avg, dmatrix = all_dist))$avg.width # Silhouette Score
    ARI_res <- append(ARI_res, ARI0)
    CI_res <- append(CI_res, CI0)
    SS_res <- append(SS_res, SS0)
  }
  names(ARI_res) <- method_vec
  names(CI_res) <- method_vec
  names(SS_res) <- method_vec
  return(list(ARI_res, CI_res, SS_res))
}

cor_vec <- function(x1, rho){
  n <- length(x1)
  theta <- acos(rho)
  x2 <- rnorm(n, 2, 0.5)
  
  X     <- cbind(x1, x2)         # matrix
  Xctr  <- scale(X, center=TRUE, scale=FALSE)   # centered columns (mean 0)
  
  Id   <- diag(n)                               # identity matrix
  Q    <- qr.Q(qr(Xctr[ , 1, drop=FALSE]))      # QR-decomposition, just matrix Q
  P    <- tcrossprod(Q)          # = Q Q'       # projection onto space defined by x1
  x2o  <- (Id-P) %*% Xctr[ , 2]                 # x2ctr made orthogonal to x1ctr
  Xc2  <- cbind(Xctr[ , 1], x2o)                # bind to matrix
  Y    <- Xc2 %*% diag(1/sqrt(colSums(Xc2^2)))  # scale columns to length 1
  
  x <- Y[ , 2] + (1 / tan(theta)) * Y[ , 1]     # final new vector
  return (x)
}

# Simulation Idea 1------------------------
sim0 <- function(SR_prop, b, rho){
  ARI_df = data.frame()
  for(i in 1:30){
    if (SR_prop == 20){
      x1 <- rnorm(2,0,1)
      x2 <- cor_vec(x1, rho)
      dat <- genMixedData(500, 8, 4, nCatLevels=4, nConWithErr=8, nCatWithErr=4,
                          popProportions=c(0.3,0.7), conErrLev=0.1, catErrLev=0.1)
      sim_data_A <- as.data.frame(cbind(dat$conVars, dat$catVars))
    }
    if (SR_prop == 60){
      x11 <- rnorm(3,0,1)
      x21 <- cor_vec(x11, rho)
      
      x12 <-  rnorm(3, 0, 1)
      x22 <- cor_vec(x12, rho)
      
      x1 <- c(x11, x12)
      x2 <- c(x21, x22)
      dat <- genMixedData(500, 4, 4, nCatLevels=4, nConWithErr=4, nCatWithErr=4,
                          popProportions=c(0.3,0.7), conErrLev=0.1, catErrLev=0.1)
      sim_data_A <- as.data.frame(cbind(dat$conVars, dat$catVars))
    }
    if (SR_prop == 100){
      x11 <- rnorm(3,0,1)
      x21 <- cor_vec(x11, rho)
      
      x12 <-  rnorm(3, 0, 1)
      x22 <- cor_vec(x12, rho)
      
      x13 <-  rnorm(4, 0, 1)
      x23 <- cor_vec(x13, rho)
      
      x1 <- c(x11, x12, x13)
      x2 <- c(x21, x22, x23)
      dat <- genMixedData(500, 1, 4, nCatLevels=4, nConWithErr=1, nCatWithErr=4,
                          popProportions=c(0.3,0.7), conErrLev=0.1, catErrLev=0.1)
      
      sim_data_A <- dat$catVars
    }
    means <- list(
      x1,   # Cluster 1 mean
      x2   # Cluster 2 mean
    )
    n = length(x1)
    G1 <- matrix(runif(n^2)*2-1, ncol=n) 
    G2 <- matrix(runif(n^2)*2-1, ncol=n) 
    cov_matrices <- list(
      t(G1) %*% G1,  # Cluster 1 covariance
      t(G2) %*% G2   # Cluster 2 covariance
    )
    sim_data_SR<- matrix(NA, nrow = 500, ncol = n)
    for (k in 1:2) {
      cluster_indices <- which(dat$trueID == k)
      n_k <- length(cluster_indices)
      # Simulate data from multivariate normal distribution for each cluster
      sim_data_SR[cluster_indices, ] <- mvrnorm(n_k, means[[k]], cov_matrices[[k]])
    }
    sim_data <- cbind(sim_data_SR, sim_data_A)
    if (SR_prop == 20){
      sds <- apply(sim_data[,1:2], 2, sd)
      e <- rnorm(500, 0, 1)
      e <- sapply(sds, function(s){e*(s/b)})
      sim_data[,1:2] <- sim_data[, 1:2] + e
      ARI_res <- ARI_sim(sim_data, dat$trueID, 3:10, 11:14, 1:2)
    }
    if (SR_prop == 60){
      sds_1 <- apply(sim_data[,1:3], 2, sd)
      sds_2 <- apply(sim_data[,4:6], 2, sd)
      e1 <- rnorm(500, 0, 1)
      e2 <- rnorm(500, 0, 1)
      e1 <- sapply(sds_1, function(s){e1*(s/b)})
      e2 <- sapply(sds_2, function(s){e2*(s/b)})
      sim_data[,1:3] <- sim_data[, 1:3] + e1
      sim_data[,4:6] <- sim_data[, 4:6] + e2
      ARI_res <- ARI_sim(sim_data, dat$trueID, 7:10, 11:14, list(c(1:3), c(4:6)))
    }
    if (SR_prop == 100){
      sds_1 <- apply(sim_data[,1:3], 2, sd)
      sds_2 <- apply(sim_data[,4:6], 2, sd)
      sds_3 <- apply(sim_data[,7:10], 2, sd)
      e1 <- rnorm(500, 0, 1)
      e2 <- rnorm(500, 0, 1)
      e3 <- rnorm(500, 0, 1)
      e1 <- sapply(sds_1, function(s){e1*(s/b)})
      e2 <- sapply(sds_2, function(s){e2*(s/b)})
      e3 <- sapply(sds_3, function(s){e3*(s/b)})
      sim_data[,1:3] <- sim_data[, 1:3] + e1
      sim_data[,4:6] <- sim_data[, 4:6] + e2
      sim_data[,7:10] <- sim_data[, 7:10] + e3
      ARI_res <- ARI_sim(sim_data, dat$trueID, NULL, 11:14, list(c(1:3), c(4:6), c(7:10)))
    }
    ARI_df <- rbind(ARI_df, ARI_res)
  }
  colnames(ARI_df) <- c("ques", "gower", "podani", "ahmad", "harikumar")
  Mean <- round(colMeans(ARI_df),3)
  SD <- round(apply(ARI_df,2,sd),3)
  Mean_SD <- paste0(Mean, ' (', SD, ')')
  ARI_df <- rbind(ARI_df, Mean, SD, Mean_SD)
  rownames(ARI_df)[31:33] <- c('Mean', 'SD', 'Mean_SD')
  file_name <- paste('SR ', SR_prop, '_Bias ', b, '_Rho ', rho, sep = "")
  file_name <- paste(file_name, '.csv', sep = "")
  write.csv(ARI_df, file_name)
}
sim0(100, 2, -0.6)

SR_props <- c(20, 60, 100)
bs <- c(1, 2, 3)
rhos <- c(0.6, 0.3, 0, -0.3, -0.6)

para_grid <- expand.grid(SR_prop = SR_props, b = bs,
                     rho = rhos)
with(para_grid, sapply(1:nrow(para_grid), function(j){sim0(SR_props[j], bs[j], rhos[j])}))

# Simulation Idea 2------------------------
## Normal Error
sim0 <- function(SR_prop, b, ol){
  ARI_df = data.frame()
  CI_df = data.frame()
  SS_df = data.frame()
  for(i in 1:30){
    dat <- genMixedData(500, 10, 4, nCatLevels=4, nConWithErr=10, nCatWithErr=4,
                        popProportions=c(0.3,0.7), conErrLev=ol, catErrLev=ol)
    sim_data <- as.data.frame(cbind(dat$conVars, dat$catVars))
    if (SR_prop == 20){
      sds <- apply(sim_data[,1:2], 2, sd)
      e <- rnorm(500, 0, 1)
      e <- sapply(sds, function(s){e*(s/b)})
      sim_data[,1:2] <- sim_data[, 1:2] + e
      ARI_res <- ARI_sim(sim_data, dat$trueID, 3:10, 11:14, list(1:2))[[1]]
      CI_res <- ARI_sim(sim_data, dat$trueID, 3:10, 11:14, list(1:2))[[2]]
      SS_res <- ARI_sim(sim_data, dat$trueID, 3:10, 11:14, list(1:2))[[3]]
    }
    if (SR_prop == 60){
      sds_1 <- apply(sim_data[,1:3], 2, sd)
      sds_2 <- apply(sim_data[,4:6], 2, sd)
      e1 <- rnorm(500, 0, 1)
      e2 <- rnorm(500, 0, 1)
      e1 <- sapply(sds_1, function(s){e1*(s/b)})
      e2 <- sapply(sds_2, function(s){e2*(s/b)})
      sim_data[,1:3] <- sim_data[, 1:3] + e1
      sim_data[,4:6] <- sim_data[, 4:6] + e2
      ARI_res <- ARI_sim(sim_data, dat$trueID, 7:10, 11:14, list(c(1:3), c(4:6)))[[1]]
      CI_res <- ARI_sim(sim_data, dat$trueID, 7:10, 11:14, list(c(1:3), c(4:6)))[[2]]
      SS_res <- ARI_sim(sim_data, dat$trueID, 7:10, 11:14, list(c(1:3), c(4:6)))[[3]]
    }
    if (SR_prop == 100){
      sds_1 <- apply(sim_data[,1:3], 2, sd)
      sds_2 <- apply(sim_data[,4:6], 2, sd)
      sds_3 <- apply(sim_data[,7:10], 2, sd)
      e1 <- rnorm(500, 0, 1)
      e2 <- rnorm(500, 0, 1)
      e3 <- rnorm(500, 0, 1)
      e1 <- sapply(sds_1, function(s){e1*(s/b)})
      e2 <- sapply(sds_2, function(s){e2*(s/b)})
      e3 <- sapply(sds_3, function(s){e3*(s/b)})
      sim_data[,1:3] <- sim_data[, 1:3] + e1
      sim_data[,4:6] <- sim_data[, 4:6] + e2
      sim_data[,7:10] <- sim_data[, 7:10] + e3
      ARI_res <- ARI_sim(sim_data, dat$trueID, NULL, 11:14, list(c(1:3), c(4:6), c(7:10)))[[1]]
      CI_res <- ARI_sim(sim_data, dat$trueID, NULL, 11:14, list(c(1:3), c(4:6), c(7:10)))[[2]]
      SS_res <- ARI_sim(sim_data, dat$trueID, NULL, 11:14, list(c(1:3), c(4:6), c(7:10)))[[3]]
    }
    ARI_df <- rbind(ARI_df, ARI_res)
    CI_df <- rbind(CI_df, CI_res)
    SS_df <- rbind(SS_df, SS_res)
  }
  colnames(ARI_df) <- c("ques", "gower", "podani", "ahmad", "harikumar")
  Mean <- round(colMeans(ARI_df),3)
  SD <- round(apply(ARI_df,2,sd),3)
  Mean_SD <- paste0(Mean, ' (', SD, ')')
  ARI_df <- rbind(ARI_df, Mean, SD, Mean_SD)
  rownames(ARI_df)[31:33] <- c('Mean', 'SD', 'Mean_SD')
  file_name <- paste('ARI_', 'SR ', SR_prop, '_Bias ', b, '_OL ', ol, sep = "")
  file_name <- paste(file_name, '.csv', sep = "")
  write.csv(ARI_df, file_name)
  
  colnames(CI_df) <- c("ques", "gower", "podani", "ahmad", "harikumar")
  Mean <- round(colMeans(CI_df),3)
  SD <- round(apply(CI_df,2,sd),3)
  Mean_SD <- paste0(Mean, ' (', SD, ')')
  CI_df <- rbind(CI_df, Mean, SD, Mean_SD)
  rownames(CI_df)[31:33] <- c('Mean', 'SD', 'Mean_SD')
  file_name <- paste('CI_','SR ', SR_prop, '_Bias ', b, '_OL ', ol, sep = "")
  file_name <- paste(file_name, '.csv', sep = "")
  write.csv(CI_df, file_name)
  
  colnames(SS_df) <- c("ques", "gower", "podani", "ahmad", "harikumar")
  Mean <- round(colMeans(SS_df),3)
  SD <- round(apply(SS_df,2,sd),3)
  Mean_SD <- paste0(Mean, ' (', SD, ')')
  SS_df <- rbind(SS_df, Mean, SD, Mean_SD)
  rownames(SS_df)[31:33] <- c('Mean', 'SD', 'Mean_SD')
  file_name <- paste('SS_','SR ', SR_prop, '_Bias ', b, '_OL ', ol, sep = "")
  file_name <- paste(file_name, '.csv', sep = "")
  write.csv(SS_df, file_name)
}

#sim0(100,2,0.2)

SR_props <- c(20)
bs <- c(2, 3, 6)
ols <- c(0.1, 0.2, 0.3, 0.4)

para_grid <- expand.grid(SR_prop = SR_props, b = bs, ol = ols)
with(para_grid, sapply(1:nrow(para_grid), function(j){sim0(SR_prop[j], b[j], ol[j])}))

## Weibull Errors
sim1 <- function(SR_prop, para, ol){
  ARI_df = data.frame()
  CI_df = data.frame()
  SS_df = data.frame()
  for(i in 1:30){
    dat <- genMixedData(500, 10, 4, nCatLevels=4, nConWithErr=10, nCatWithErr=4,
                        popProportions=c(0.3,0.7), conErrLev=ol, catErrLev=ol)
    sim_data <- as.data.frame(cbind(dat$conVars, dat$catVars))
    if (SR_prop == 20){
      sds <- apply(sim_data[,1:2], 2, sd)
      e <- rweibull(500, para[1], para[2]) - 3
      sim_data[,1:2] <- sim_data[, 1:2] + e
      ARI_res <- ARI_sim(sim_data, dat$trueID, 3:10, 11:14, list(1:2))[[1]]
      CI_res <- ARI_sim(sim_data, dat$trueID, 3:10, 11:14, list(1:2))[[2]]
      SS_res <- ARI_sim(sim_data, dat$trueID, 3:10, 11:14, list(1:2))[[3]]
    }
    if (SR_prop == 60){
      sds_1 <- apply(sim_data[,1:3], 2, sd)
      sds_2 <- apply(sim_data[,4:6], 2, sd)
      e1 <- rweibull(500, para[1], para[2]) - 3
      e2 <- rweibull(500, para[1], para[2]) - 3
      sim_data[,1:3] <- sim_data[, 1:3] + e1
      sim_data[,4:6] <- sim_data[, 4:6] + e2
      ARI_res <- ARI_sim(sim_data, dat$trueID, 7:10, 11:14, list(c(1:3), c(4:6)))[[1]]
      CI_res <- ARI_sim(sim_data, dat$trueID, 7:10, 11:14, list(c(1:3), c(4:6)))[[2]]
      SS_res <- ARI_sim(sim_data, dat$trueID, 7:10, 11:14, list(c(1:3), c(4:6)))[[3]]
    }
    if (SR_prop == 100){
      sds_1 <- apply(sim_data[,1:3], 2, sd)
      sds_2 <- apply(sim_data[,4:6], 2, sd)
      sds_3 <- apply(sim_data[,7:10], 2, sd)
      e1 <- rweibull(500, para[1], para[2]) - 3
      e2 <- rweibull(500, para[1], para[2]) - 3
      e3 <- rweibull(500, para[1], para[2]) - 3
      sim_data[,1:3] <- sim_data[, 1:3] + e1
      sim_data[,4:6] <- sim_data[, 4:6] + e2
      sim_data[,7:10] <- sim_data[, 7:10] + e3
      ARI_res <- ARI_sim(sim_data, dat$trueID, NULL, 11:14, list(c(1:3), c(4:6), c(7:10)))[[1]]
      CI_res <- ARI_sim(sim_data, dat$trueID, NULL, 11:14, list(c(1:3), c(4:6), c(7:10)))[[2]]
      SS_res <- ARI_sim(sim_data, dat$trueID, NULL, 11:14, list(c(1:3), c(4:6), c(7:10)))[[3]]
    }
    ARI_df <- rbind(ARI_df, ARI_res)
    CI_df <- rbind(CI_df, CI_res)
    SS_df <- rbind(SS_df, SS_res)
  }
  colnames(ARI_df) <- c("ques", "gower", "podani", "ahmad", "harikumar")
  Mean <- round(colMeans(ARI_df),3)
  SD <- round(apply(ARI_df,2,sd),3)
  Mean_SD <- paste0(Mean, ' (', SD, ')')
  ARI_df <- rbind(ARI_df, Mean, SD, Mean_SD)
  rownames(ARI_df)[31:33] <- c('Mean', 'SD', 'Mean_SD')
  file_name <- paste('ARI_', 'SR ', SR_prop, '_Para ', para[1], " ", para[2], '_OL ', ol, sep = "")
  file_name <- paste(file_name, '.csv', sep = "")
  write.csv(ARI_df, file_name)
  
  colnames(CI_df) <- c("ques", "gower", "podani", "ahmad", "harikumar")
  Mean <- round(colMeans(CI_df),3)
  SD <- round(apply(CI_df,2,sd),3)
  Mean_SD <- paste0(Mean, ' (', SD, ')')
  CI_df <- rbind(CI_df, Mean, SD, Mean_SD)
  rownames(CI_df)[31:33] <- c('Mean', 'SD', 'Mean_SD')
  file_name <- paste('CI_','SR ', SR_prop, '_Para ', para[1], " ", para[2], '_OL ', ol, sep = "")
  file_name <- paste(file_name, '.csv', sep = "")
  write.csv(CI_df, file_name)
  
  colnames(SS_df) <- c("ques", "gower", "podani", "ahmad", "harikumar")
  Mean <- round(colMeans(SS_df),3)
  SD <- round(apply(SS_df,2,sd),3)
  Mean_SD <- paste0(Mean, ' (', SD, ')')
  SS_df <- rbind(SS_df, Mean, SD, Mean_SD)
  rownames(SS_df)[31:33] <- c('Mean', 'SD', 'Mean_SD')
  file_name <- paste('SS_','SR ', SR_prop, '_Para ', para[1], " ", para[2], '_OL ', ol, sep = "")
  file_name <- paste(file_name, '.csv', sep = "")
  write.csv(SS_df, file_name)
}

#sim1(100, c(10,6), 0.4)

SR_props <- c(20)
paras <- list(c(3,4), c(8,4), c(10,6))
ols <- c(0.1, 0.2, 0.3, 0.4)

para_grid <- expand.grid(SR_prop = SR_props, para = paras, ol = ols)
with(para_grid, sapply(1:nrow(para_grid), function(j){sim1(SR_prop[j], para[[j]], ol[j])}))

## No SR-------------------------------
sim_nsr <- function(ol){
  ARI_df = data.frame()
  CI_df = data.frame()
  SS_df = data.frame()
  for(i in 1:30){
  dat <- genMixedData(500, 10, 4, nCatLevels=4, nConWithErr=10, nCatWithErr=4,
                      popProportions=c(0.3,0.7), conErrLev=ol, catErrLev=ol)
  sim_data_A <- as.data.frame(cbind(dat$conVars, dat$catVars))
  ARI_res <- ARI_sim(sim_data_A, dat$trueID, 1:10, 11:14, NULL)[[1]]
  CI_res <- ARI_sim(sim_data_A, dat$trueID, 1:10, 11:14, NULL)[[2]]
  SS_res <- ARI_sim(sim_data_A, dat$trueID, 1:10, 11:14, NULL)[[3]]
  ARI_df <- rbind(ARI_df, ARI_res)
  CI_df <- rbind(CI_df, CI_res)
  SS_df <- rbind(SS_df, SS_res)
}
colnames(ARI_df) <- c("ques", "gower", "podani", "ahmad", "harikumar")
Mean <- round(colMeans(ARI_df),3)
SD <- round(apply(ARI_df,2,sd),3)
Mean_SD <- paste0(Mean, ' (', SD, ')')
ARI_df <- rbind(ARI_df, Mean, SD, Mean_SD)
rownames(ARI_df)[31:33] <- c('Mean', 'SD', 'Mean_SD')
file_name <- paste('ARI_NOSR_OL ', ol, '.csv', sep = "")
write.csv(ARI_df, file_name)

colnames(CI_df) <- c("ques", "gower", "podani", "ahmad", "harikumar")
Mean <- round(colMeans(CI_df),3)
SD <- round(apply(CI_df,2,sd),3)
Mean_SD <- paste0(Mean, ' (', SD, ')')
CI_df <- rbind(CI_df, Mean, SD, Mean_SD)
rownames(CI_df)[31:33] <- c('Mean', 'SD', 'Mean_SD')
file_name <- paste('CI_NOSR_OL ', ol, '.csv', sep = "")
write.csv(CI_df, file_name)

colnames(SS_df) <- c("ques", "gower", "podani", "ahmad", "harikumar")
Mean <- round(colMeans(SS_df),3)
SD <- round(apply(SS_df,2,sd),3)
Mean_SD <- paste0(Mean, ' (', SD, ')')
SS_df <- rbind(SS_df, Mean, SD, Mean_SD)
rownames(SS_df)[31:33] <- c('Mean', 'SD', 'Mean_SD')
file_name <- paste('SS_NOSR_OL ', ol, '.csv', sep = "")
write.csv(SS_df, file_name)
}
ols <- c(0.1, 0.2, 0.3, 0.4)
sapply(ols, sim_nsr)

