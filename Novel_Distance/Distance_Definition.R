# Function for Calculating the Co-occurrence Distance
cooccur1 <- function (data, cat_col) 
{
  require(kmed)
  # cat_col: index of categorical column
  if ((is.matrix(data) || is.data.frame(data)) == FALSE) 
    stop("The data must be a matrix or a data frame!")
  rn <- rownames(data)
  if (is.logical(data) == TRUE) 
    data <- data.matrix(as.data.frame(data)) + as.integer(1)
  if (is.numeric(data) == TRUE) 
    data <- apply(data, 2, function(x) as.integer(x))
  if (is.data.frame(data) == TRUE) 
    data <- apply(data, 2, function(x) as.integer(as.factor(x)))
  col <- ncol(data)
  ncatcol <- length(cat_col) #
  if (col != 1) {
    newdist <- function(data, col = col, colnum) {
      nvar <- 1:col
      n <- length(levels(as.factor(data[, colnum])))
      var <- prob <- vector("list", (col - 1))
      for (i in 1:(col - 1)) {
        var[[i]] <- table(data[, nvar[-colnum][i]], data[, 
                                                         colnum])
        prob[[i]] <- var[[i]]/matrix(colSums(var[[i]]), 
                                     nrow = nrow(var[[i]]), ncol = ncol(var[[i]]), 
                                     byrow = TRUE)
      }
      probmat <- do.call(rbind, prob)
      matnew <- matrix(0, nrow = n, ncol = n)
      rownames(matnew) <- colnames(matnew) <- 1:n
      for (i in 1:n) {
        for (j in 1:n) {
          matnew[i, j] <- (sum(apply(probmat[, c(i, j)], 
                                     1, max)) - (col - 1))/(col - 1)
        }
      }
      return(matnew)
    }
    newdata <- vector("list", ncatcol)
    for (i in 1:ncatcol) {
      # Only calculate distance for categorical variables but also consider other numeric variables
      newdata[[i]] <- newdist(data, col = col, cat_col[i])
    }
    distmat <- matrix(0, nrow(data), nrow(data))
    for (i in 1:nrow(data)) {
      for (j in 1:nrow(data)) {
        distsum <- numeric(ncatcol) #
        for (k in 1:ncatcol) {#
          distsum[k] <- newdata[[k]][data[i, cat_col[k]], data[j, cat_col[k]]]
        }
        distmat[i, j] <- sum(distsum)
      }
    }
    rownames(distmat) <- colnames(distmat) <- rn
  }
  else {
    distmat <- kmed::matching(data, data)
    rownames(distmat) <- colnames(distmat) <- rn
    warning("Due to only 1 variable, simple matching distance is calculated instead!\n            To produce coocurrence distance, it requires at least 2 variables.")
  }
  return(distmat)
}

# The function for calculating the proposed distance for questionnaire data
ques_dist <- function(df, num_col = NULL, cat_col = NULL, sr_col = NULL){
  # df: The dataframe with the data to cluster
  # num_col: The index of numeric and ordinal (tranformed to ranks) columns
  # cat_col: The index of categorical columns
  # ord_col: The index of ordinal columns
  # sr_col: The list of self-reported (subjective) columns
  require(StatMatch)
  require(arules)
  require(kmed)
  dist_num = 0
  dist_cat = 0
  dist_sr = 0
  n1 <- length(num_col)
  n2 <- length(cat_col)
  n3 <- length(sr_col)
  n_sr <- 0
  
  # Calculate the Gower's distance for numeric variables
  if (n1>0){
    df_num <- scale(df[, num_col])
    dist_num <- gower.dist(df_num) 
  }
  
  # Calculate the Co-occurrence distance for categorical variables
  if (n2>0){
    df_discrete <- discretizeDF(df, default = list(method = "interval")) # Discretize the dataframe
    dist_cat <- cooccur1(df_discrete, cat_col)
  }
  
  if (n3>0){
    for (i in 1:n3){
      df_sr <- df[, sr_col[[i]]]
      
      # Calculate the Correlation Distance for self-reported variables
      dist_sr1 <- sin(acos(cor(t(df_sr)))/2)
      na_index <- which(is.na(dist_sr1), arr.ind=TRUE)
      if (nrow(na_index) > 0){
        for (k in 1:nrow(na_index)){
          if (sd(df_sr[na_index[k, 1],]) == 0 & (sd(df_sr[na_index[k,2],]) == 0)){
            dist_sr1[na_index[k,1], na_index[k,2]] <- 0
          }
          else{
            dist_sr1[na_index[k,1], na_index[k,2]] <- sqrt(2)/2
          }  
        }
      }
      
      #dist_sr1[is.na(dist_sr1)] <- sqrt(2)/2
      dist_sr1_lt <- dist_sr1[lower.tri(dist_sr1, diag = FALSE)]
      w = sd(dist_sr1_lt)
      
      # Calculate the Gower Distance for self-reported variables
      dist_sr2 <- gower.dist(df_sr)
      
      dist_sr <- dist_sr + w * dist_sr1 + dist_sr2 * length(sr_col[[i]])
      
      n_sr <- n_sr + length(sr_col[[i]])
    }
  }
  
  all_dist <- (dist_num * n1 + dist_cat + dist_sr)/(n1 + n2 + n_sr + n3)
  return (all_dist)
}

