#############
# R functions used in Rincent et al. for evaluating the predictive abilities of SNP and NIRS through cross-validations
####
# There are 3 functions : 
# 1. RidgeBLUP.r
# 2. GBLUP.r
# 3. BL.r
# The first function (RidgeBLUP.r) was used to predict phenotypes with NIRS data using ridge regression, while the two others (GBLUP.r & BL.r) were used to predict phenotypes with SNP data using GBLUP and bayesian lasso
# Dependencies : doParallel, MASS, emma, corpcor, BGLR
# All 3 functions have the possibility to run in parallel on multicore computers and thus require the R package "doParallel"
# RidgeBLUP.r make use of the function lm.ridge from MASS package
# GBLUP.r make use of the emma.REMLE function from emma package and from two function from corpcor package
# BL_BGLR.r make use of the BGLR function from BGLR
#############
# RidgeBLUP.r
####
# Y is a vector of phenotypic records
# X is a matrix of predictors
# lambda is a vector of ridge constants
# fold is an integer defining the number of folds from the cross-validation
# iter is an integer defining the number of times the cross-validation is repeated
# cores is an integer defing the number of computer cores to use for the analysis
RidgeBLUP <- function(Y, X, lambda = 10^seq(4, -2, by = -.1), fold = 5, iter = 50, cores = 1){
  require(doParallel)
  require(MASS)
  #cores
  cl <- makeCluster(cores)
  #checks
  stopifnot(length(Y) == nrow(X))
  stopifnot(cor(match(rownames(X), names(Y)), 1:length(Y)) == 1)
  #na.omit
  data <- na.omit(data.frame(Y, X = I(scale(X))))
  n <- nrow(data)
  #fullmod
  ridgemod <- lm.ridge(Y ~ X, data = data, lambda = lambda)
  lambda_full <- ridgemod$lambda[which.min(ridgemod$GCV)]
  rm(ridgemod)
  #repeated cross-validation
  registerDoParallel(cl)
  rr_analysis <- foreach(i = 1:iter, .combine = cbind, .packages = c("MASS")) %dopar% {
    sampling <- sample(1:n)
    y_hat <- list()
    lambda_cv <- list()
    #cross-validation
    for (k in 1:fold) {
      set <- sort(sampling[floor((k - 1) * n / fold + 1):floor(k * n / fold)])
      #training
      data_cv <- data[-set, ]
      mod_ridge <- lm.ridge(Y ~ X, data = data_cv, lambda = lambda)
      #prediction
      y_hat[[k]] <- cbind(rep(1, length(set)), data$X[set, ]) %*% c(mod_ridge$ym, mod_ridge$coef[, which.min(mod_ridge$GCV)])
      rownames(y_hat[[k]]) <- rownames(data)[set]
      #storing cv lambdas
      lambda_cv[[k]] <- mod_ridge$lambda[which.min(mod_ridge$GCV)]
      rm(set, data_cv, mod_ridge)
    }
    #gather predictions of each cv repetition
    y_hat <- do.call("rbind", y_hat)
    #gather cv lambdas
    lambda_cv <- do.call("c", lambda_cv)
    return(list("y_hat" = y_hat[match(rownames(data), rownames(y_hat)), 1], "lambda" = lambda_cv))
  }
  stopCluster(cl)
  Y_pred <- do.call("cbind", rr_analysis[which(1:length(rr_analysis) %% 2 == 1)])
  lambda_cv <- do.call("c", rr_analysis[which(1:length(rr_analysis) %% 2 == 0)])
  rm(rr_analysis)
  #CV Statistics
  PRESS <- colSums((Y_pred - data$Y)^2)
  R2 <- 1 - (PRESS / (var(data$Y) * (n - 1)))
  RMSEP <- sqrt(PRESS / n)
  Rho <- cor(Y_pred, data$Y)[, 1]
  rm(n, fold, iter)
  #output
  output <- list("lamdafull" = lambda_full, "lambdacv" = lambda_cv, "CVstats" = cbind(R2, RMSEP, Rho))
  return(output)
}
#############
# GBLUP.r
####
# Y is a vector of phenotypic records
# X is a matrix of predictors (SNPs coded 0,1,2)
# fold is an integer defining the number of folds from the cross-validation
# iter is an integer defining the number of times the cross-validation is repeated
# cores is an integer defing the number of computer cores to use for the analysis
GBLUP <- function(Y, X, fold = 5, iter = 50, cores = 1){
  require(doParallel)
  require(corpcor)
  require(emma)
  #cores
  cl <- makeCluster(cores)
  #checks
  stopifnot(length(Y) == nrow(X))
  stopifnot(cor(match(rownames(X), names(Y)), 1:length(Y)) == 1)
  #na.omit
  data <- na.omit(data.frame(Y, I(X)))
  n <- nrow(data)
  #Kinship matrix on all data
  ##allele freqs
  p_freq <- colMeans(data$X / 2)
  q_freq <- 1 - p_freq
  K <- tcrossprod(scale(data$X, center = TRUE, scale = FALSE)) / (2 * sum(p_freq * q_freq))
  n <- ncol(K)
  ##k matrix normalization
  K <- (n - 1) / sum((diag(n) - matrix(1, n, n) / n) * K) * K
  ## check if it is positive definite
  if (!is.positive.definite(K)) {K <- make.positive.definite(K)}
  #mixed-model on all data
  ##define an intercept
  data$int <- 1
  mixmod <- emma.REMLE(y = data$Y, X = as.matrix(data$int), K = K)
  ##genomic h2
  h2 <- mixmod$vg / (mixmod$vg + mixmod$ve)
  #repeated cross-validation
  registerDoParallel(cl)
  Y_pred <- foreach(i = 1:iter, .combine = cbind, .packages = c("corpcor", "emma")) %dopar% {
    sampling <- sample(1:n)
    y_hat <- list()
    #cross-validation
    for (k in 1:fold) {
      set <- sort(sampling[floor((k - 1) * n / fold + 1):floor(k * n / fold)])
      #training
      data_cv <- data[-set, ]
      #Kinship matrix in the training set
      p_freq_cv <- colMeans(data_cv$X / 2)
      q_freq_cv <- 1 - p_freq_cv
      k_cv <- tcrossprod(scale(data_cv$X, center = TRUE, scale = FALSE)) / (2 * sum(p_freq_cv * q_freq_cv))
      n_cv <- ncol(k_cv)
      k_cv <- (n_cv - 1) / sum((diag(n_cv) - matrix(1, n_cv, n_cv) / n_cv) * k_cv) * k_cv
      if (!is.positive.definite(k_cv)) {k_cv <- make.positive.definite(k_cv)}
      #mixed-model in the training set
      mixmod <- emma.REMLE(y = data_cv$Y, X = as.matrix(data_cv$int), K = k_cv)
      delta <- mixmod$delta
      #prediction
      if (is.infinite(delta)) {
        y_hat[[k]] <- as.matrix(rep(NA, length(set)))
      } else {
        int <- as.matrix(data_cv$int)
        z <- contrasts(as.factor(as.character(rownames(data))), contrasts = FALSE)[-set, ]
        lhs <- cbind(rbind(crossprod(int), crossprod(z, int)),
                     rbind(crossprod(int, z), crossprod(z) + delta * solve(K)))
        if (!is.positive.definite(lhs)) {lhs <- make.positive.definite(lhs)}
        rhs <- rbind(crossprod(int, as.matrix(data_cv$Y)), crossprod(z, as.matrix(data_cv$Y)))
        results <- solve(lhs) %*% rhs
        beta_hat <- results[1]
        a_hat <- results[2:length(results)]
        y_hat[[k]] <- as.matrix(rep(1, length(set))) %*% beta_hat + a_hat[set]
        rm(int, z, lhs, rhs, results, beta_hat, a_hat)
      }
      rownames(y_hat[[k]]) <- rownames(data)[set]
      rm(set, data_cv, k_cv)
    }
    #gather predictions of each cv repetition
    y_hat <- do.call("rbind", y_hat)
    return(y_hat[match(rownames(data), rownames(y_hat)), 1])
  }
  stopCluster(cl)
  #CV Statistics
  PRESS <- colSums((Y_pred - data$Y)^2)
  R2 <- 1 - (PRESS / (var(data$Y) * (n - 1)))
  RMSEP <- sqrt(PRESS / n)
  Rho <- cor(Y_pred, data$Y)[, 1]
  rm(n, fold, iter)
  #output
  output <- list("h2" = h2, "CVstats" = cbind(R2, RMSEP, Rho))
  return(output)
}
#############
# BL.r
####
# Y is a vector of phenotypic records
# X is a matrix of predictors (SNPs coded 0,1,2)
# fold is an integer defining the number of folds from the cross-validation
# iter is an integer defining the number of times the cross-validation is repeated
# cores is an integer defing the number of computer cores to use for the analysis
# nIterBGLR, burnInBGLR, thinBGLR are inetegers giving the number of iterations, burn-in and thinning of the bayesian chain
# h2 is a prior for the trait heritability
BL <- function(Y, X, fold = 5, iter = 50, cores = 1, nIterBGLR = 2000, burnInBGLR = 1000, thinBGLR = 1, h2 = 0.5){
  require(doParallel)
  require(BGLR)
  #cores
  cl <- makeCluster(cores)
  #checks
  stopifnot(length(Y) == nrow(X))
  stopifnot(cor(match(rownames(X), names(Y)), 1:length(Y)) == 1)
  #na.omit
  data <- na.omit(data.frame(Y, X))
  n <- nrow(data)
  #repeated cross-validation
  registerDoParallel(cl)
  Y_pred <- foreach(i = 1:iter, .combine = cbind, .packages = c("BGLR")) %dopar% {
    folds <- sample(1:fold, size = n, replace = T)
    y_hat <- rep(NA, n)
    #cross-validation
    for (k in 1:max(folds)) {
      tst <- which(folds == k)
      #training set
      yNA <- data$Y
      yNA[tst] <- NA
      #specifying the regression function
      #lambda is defined according to h2 and MSx as defined in table 1 of de los Campos et al. 2013
      MSx <- mean(rowSums(scale(data$X, center = TRUE, scale = FALSE)^2))
      lbd <- sqrt(2 * MSx * (1 - h2) / h2)
      ETA_BL <- list(list(X = data$X, model = 'BL', lambda = lbd, type = 'FIXED'))
      #BL regression
      BLmod <- BGLR(y = yNA, ETA = ETA_BL,
                    nIter = nIterBGLR, burnIn = burnInBGLR, thin = thinBGLR,
                    S0 = NULL, df0 = 5, R2 = 0.5,
                    verbose = FALSE, weights = NULL, rmExistingFiles = TRUE, groups = NULL)
      #prediction within test set
      y_hat[tst] <- BLmod$yHat[tst]
      rm(tst, yNA, BLmod)
    }
    return(y_hat)
  }
  stopCluster(cl)
  #CV Statistics
  PRESS <- colSums((Y_pred - data$Y)^2)
  R2 <- 1 - (PRESS / (var(data$Y) * (n - 1)))
  RMSEP <- sqrt(PRESS / n)
  Rho <- cor(Y_pred, data$Y)[, 1]
  rm(n, fold, iter)
  #output
  output <- list("CVstats" = cbind(R2, RMSEP, Rho))
  return(output)
}
