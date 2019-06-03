### EM wold cross validation code

## packages used in this script
library(pcaMethods)
library(logisticPCA)

## select thresholds to binarize continuous estimation
opt_thres_selection = function(X, Xhat){
  ###
  # select thresholds by minimizing the total and balanced error
  # Input: 
  #       X: binary data X; Xhat: continuous estimation of X
  # Output: 
  #       total_thre:     threshold for total error  
  #       total_error:    optimal total error
  #       balanced_thre:  threshold for balanced error
  #       balanced_error: optimal balanced error
  ###
  
  mn     = length(X)     # length of X 
  p_zero = 1-sum(X)/(mn) # proportion of 0 in X
  search_length = 100;   # searching length
  search_seq = seq(0,1, length.out=search_length) # seaching space
  f_positive = numeric(search_length) # positive error rate
  f_negative = numeric(search_length) # negative error rate
  
  ## compute positive and negative error rates for different thresholds
  for (i in 1:search_length){         
    X_binary = Xhat                  
    X_binary[Xhat > search_seq[i]]  = 1
    X_binary[Xhat <= search_seq[i]] = 0
    X_minus = X-X_binary
    f_positive[i] = sum(X_minus == -1)/(mn)
    f_negative[i] = sum(X_minus == 1)/(mn)
  }
  balanced_error = 0.5*(f_positive/p_zero + f_negative/(1-p_zero)) # balanced error
  total_error    = f_positive + f_negative
  
  opt_thres_balanced = mean(search_seq[which(balanced_error==min(balanced_error))]) # optimal threshold
  opt_error_balanced = balanced_error[which.min(balanced_error)] # optimal balanced error
  
  opt_thres_total = mean(search_seq[which(total_error==min(total_error))]) # optimal threshold
  opt_error_total = total_error[which.min(total_error)] # optimal total error
  
  results = list("total_thres"    = opt_thres_total,
                 "total_error"    = opt_error_total,
                 "balanced_thres" = opt_thres_balanced,
                 "balanced_error" = opt_error_balanced)
  return(results)}

## compute balanced error for cross validation
opt_cv_error = function(X, Xhat, thresholds){
  ###
  # mainly used to compute the total and balanced error for CV sets
  # Input: 
  #       X: CV data sets
  #       Xhat: continuous estimation of CV data sets 
  #       thresholds:
  #           thresholds[1]: threshold for total error computed using training data
  #           thresholds[2]: threshold for balanced error computed using training data
  # Output: 
  #       cv_total_error:    total error for CV data sets
  #       cv_balanced_error: balanced error for CV data sets  
  ###
  
  mn     = length(X)
  p_zero = 1-sum(X)/mn
  
  ## compute total cv error 
  Xhat_total = Xhat
  Xhat_total[Xhat>thresholds[1]]  = 1 # using threshold for total error
  Xhat_total[Xhat<=thresholds[1]] = 0
  X_minus_Xhat_total = X-Xhat_total
  cv_total_error     = sum(X_minus_Xhat_total != 0)/mn
  
  ## compute balanced cv error
  # binarize Xhat to binary data according to provided threshold
  Xhat_balanced = Xhat
  Xhat_balanced[Xhat>thresholds[2]]  = 1
  Xhat_balanced[Xhat<=thresholds[2]] = 0
  
  X_minus_Xhat_balanced = X - Xhat_balanced
  f_poistive = sum(X_minus_Xhat_balanced==-1)/mn
  f_negative = sum(X_minus_Xhat_balanced==1)/mn
  
  cv_balanced_error = 0.5*(f_poistive/p_zero + f_negative/(1-p_zero))
  
  result = list("cv_total_error"    = cv_total_error,
                "cv_balanced_error" = cv_balanced_error)
  return(result)}

## select sets for K-fold EM wold cross validation 
select_sets = function(n, p, K){
  ###
  # used to select K folds for EM wold CV
  # Input:
  #       n: number of samples; p: number of variables; K: K-fold CV
  # Output: 
  #       sets_list: a list contains K sets index         
  ###
  
  sets_list=list()
  for (i in 1:K){
    sets = c()
    temp = i
    while(temp <= n*p){
      sets = c(sets, temp)
      temp = temp + K
    }
    temp_name      = paste("sets", i, sep = "")
    sets_list[[i]] = sets
  }
  return(sets_list)
}

## EM-wold cross validation for PCA model
CV_pca_wold = function(X, K, r){
  ###
  # doing EM-wold CV using pca method 
  # Input:
  #       X: binary data
  #       K: K-fold CV;   r: number of PCs
  # Output: 
  #       cv_error:
  #                cv_error[1]: mean total error for CV sets
  #                cv_error[2]: mean balanced error for CV sets
  #       opt_thres:
  #                opt_thres[1]: threshold for total error 
  #                opt_thres[2]: threshold for balanced error
  #       train_error: 
  #                train_error[1]: total error for training set
  #                train_error[2]: balanced error for training set
  ###
  
  X = as.matrix(X)
  misclass_error_total    = numeric(K)
  misclass_error_balanced = numeric(K)
  
  ## select K-fold CV data sets
  m = dim(X)[1]; n = dim(X)[2]
  CV_sets = select_sets(m, n, K)
  
  ## infer optmal thresholds using traing data set
  if(r>0){
    full_pca      = pca(X, method="svd", 
	                    scale="none", center=TRUE, nPcs=r) # pca model  
    X_fitted_full = fitted(full_pca, nPcs=r)
    misclass_error_thres = opt_thres_selection(X, X_fitted_full)
  } else{
    full_pca = as.matrix(rep(1, dim(X)[1])) %*% t(as.matrix(colMeans(X, na.rm=TRUE)))
    X_fitted_full = full_pca
    misclass_error_thres = opt_thres_selection(X, X_fitted_full)
  }
  
  ## index out thresholds
  opt_thres   = c(misclass_error_thres$total_thres, misclass_error_thres$balanced_thres)
  
  ## index out training errors
  train_error = c(misclass_error_thres$total_error, misclass_error_thres$balanced_error)
  
  ## doing K-fold CV
  for (i in 1:K){
    missing_sets = CV_sets[[i]] 
    X_minus = X
    X_minus[missing_sets] = NA
	
    if (sum(is.na(X_minus)) <= 1){
      stop("missing values are not added")}
	  
    if(r>0){
      temp_pca    = pca(X_minus, method="svdImpute",
                	    scale="none", center=TRUE, nPcs=r) # pca model  
      X_estimated = completeObs(temp_pca) #input missing value
    } else {
      temp_pca = as.matrix(rep(1, dim(X)[1])) %*% t(as.matrix(colMeans(X_minus, na.rm=TRUE)))
      X_estimated = temp_pca
    }

    ## misclassfication error
    missing_estimated = X_estimated[missing_sets]
    missing_real      = X[missing_sets]
    cv_error = opt_cv_error(missing_real, missing_estimated, opt_thres)
    misclass_error_total[i]    = cv_error$cv_total_error
    misclass_error_balanced[i] = cv_error$cv_balanced_error
  }
  misclass_error = c(mean(misclass_error_total), mean(misclass_error_balanced))
  names(misclass_error) = c("total", "balanced")
  names(train_error)    = c("total", "balanced")
  names(opt_thres)      = c("total", "balanced")
  
  return(list("cv_error"    = misclass_error, 
              "train_error" = train_error, 
              "opt_thres"   = opt_thres))
}

## EM-wold cross validation for Gifi model 
CV_gifi_wold = function(X, K, r){
  ###
  # doing EM-wold CV using gifi method 
  # Input:
  #       X: binary data;
  #       K: K-fold CV;   r: number of PCs
  # Output: 
  #       cv_error:
  #                cv_error[1]: mean total error for CV sets
  #                cv_error[2]: mean balanced error for CV sets
  #       opt_thres:
  #                opt_thres[1]: threshold for total error 
  #                opt_thres[2]: threshold for balanced error
  #       train_error: 
  #                train_error[1]: total error for training set
  #                train_error[2]: balanced error for training set
  ###
  
  X = as.matrix(X)
  misclass_error_total    = numeric(K)
  misclass_error_balanced = numeric(K)
  
  ## select K-fold CV data sets
  m = dim(X)[1]; n = dim(X)[2]
  CV_sets = select_sets(m, n, K)
  
  ## infer optimal thresholds
  if(r>0){
    full_gifi     = pca(X, method = "svd",
	                    scale = "uv", center = TRUE, nPcs = r) # gifi model  
    X_fitted_full = fitted(full_gifi, nPcs=r)
    misclass_error_thres = opt_thres_selection(X, X_fitted_full)
  } else{
    full_gifi = as.matrix(rep(1, dim(X)[1])) %*% t(as.matrix(colMeans(X, na.rm=TRUE)))
    X_fitted_full = full_gifi
    misclass_error_thres = opt_thres_selection(X, X_fitted_full)
  }
  
  ## index out thresholds
  opt_thres   = c(misclass_error_thres$total_thres, misclass_error_thres$balanced_thres)
  
  ## index out training errors
  train_error = c(misclass_error_thres$total_error, misclass_error_thres$balanced_error)
  
  ## doing K-fold CV
  for (i in 1:K){
    missing_sets = CV_sets[[i]]
    X_minus      = X
    X_minus[missing_sets] = NA
	
    if (sum(is.na(X_minus)) != length(missing_sets)){
      stope("missing values are not added")}
	  
    # using scaled PCA for gifi's method
    if(r>0){
      temp_gifi   = pca(X_minus, method="svdImpute", 
	                    scale="uv", center=TRUE, nPcs=r)  
      X_estimated = completeObs(temp_gifi) # input missing value
    } else {
      temp_gifi   = as.matrix(rep(1, dim(X)[1])) %*% t(as.matrix(colMeans(X_minus, na.rm=TRUE)))
      X_estimated = temp_gifi
    }
    
    ## misclassfication error
    missing_estimated = X_estimated[missing_sets]
    missing_real      = X[missing_sets]
    cv_error = opt_cv_error(missing_real, missing_estimated, opt_thres)
    misclass_error_total[i]    = cv_error$cv_total_error
    misclass_error_balanced[i] = cv_error$cv_balanced_error
  }
  misclass_error = c(mean(misclass_error_total), mean(misclass_error_balanced))
  
  names(misclass_error) = c("total", "balanced")
  names(train_error)    = c("total", "balanced")
  names(opt_thres)      = c("total", "balanced")
  
  return(list("cv_error"    = misclass_error, 
              "train_error" = train_error, 
              "opt_thres"   = opt_thres))
}

## EM-wold cross validation for logistic PCA model 
CV_logpca_wold = function(X, K, r){
  ###
  # doing EM-wold CV using logistic PCA method 
  # Input:
  #       X: binary data;
  #       K: K-fold CV;   r: number of PCs
  # Output: 
  #       cv_error:
  #                cv_error[1]: mean total error for CV sets
  #                cv_error[2]: mean balanced error for CV sets
  #       opt_thres:
  #                opt_thres[1]: threshold for total error 
  #                opt_thres[2]: threshold for balanced error
  #       train_error: 
  #                train_error[1]: total error for training set
  #                train_error[2]: balanced error for training set
  ###
  
  X = as.matrix(X)
  misclass_error_total    = numeric(K)
  misclass_error_balanced = numeric(K)
  m = 2.95   # approximation of natural parameter from staturated model
  
  ## select K-fold CV data sets
  m = dim(X)[1]; n = dim(X)[2]
  CV_sets = select_sets(m, n, K)
  
  ## infer optimal thresholds
  if(r>0){
    full_logpca   = logisticPCA(X, k=r, m=m,
                                main_effects=TRUE, random_start=FALSE)
    X_fitted_full = fitted(full_logpca, type="response") # input missing value
    misclass_error_thres = opt_thres_selection(X, X_fitted_full)
  }else{
    full_logpca   = as.matrix(rep(1, dim(X)[1])) %*% t(as.matrix(colMeans(X, na.rm=TRUE)))
    X_fitted_full = full_logpca
    misclass_error_thres = opt_thres_selection(X, X_fitted_full)
  }
  
  ## index out thresholds
  opt_thres   = c(misclass_error_thres$total_thres, misclass_error_thres$balanced_thres)
  
  ## index out training errors
  train_error = c(misclass_error_thres$total_error, misclass_error_thres$balanced_error)
  
  ## doing K-fold CV
  for (i in 1:K){
    missing_sets = CV_sets[[i]] 
    X_minus      = X
    X_minus[missing_sets] = NA
	
    if (sum(is.na(X_minus)) != length(missing_sets)){
      stop("missing values are not added")}
    
    if(r>0){
      temp_logpca = logisticPCA(X_minus, k=r, m=m,
                                main_effects=TRUE, random_start=FALSE)
      X_estimated = fitted(temp_logpca, type="response") # input missing value
    } else {
      temp_logpca = as.matrix(rep(1, dim(X)[1])) %*% t(as.matrix(colMeans(X_minus, na.rm=TRUE)))
      X_estimated = temp_logpca
    }
    
    ## misclassfication error
    missing_estimated = X_estimated[missing_sets]
    missing_real      = X[missing_sets]
    cv_error = opt_cv_error(missing_real, missing_estimated, opt_thres)
    misclass_error_total[i]    = cv_error$cv_total_error
    misclass_error_balanced[i] = cv_error$cv_balanced_error
  }
  misclass_error = c(mean(misclass_error_total), mean(misclass_error_balanced))
  
  names(misclass_error) = c("total", "balanced")
  names(train_error)    = c("total", "balanced")
  names(opt_thres)      = c("total", "balanced")
  
  return(list("cv_error"    = misclass_error, 
              "train_error" = train_error, 
              "opt_thres"   = opt_thres))
}

## EM-wold cross validation for logistic SVD model 
CV_logsvd_wold = function(X, K, r){
  ###
  # doing EM-wold CV using logistic SVD method 
  # Input:
  #       X: binary data; 
  #       K: K-fold CV;   r: number of PCs
  # Output: 
  #       cv_error:
  #                cv_error[1]: mean total error for CV sets
  #                cv_error[2]: mean balanced error for CV sets
  #       opt_thres:
  #                opt_thres[1]: threshold for total error 
  #                opt_thres[2]: threshold for balanced error
  #       train_error: 
  #                train_error[1]: total error for training set
  #                train_error[2]: balanced error for training set
  ###
  
  X = as.matrix(X)
  misclass_error_total    = numeric(K)
  misclass_error_balanced = numeric(K)
  
  ## select K-fold CV data sets
  m = dim(X)[1]; n = dim(X)[2]
  CV_sets = select_sets(m, n, K)
  
  ## infer optimal threshold
  if(r>0){
    full_logsvd   = logisticSVD(X, k=r,
                                main_effects=TRUE, random_start=FALSE)
    X_fitted_full = fitted(full_logsvd, type="response") # input missing value
    misclass_error_thres = opt_thres_selection(X, X_fitted_full)
  }else{
    full_logsvd   = as.matrix(rep(1, dim(X)[1])) %*% t(as.matrix(colMeans(X, na.rm=TRUE)))
    X_fitted_full = full_logsvd
    misclass_error_thres = opt_thres_selection(X, X_fitted_full)
  }
  
  ## index out thresholds
  opt_thres   = c(misclass_error_thres$total_thres, misclass_error_thres$balanced_thres)
  
  ## index out training errors
  train_error = c(misclass_error_thres$total_error, misclass_error_thres$balanced_error)
  
  ## doing K-fold CV
  for (i in 1:K){
    missing_sets = CV_sets[[i]]
    X_minus = X
    X_minus[missing_sets] = NA
	
    if (sum(is.na(X_minus)) != length(missing_sets)){
      stop("missing values are not added")}
    
    if(r>0){
        temp_logsvd = logisticSVD(X_minus, k=r, 
                                  main_effects=TRUE, random_start=FALSE)
        X_estimated = fitted(temp_logsvd, type = "response") #input missing value
    } else {
      temp_logsvd = as.matrix(rep(1, dim(X)[1])) %*% t(as.matrix(colMeans(X_minus, na.rm=TRUE)))
      X_estimated = temp_logsvd
    }
    
    ## misclassfication error
    missing_estimated = X_estimated[missing_sets]
    missing_real      = X[missing_sets]
    cv_error = opt_cv_error(missing_real, missing_estimated, opt_thres)
    misclass_error_total[i]    = cv_error$cv_total_error
    misclass_error_balanced[i] = cv_error$cv_balanced_error
  }
  misclass_error = c(mean(misclass_error_total), mean(misclass_error_balanced))
  
  names(misclass_error) = c("total", "balanced")
  names(train_error)    = c("total", "balanced")
  names(opt_thres)      = c("total", "balanced")
  
  return(list("cv_error"    = misclass_error, 
              "train_error" = train_error, 
              "opt_thres"   = opt_thres))
}










