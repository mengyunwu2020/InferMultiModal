library(glmnet)
library(flare)
library(Rglpk)
library(expm)
li_estimate_f<-function(X,XX,K,n,p,Y,Y2,Y3,Y4,Y5,Y6,Y7){
  ##Li 2022
  hatf_all=list()
  hatU2=NULL
  rr_hat=rep(0,K)
  for(k in 1:K){
    cov<-XX[[k]]%*%t(XX[[k]])/n
    eigvec<-eigen(cov)$vectors #eigenvectors
    eigval<-eigen(cov)$values  #eigenvalues
    K_est<-which.min(diff(log(eigval[1:10]))) #estimate factor numbers using eigenvalue ratios.
    rr_hat[k]=K_est
    hatf2<-sqrt(n)*eigvec[,1:K_est] #estimated factors
    hatf_all[[k]]=hatf2
    hatB2<-t(1/n*t(hatf2)%*%XX[[k]]) #Estimated Factor Loading
    tmp=XX[[k]]-hatf2%*%t(hatB2)
    hatU2<-cbind(hatU2,tmp)   #Estimated Idiosyncratic component
  }
  #y <- y - mean(y)
  X_1=X
  fit = glmnet(X_1,Y, intercept=FALSE,
                lambda=cv.glmnet(X_1,Y,intercept=FALSE)$lambda.1se)
  lambda_beta<-fit$lambda
  beta_hat<-as.vector(fit$beta)
  y_pred <- X_1 %*% beta_hat
  square_loss <- sqrt(mean((Y - y_pred)^2))


  fit = glmnet(X_1,Y2, intercept=FALSE,
               lambda=cv.glmnet(X_1,Y2,intercept=FALSE)$lambda.1se)
  lambda_beta<-fit$lambda
  beta_hat<-as.vector(fit$beta)
  y_pred <- X_1 %*% beta_hat
  square_loss2 <- sqrt(mean((Y2 - y_pred)^2))


  fit = glmnet(X_1,Y3, intercept=FALSE,
               lambda=cv.glmnet(X_1,Y3,intercept=FALSE)$lambda.1se)
  lambda_beta<-fit$lambda
  beta_hat<-as.vector(fit$beta)
  y_pred <- X_1 %*% beta_hat
  square_loss3 <- sqrt(mean((Y3 - y_pred)^2))



  fit = glmnet(X_1,Y4, intercept=FALSE,
               lambda=cv.glmnet(X_1,Y4,intercept=FALSE)$lambda.1se)
  lambda_beta<-fit$lambda
  beta_hat<-as.vector(fit$beta)
  y_pred <- X_1 %*% beta_hat
  square_loss4 <- sqrt(mean((Y4 - y_pred)^2))



  fit = glmnet(X_1,Y5, intercept=FALSE,
               lambda=cv.glmnet(X_1,Y5,intercept=FALSE)$lambda.1se)
  lambda_beta<-fit$lambda
  beta_hat<-as.vector(fit$beta)
  y_pred <- X_1 %*% beta_hat
  square_loss5 <- sqrt(mean((Y5 - y_pred)^2))



  fit = glmnet(X_1,Y6, intercept=FALSE,
               lambda=cv.glmnet(X_1,Y6,intercept=FALSE)$lambda.1se)
  lambda_beta<-fit$lambda
  beta_hat<-as.vector(fit$beta)
  y_pred <- X_1 %*% beta_hat
  square_loss6 <- sqrt(mean((Y6 - y_pred)^2))


  fit = glmnet(X_1,Y7, intercept=FALSE,
               lambda=cv.glmnet(X_1,Y7,intercept=FALSE)$lambda.1se)
  lambda_beta<-fit$lambda
  beta_hat<-as.vector(fit$beta)
  y_pred <- X_1 %*% beta_hat
  square_loss7 <- sqrt(mean((Y7 - y_pred)^2))

  return(list(hatf_all=hatf_all,sigma=c(square_loss,square_loss2,square_loss3,square_loss4,square_loss5,square_loss6,square_loss7)))


}
li_betam<-function(y,X_minus_m,fhat_m,lambda2,sigmahat){
  X_full <- cbind(X_minus_m, fhat_m)
  penalty <- c(rep(1, ncol(X_minus_m)), rep(0, ncol(fhat_m)))
  fit <- glmnet(X_full, y, alpha = 1, penalty.factor = penalty, intercept = FALSE)

  cvfit <- cv.glmnet(X_full, y, alpha = 1, penalty.factor = penalty, intercept = FALSE)
  lambda1_opt <- cvfit$lambda.min


  beta_gamma_hat <- coef(fit, s = lambda1_opt)[-1]
  beta_hat_minus_m <- beta_gamma_hat[1:ncol(X_minus_m)]
  gamma_hat_m <- beta_gamma_hat[(ncol(X_minus_m) + 1):length(beta_gamma_hat)]


  n <- nrow(X_minus_m)
  p <- ncol(X_minus_m)
  q <- ncol(fhat_m)

  W_hat <- matrix(NA, nrow = p, ncol = q)



  W_hat <- danz(X_minus_m,fhat_m,lambda_seq=lambda2)


  sum_xf=matrix(0,p,q)
  tmp=rep(0,q)
  sum_fhat_m=matrix(0,q,q)
  for(i in 1:n){
    tmp=tmp+(y[i]-sum(X_minus_m[i,]*beta_hat_minus_m))*(matrix(fhat_m[i,],ncol=1)-t(W_hat)%*%matrix(X_minus_m[i,],ncol=1))
    sum_fhat_m=sum_fhat_m+matrix(fhat_m[i,],ncol=1)%*%matrix(fhat_m[i,],nrow=1)
    sum_xf=sum_xf+matrix(X_minus_m[i,],ncol=1)%*%matrix(fhat_m[i,],nrow=1)
  }
  sum_fhat_m=1/n*sum_fhat_m
  sum_xf=t(W_hat)%*%sum_xf/n

  Shat=1/(n*sigmahat^2)*tmp
  Ihat=1/(sigmahat^2)*(sum_fhat_m-sum_xf)


  Tn_2=n*t(Shat)%*%solve(Ihat)%*%Shat
  cc=qchisq(0.95, df=q)

  test_li=as.numeric(Tn_2)>cc
  test_li

}


solve_dantzig_lp <- function(X, f_j, lambda2) {


  n <- nrow(X)
  p <- ncol(X)

  obj <- rep(1, 2 * p)
  XtX <- crossprod(X) / n
  Xtf <- crossprod(X, f_j) / n


  A1 <- cbind(XtX, -XtX)
  A2 <- cbind(-XtX, XtX)
  A <- rbind(A1, A2)
  b <- c(Xtf + lambda2, -Xtf + lambda2)

  bounds <- list(lower = list(ind = 1:(2*p), val = rep(0, 2*p)),
                 upper = list(ind = 1:(2*p), val = rep(Inf, 2*p)))


  dir <- rep("<=", 2 * p)

  res <- Rglpk_solve_LP(obj = obj, mat = A, dir = dir, rhs = b,
                        bounds = bounds, max = FALSE)

  u <- res$solution[1:p]
  v <- res$solution[(p+1):(2*p)]
  return(u - v)
}


fit_dantzig_column <- function(X, Y_col, lambda_val) {
  model <- slim(X, Y_col, method = "dantzig", lambda = lambda_val)
  return(as.vector(model$beta))
}

cv_dantzig <- function(X, Y, lambda_seq, nfolds = 3) {
  n <- nrow(X)
  p <- ncol(Y)
  folds <- sample(rep(1:nfolds, length.out = n))
  val_errors <- matrix(0, length(lambda_seq), nfolds)

  for (fold in 1:nfolds) {
    test_idx <- which(folds == fold)
    train_idx <- setdiff(1:n, test_idx)

    X_train <- X[train_idx, ]
    X_test <- X[test_idx, ]
    Y_train <- Y[train_idx, ]
    Y_test <- Y[test_idx, ]

    for (l in seq_along(lambda_seq)) {
      lambda <- lambda_seq[l]
      W_hat <- matrix(0, ncol(X), ncol(Y))

      for (j in 1:ncol(Y)) {
        W_hat[, j] <- fit_dantzig_column(X_train, Y_train[, j], lambda)
      }

      Y_pred <- X_test %*% W_hat
      val_errors[l, fold] <- mean((Y_test - Y_pred)^2)
    }
  }

  mean_errors <- rowMeans(val_errors)
  best_lambda <- lambda_seq[which.min(mean_errors)]
  return(list(lambda = best_lambda, errors = mean_errors))
}





danz<-function(X,Y,lambda_seq= exp(seq(log(1), log(0.01), length.out = 10))){

  lambda=lambda_seq
  if(length(lambda_seq)!=1){
    res <- cv_dantzig(X, Y, lambda_seq)
    lambda <- res$lambda
  }
  W_hat <- matrix(0, ncol(X), ncol(Y))

  for (j in 1:ncol(Y)) {
    W_hat[, j] <- fit_dantzig_column(X, Y[, j], lambda)
  }
  W_hat

}


li_estimate_f1<-function(X,XX,K,n,p,Y7){
  ##Li 2022
  hatf_all=list()
  hatU2=NULL
  rr_hat=rep(0,K)
  for(k in 1:K){
    cov<-XX[[k]]%*%t(XX[[k]])/n
    eigvec<-eigen(cov)$vectors #eigenvectors
    eigval<-eigen(cov)$values  #eigenvalues
    K_est<-which.min(diff(log(eigval[1:10]))) #estimate factor numbers using eigenvalue ratios.
    rr_hat[k]=K_est
    hatf2<-sqrt(n)*eigvec[,1:K_est] #estimated factors
    hatf_all[[k]]=hatf2
    hatB2<-t(1/n*t(hatf2)%*%XX[[k]]) #Estimated Factor Loading
    tmp=XX[[k]]-hatf2%*%t(hatB2)
    hatU2<-cbind(hatU2,tmp)   #Estimated Idiosyncratic component
  }
  #y <- y - mean(y)
  X_1=X



  fit = glmnet(X_1,Y7, intercept=FALSE,
               lambda=cv.glmnet(X_1,Y7,intercept=FALSE)$lambda.1se)
  lambda_beta<-fit$lambda
  beta_hat<-as.vector(fit$beta)
  y_pred <- X_1 %*% beta_hat
  square_loss7 <- sqrt(mean((Y7 - y_pred)^2))

  return(list(hatf_all=hatf_all,sigma=c(square_loss7)))


}
