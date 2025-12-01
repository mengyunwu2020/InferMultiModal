#' @title Prediction method for estimation objects
#' @description Defines the predict() S3 method for objects of class "estimation".
#'
#' @param object An object of class "estimation" (returned by estimation()).
#' @param newdata A data frame containing the new predictor values.
#'
#' @return A numeric vector of predicted values.
#' @method predict estimation
#' @import Matrix
#' @export
#'
#' @examples
#'
#' set.seed(1)
#' K=4
#' p=1200
#' library(Matrix)
#' r_f=c(2,rep(2,K))
#' r_all=sum(r_f)
#' pm=rep(p/K,K)
#' gamma1=rep(0.5,r_all)#'c(0,0.3,0.5,0.5,0.5,0,0,0,0,0.5)
#' beta=c(rep(0.5,2),rep(0,pm[1]-2),rep(0.5,2),rep(0,pm[2]-2),rep(0.5,2),rep(0,pm[3]-2),rep(0.5,2),rep(0,pm[4]-2))
#' B_c<-matrix(rnorm(p*r_f[1]),nrow=p)
#' B=B_c
#' B_s<-list()
#' for(k in 1:K){
#'   B_s[[k]]<-matrix(rnorm(pm[k]*r_f[k+1]),nrow=pm[k])
#' }
#'
#' block_diag <- as.matrix(bdiag(B_s[[1]], B_s[[2]]))
#' block_diag <- as.matrix(bdiag(block_diag, B_s[[3]]))
#' block_diag <- as.matrix(bdiag(block_diag, B_s[[4]]))
#' B=cbind(B_c,as.matrix(block_diag))
#' n=ceiling(log(p)/(0.5^2/4^2))
#' F<-matrix(rnorm(n*r_all),nrow=n)
#' V<-matrix(rnorm(n*p),nrow=n)
#' X=F%*%t(B)+V
#' Y=F%*%gamma1+V%*%beta+rnorm(n,mean=0,sd=1)
#'
#' ss=1
#' XX=list()
#' for(k in 1:K){
#'   XX[[k]]=X[,ss:(ss+pm[k]-1)]
#'   ss=ss+pm[k]
#' }
#' fit<-estimation(XX,Y,K,rmax=8,method='ic',weight=TRUE)
#' n_test <- 100
#' F_test <- matrix(rnorm(n_test * r_all), nrow = n_test)
#' V_test <- matrix(rnorm(n_test * p), nrow = n_test)
#' X_test <- F_test %*% t(B) + V_test
#' Y_test <- F_test%*%gamma1+V_test%*%beta+rnorm(n_test,mean=0,sd=1)
#' y_pred <- predict(fit,X_test)

predict.estimation <- function(object, newdata, ...) {

  beta_hat1<-as.vector(object$beta)
  gamma_hat1<-as.vector(object$gamma)
  loading_c=object$loading_c
  loading_s=object$loading_s
  loading_c_mat=loading_c[[1]]
  loading_s_mat=loading_s[[1]]
  K=length(loading_s)
  for(k in 1:(K-1)){
    loading_c_mat=rbind(loading_c_mat,loading_c[[k+1]])
    loading_s_mat=as.matrix(bdiag(loading_s_mat,loading_s[[k+1]]))
  }


  if (missing(newdata)) {
    stop("Please provide 'newdata' for prediction.")
  }


  loading_mat=cbind(loading_c_mat,loading_s_mat)
  BtB_inv <- solve(t(loading_mat) %*% loading_mat)       # (B^T B)^{-1}
  projection <- BtB_inv %*% t(loading_mat)     # (B^T B)^{-1} B^T
  F_hat <- newdata %*% t(projection)       #   f_t^T
  V_hat <- newdata - F_hat%*%t(loading_mat)
  y_predict=V_hat%*%beta_hat1+F_hat%*%gamma_hat1

  return(y_predict)
}

#---------------Main function:selectfeature: use first part of data to conduct sure screening and use the second part of data to construct test statistics---------
#input Y1: response
#X: covariate
#n: number of observations
#ratio1: the ratio in the first stage (sure screening)
#ratio2: the ratio in the second stage (anova test)
selectfeature<-function(Y1,X,n,K,pm,ratio1=0.8,ratio2=0.9){
  ###In the following step, we use PCA of sample covariance matrix to estimate latent factors via the first part of data################
  ss=1
  XX=list()
  for(k in 1:K){
    XX[[k]]=X[1:floor(n^(ratio1)),ss:(ss+pm[k]-1)]
    ss=ss+pm[k]
  }
  p=sum(pm)

  fit1<-estimation(XX,Y1[1:floor(n^(ratio1))],K,rmax=8,method='ic',weight=TRUE)
  r_all<-dim(fit1$FV)[2]-p;
  hatV = fit1$FV[,(r_all+1):(r_all+p)];
  hatF = fit1$FV[,1:(r_all)];
  Y_new = Y1[1:floor(n^(ratio1))] -  hatF%*%fit1$gamma
  sis<-SIS(hatV,Y_new) #use sure screening to select variables

  select1<-sis$ix #ISIS


  rest_n=n-floor(n^(ratio1))
  ss=1
  XX=list()
  for(k in 1:K){
    XX[[k]]=X[(floor(n^(ratio1))+1):n,ss:(ss+pm[k]-1)]
    ss=ss+pm[k]
  }

  factor_res<-APM(XX, rmax = 8, localfactor = TRUE,method='ic',weight=TRUE)

  hatfc<-factor_res$Ghat#estimated common factors
  hatfs<-factor_res$Fhat#estimated specific factors
  hatV<-factor_res$residual#Estimated Idiosyncratic component


  rr_hat=rep(0,length(hatfs)+1)
  rr_hat[1]=dim(hatfc)[2]
  V=fs=NULL
  for(k in 1:length(hatfs)){
    if(!anyNA(hatfs[[k]])){
      fs=cbind(fs,hatfs[[k]])
      rr_hat[k+1]=dim(hatfs[[k]])[2]
    }
    V=cbind(V,hatV[[k]])
  }

  if(!anyNA(hatfc)){
    FV=cbind(hatfc,fs)
  }else{
    FV=fs
  }
  FV=cbind(FV,V[,select1])


  # -------------------------------------------------------------------------
  penalty_factors=c(rep(0,sum(rr_hat)),rep(1,length(select1)))
  cvfit <- cv.ncvreg(FV,Y1[(floor(n^(ratio1))+1):n],penalty.factor = penalty_factors, penalty="MCP")#,intercept = FALSE
  fit <- cvfit$fit
  cv_res<-fit$beta[,cvfit$min]
  beta_hat=cv_res[(sum(rr_hat)+2):(sum(rr_hat)+length(select1)+1)]
  gamma_hat=cv_res[2:(sum(rr_hat)+1)]
  print(beta_hat)
  print(gamma_hat)

  theta1=solve(cov(FV[,1:sum(rr_hat)]))
  #------------------Conduct bootstrap to compute the critical values----------
  quant<-c() #threshold
  for (m in 1:5){
    boostrap<-c()
    for(j in 1:500){
      boostrap<-c(boostrap,max(abs(theta1%*%t(FV[,1:sum(rr_hat)])%*%rnorm(dim(FV)[1],mean=0,sd=1))))
    }
    quant<-c(quant,quantile(boostrap,0.95))
  }
  c_alpha<-1/sqrt(dim(FV)[1])*mean(quant)#0.95 threshold

  n_anova=length(c((floor(n^(ratio1))+1):(floor(n^(ratio1)+rest_n^ratio2))))
  P_all<-FV[1:n_anova,]%*%solve(t(FV[1:n_anova,])%*%(FV[1:n_anova,]))%*%t(FV[1:n_anova,])
  tmp=FV[1:n_anova,-c(1:rr_hat[1])]
  P_noc<-tmp%*%solve(t(tmp)%*%(tmp))%*%t(tmp)
  Q_stat=rep(0,K+1)
  P_diff=P_all-P_noc
  Q_stat[1]<-t(Y1[c((floor(n^(ratio1))+1):(floor(n^(ratio1)+rest_n^ratio2)))])%*%P_diff%*%Y1[c((floor(n^(ratio1))+1):(floor(n^(ratio1)+rest_n^ratio2)))] #test statistics
  #Output: test statistics
  adjust_v<-log(p)/(n^(1-ratio2))
  ss=rr_hat[1]+1
  for(k in 1:K){
    ind1=c(ss:(ss+rr_hat[k+1]-1))
    tmp=FV[1:n_anova,-ind1]
    P_no_sk<-tmp%*%solve(t(tmp)%*%(tmp))%*%t(tmp)
    P_diff=P_all-P_no_sk
    Q_stat[k+1]<-t(Y1[c((floor(n^(ratio1))+1):(floor(n^(ratio1)+rest_n^ratio2)))])%*%P_diff%*%Y1[c((floor(n^(ratio1))+1):(floor(n^(ratio1)+rest_n^ratio2)))] #test statistics
    ss=ss+rr_hat[k+1]
    tmp=log(p)/(n^(1-ratio2))
    adjust_v<-c(adjust_v,tmp)
  }
  Q_res=rbind(Q_stat,adjust_v)

  ss=1
  c_res=NULL
  for(k in 1:(K+1)){
    tmp=sqrt(dim(FV)[1])*max(abs(gamma_hat[ss:(ss+rr_hat[k]-1)]))/c_alpha
    ss=ss+rr_hat[k]
    c_res=c(c_res,tmp)
  }
  Q_res=rbind(Q_stat,c_res)


  return(Q_res)
}

#-------------Helper function: rcv: compute the variance of Y given X-----------
#input Y: response
#X: covariate of all modalities
#n: number of observations
#pm: a vector of the dimensions of all modalities
rcv<-function(Y,X,n,pm,K){
  Idx = 1:(n/2)
  X0 = X[Idx,]; X2 = X[-Idx,]
  Y0 = Y[Idx]; Y2 = Y[-Idx]

  ss=1
  XX=list()
  for(k in 1:K){
    XX[[k]]=X0[,ss:(ss+pm[k]-1)]
    ss=ss+pm[k]
  }
  XX0=XX

  ss=1
  XX=list()
  for(k in 1:K){
    XX[[k]]=X2[,ss:(ss+pm[k]-1)]
    ss=ss+pm[k]
  }
  XX2=XX


  fit1<-estimation(XX0,Y0,K,rmax=8,method='ic',weight=TRUE)
  fit2<-estimation(XX2,Y2,K,rmax=8,method='ic',weight=TRUE)

  r_all<-dim(fit1$FV)[2]-p;  r_all2<-dim(fit2$FV)[2]-p
  V.hat0 = fit1$FV[,(r_all+1):(r_all+p)]; V.hat2 = fit2$FV[,(r_all2+1):(r_all2+p)]
  F.hat0 = fit1$FV[,1:(r_all)]; F.hat2 = fit2$FV[,1:(r_all2)]
  Y0.new = Y0 -  F.hat0%*%fit1$gamma; Y2.new = Y2 -  F.hat2%*%fit2$gamma

  SIS0 = SIS(V.hat0,Y0.new); SIS2 = SIS(V.hat2,Y2.new)
  S.hat0.ISIS = SIS0$ix;    S.hat0.SIS = SIS0$ix0
  S.hat2.ISIS = SIS2$ix;    S.hat2.SIS = SIS2$ix0

  X0=cbind(F.hat0,V.hat0); X2=cbind(F.hat2,V.hat2)
  X0.ISIS = X0[,c(1:r_all,S.hat2.ISIS+r_all)]; X0.SIS = X0[,c(1:r_all,r_all+S.hat2.SIS)]
  X2.ISIS = X2[,c(1:r_all2,r_all2+S.hat0.ISIS)]; X2.SIS = X2[,c(1:r_all2,r_all2+S.hat0.SIS)]

  lm0.ISIS = lm(Y0~X0.ISIS-1); lm0.SIS = lm(Y0~X0.SIS-1)
  lm2.ISIS = lm(Y2~X2.ISIS-1); lm2.SIS = lm(Y2~X2.SIS-1)

  Q0.ISIS.H1 = sum((resid(lm0.ISIS))^2)
  sigma.hat0.ISIS = Q0.ISIS.H1/(n/2 - r_all - length(S.hat2.ISIS))
  Q0.SIS.H1 = sum((resid(lm0.SIS))^2)
  sigma.hat0.SIS = Q0.SIS.H1/(n/2 - r_all - length(S.hat2.SIS))
  Q2.ISIS.H1 = sum((resid(lm2.ISIS))^2)
  sigma.hat2.ISIS = Q2.ISIS.H1/(n/2 - r_all2 - length(S.hat0.ISIS))
  Q2.SIS.H1 = sum((resid(lm2.SIS))^2)
  sigma.hat2.SIS = Q2.SIS.H1/(n/2 - r_all2 - length(S.hat0.SIS))
  sigma.hat.ISIS = mean(c(sigma.hat0.ISIS, sigma.hat2.ISIS))
  sigma.hat.SIS = mean(c(sigma.hat0.SIS, sigma.hat2.SIS))
  c(sigma.hat.ISIS,sigma.hat.SIS)} #output estimated variance using (ISIS method, SIS method), respectively.
