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
#' res_fac<-estimation_fac(XX,K,method='gap')
#' fit<-estimation(XX,Y,K,res_fac=res_fac)
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
