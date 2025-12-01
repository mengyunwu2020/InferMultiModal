#' Estimation function for the proposed model.
#'
#' @param X A list of the multi-modality data, each element is a data matrix of each group with dimension n*p_m.
#' @param y A response vector.
#' @param K The number of modalities.
#' @param rmax The maximum factor numbers of all modalities.
#' @param method The method used in the APM algorithm, default is ic, it can also be gap.
#' @param weight The weight of each projection matrix in the APM algorithm, default is TRUE. if weight = FALSE, then simply calculate the mean of all projection matrices.
#' @param type The method used in estimating the factor numbers in each modality initially, default is IC3.
#' @param lambda.seq A user supplied non-negative tuning parameters.
#'
#' @return An object of class "estimation" containing the fitted model.
#' @import ncvreg
#' @import SIS
#' @import Matrix
#' @import GrFA
#' @examples
#' set.seed(1)
#' K=4
#' p=1200
#' library(Matrix)
#' r_f=c(2,rep(2,K))
#' r_all=sum(r_f)
#' pm=rep(p/K,K)
#' gamma1=rep(0.5,r_all)
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
#' @export
#'
#'
estimation<-function(X,y,K,rmax=8,method='ic',weight=TRUE, type = "IC3",lambda.seq=NULL){
  factor_res<-APM(X, rmax = rmax, localfactor = TRUE, method = method,weight=weight, type = type)


  loading_c=factor_res$loading_G
  loading_s=factor_res$loading_F



  fc<-factor_res$Ghat
  fs<-factor_res$Fhat
  v<-factor_res$residual
  if(!anyNA(fc)){
    FV=fc
  }else{
    FV=NULL
  }
  V=NULL
  rr_hat=rep(0,length(fs)+1)
  if(!anyNA(fc)){
    rr_hat[1]=dim(fc)[2]
  }else{
  rr_hat[1]=0
  }
  for(i in 1:length(fs)){
    if(!anyNA(fs[[i]])){
      FV=cbind(FV,fs[[i]])
      rr_hat[i+1]=dim(fs[[i]])[2]
    }
    V=cbind(V,v[[i]])
  }
  r=dim(FV)[2]
  p=dim(V)[2]
  Fhat=FV
  FV=cbind(FV,V)

  penalty_factors=c(rep(0,r),rep(1,p))

  if(is.null(lambda.seq)){
    cvfit <- cv.ncvreg(Re(FV), y,penalty.factor = penalty_factors, penalty="MCP")

  }else{
    cvfit <- cv.ncvreg(Re(FV), y,penalty.factor = penalty_factors, penalty="MCP",lambda = lambda.seq)
  }

  fit <- cvfit$fit
  cv_res<-fit$beta[,cvfit$min]

  res<-list(Fhat=Fhat,V=V,FV=FV,gamma=cv_res[2:(r+1)],beta=cv_res[(r+2):(r+p+1)],rr_hat=rr_hat,loading_c=loading_c,loading_s=loading_s)
  class(res) <- "estimation"
  return(res)
}

