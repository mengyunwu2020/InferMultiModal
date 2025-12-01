#' Hypothesis Test of the Effects of Common and Modality-specific Factors.
#'
#' @param X A list of the multi-modality data, each element is a data matrix of each group with dimension n*p_m.
#' @param y A response vector.
#' @param K The number of modalities.
#' @param ratio1 The exponent used in the sure screening step.
#' @param ratio2 The exponent used in the ANOVA test step.
#' @param rmax The maximum number of factors across all modalities.
#' @param method The method used in the APM algorithm. The default is "ic"; it can also be "gap".
#' @param weight Logical. Whether to use weighted projection matrices in the APM algorithm.
#'   The default is TRUE. If FALSE, the algorithm simply calculates the mean of all projection matrices.
#' @param type The method used to estimate the number of factors in each modality initially.
#'   The default is "IC3".
#' @param lambda.seq A user-supplied sequence of non-negative tuning parameters.
#' @param r_f A vector specifying the number of factors in each modality. The default is NULL.
#' @param alpha The significance level used in hypothesis testing. The default is 0.05.
#'
#' @return A numeric vector containing the p-values testing the significance of
#'   the effects of common and modality-specific factors.
#' @import ncvreg
#' @import SIS
#' @import Matrix
#' @import GrFA
#' @import ACAT
#' @import metap
#' @examples
#' p <- 1200
#' n=300
#' set.seed(1)
#' K=4
#' r_f=c(2,rep(2,K))
#' r_all=sum(r_f)
#' pm=rep(p/K,K)
#' beta=c(rep(1,2),rep(0,pm[1]-2),rep(1,2),rep(0,pm[2]-2),rep(1,2),rep(0,pm[3]-2),rep(1,2),rep(0,pm[4]-2))
#' library(Matrix)
#' B_c<-matrix(rnorm(p*r_f[1],sd=1),nrow=p)
#' B=B_c
#' B_s<-list()
#' for(k in 1:K){
#'   B_s[[k]]<-matrix(rnorm(pm[k]*r_f[k+1],sd=1),nrow=pm[k])
#' }
#' block_diag <- as.matrix(bdiag(B_s[[1]], B_s[[2]]))
#' block_diag <- as.matrix(bdiag(block_diag, B_s[[3]]))
#' block_diag <- as.matrix(bdiag(block_diag, B_s[[4]]))
#' B=cbind(B_c,as.matrix(block_diag))
#' gamma1=c(rep(2,r_all))
#' set.seed(1)
#' F<-matrix(rnorm(n*r_all),nrow=n)
#' V<-matrix(rnorm(n*p),nrow=n)
#' X=F%*%t(B)+V #'design matrix
#' y=V%*%beta+rnorm(n,mean=0,sd=1)
#' y1=F%*%gamma1+y
#'
#' ss=1
#' XX=list()
#' for(k in 1:K){
#'   XX[[k]]=X[,ss:(ss+pm[k]-1)]
#'   ss=ss+pm[k]
#' }
#' res=gammatest(XX,y,K,ratio1=0.9,ratio2=0.6,rmax=8,method='ic',weight=TRUE,lambda.seq=NULL,r_f=NULL,alpha=0.05,type="IC3")
#' @export
gammatest<-function(X,y,K,ratio1,ratio2,rmax=8,method='ic',weight=TRUE,type="IC3",lambda.seq=NULL,r_f=NULL,alpha=0.05){


  X_1=X[[1]]
  pm=rep(0,K)
  pm[1]=ncol(X_1)
  for(k in 2:K){
    X_1=cbind(X_1,X[[k]])
    pm[k]=ncol(X[[k]])
  }
  n=nrow(X_1)
  sig_est<- try(
    suppressMessages(suppressWarnings(
      {
        tmp <- capture.output(out <- rcv(y,X_1,n,pm,K,rmax=rmax,method=method,weight=weight,lambda.seq = lambda.seq, type = type), file = NULL)
        out
      }
    )),
    silent = TRUE
  ) #estimate variance through function rcv

  var_est1<-sig_est[1] #estimated variance via isis.

  output11<- selectfeature2(y,X_1,n,K,pm,ratio1=ratio1,ratio2=ratio2,rf=r_f,rmax=rmax,method=method,weight=weight,lambda.seq = lambda.seq,alpha=alpha, type = type)
  output1 <- output11$Q_stat
  p=sum(pm)

  pvals=rep(0,(K+1))

  r_f=output11$rf
  for(kk in 1:(K+1)){
    tmp=pchisq(q = output1[,kk]/var_est1, df = r_f[kk], lower.tail = FALSE)


    restmp=c(ACAT(tmp),maximump(tmp)$p)

    pvals[kk]=mean(restmp)


  }

names(pvals)=c('common',paste0('specific',1:K))
      return(pvals)

}





selectfeature2<-function(Y1,X,n,K,pm,ratio1=0.8,ratio2=0.9,rf,rmax=8,method='ic',weight=TRUE,lambda.seq =NULL,alpha=0.05, type = "IC3"){

  time_test= floor((n-floor(n^(ratio1)))/n^ratio2)
  ss=1
  XX=list()
  for(k in 1:K){
    XX[[k]]=X[1:floor(n^(ratio1)),ss:(ss+pm[k]-1)]
    ss=ss+pm[k]
  }
  p=sum(pm)

  fit1<-estimation(XX,Y1[1:floor(n^(ratio1))],K,rmax=rmax,method=method,weight=weight,lambda.seq = lambda.seq, type = type)
  r_all<-dim(fit1$FV)[2]-p;
  hatV = fit1$FV[,(r_all+1):(r_all+p)];
  hatF = fit1$FV[,1:(r_all)];
  Y_new = Y1[1:floor(n^(ratio1))] -  hatF%*%(fit1$gamma)
  sis<-try(
    suppressMessages(suppressWarnings(
      {
        tmp <- capture.output(out <- SIS(hatV, Y_new), file = NULL)
        out
      }
    )),
    silent = TRUE
  ) #use sure screening to select variables

  select1<-sis$ix #ISIS

  if(length(select1)>n^ratio2-r_all){
    tmp=order(abs(sis$coef.est)[-1],decreasing=TRUE)
    select1=select1[tmp[1:(floor(n^ratio2-r_all)-1)]]
  }


  rest_n=n-floor(n^(ratio1))
  ss=1
  XX=list()
  for(k in 1:K){
    XX[[k]]=X[(floor(n^(ratio1))+1):n,ss:(ss+pm[k]-1)]
    ss=ss+pm[k]
  }

  factor_res<-APM(XX, localfactor = TRUE,rmax=rmax,method=method,weight=weight,type = type)

  hatfc<-factor_res$Ghat#estimated common factors
  hatfs<-factor_res$Fhat#estimated specific factors
  hatV<-factor_res$residual#Estimated Idiosyncratic component


  rr_hat=rep(0,length(hatfs)+1)

  if(!anyNA(hatfc)){
    rr_hat[1]=dim(hatfc)[2]
  }else{
    rr_hat[1]=0
  }


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

  if(is.null(rf)) rf=rr_hat

  n_anova=length(c((floor(n^(ratio1))+1):(floor(n^(ratio1)+n^ratio2))))
  hhh=1
  Q_stat=matrix(0,nrow=time_test,ncol=K+1)

  for(iii in 1:time_test){

    P_all<-projection_matrix(FV[hhh:(hhh+n_anova-1),], tol = 1e-8)


    tmp=FV[hhh:(hhh+n_anova-1),(rr_hat[1]+1):dim(FV)[2]]
    P_noc<-projection_matrix(tmp, tol = 1e-8)


    P_diff=P_all-P_noc

    tmp_ind=c((floor(n^(ratio1))+hhh):(floor(n^(ratio1)+n_anova+hhh-1)))
    Q_stat[iii,1]<- matrix(Y1[tmp_ind],nrow=1)%*%P_diff%*%matrix(Y1[tmp_ind],ncol=1) #test statistics
    ss=rr_hat[1]+1
    for(k in 1:K){


      if(rr_hat[k+1]!=0){

        ind1=c(ss:(ss+rr_hat[k+1]-1))

        tmp=FV[hhh:(hhh+n_anova-1),-ind1]
      }else{
        tmp=FV[hhh:(hhh+n_anova-1),]

      }
      P_no_sk<-projection_matrix(tmp, tol = 1e-8)
      P_diff=P_all-P_no_sk
      tmp_ind=c((floor(n^(ratio1))+hhh):(floor(n^(ratio1)+n_anova+hhh-1)))
      Q_stat[iii,k+1]<- matrix(Y1[tmp_ind],nrow=1)%*%P_diff%*%matrix(Y1[tmp_ind],ncol=1) #test statistics

      ss=ss+rr_hat[k+1]
    }
    hhh=hhh+n_anova
  }
  return(list(Q_stat=Q_stat,rf=rf))
}





projection_matrix <- function(A, tol = 1e-8) {
  svd_A <- svd(A)
  U <- svd_A$u
  V <- svd_A$v
  d <- svd_A$d

  d_inv <- ifelse(d > tol, 1/d, 0)
  Sigma_inv <- diag(d_inv, length(d_inv))

  A_pinv <- V %*% Sigma_inv %*% t(U)
  P <- A %*% A_pinv
  return(P)
}

