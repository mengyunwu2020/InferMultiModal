#' Hypothesis Test of the Effects of Idiosyncratic Factors
#'
#' @param X A list of the multi-modality data, each element is a data matrix of each group with dimension n*p_m.
#' @param y A response vector.
#' @param K The number of modalities.
#' @param rmax The maximum factor numbers of all modalities.
#' @param method The method used in the APM algorithm, default is ic, it can also be gap.
#' @param weight The weight of each projection matrix in the APM algorithm, default is TRUE. if weight = FALSE, then simply calculate the mean of all projection matrices.
#' @param type The method used in estimating the factor numbers in each modality initially, default is IC3.
#' @param lambda.seq A user supplied non-negative tuning parameters.
#' @param alpha The significance level used in hypothesis testing. The default is 0.05.
#' @return A list containing the following components:
#'   \describe{
#'     \item{res}{A numeric indicator vector of length \eqn{1 + K + p}.
#'       The first element corresponds to the overall significance test,
#'       the next \eqn{K} elements correspond to the significance results
#'       for each modality, and the remaining \eqn{p} elements correspond
#'       to the significance results for each individual variable.}
#'
#'     \item{pvals}{A numeric vector of p-values with the same length
#'       (\eqn{1 + K + p}) as \code{res}. Each element corresponds to the
#'       respective test in \code{res}.}
#'   }
#' @import ncvreg
#' @import SIS
#' @import Matrix
#' @import GrFA
#' @import glmnet
#' @examples
#' n=300
#' set.seed(2025)
#' K=4
#' p=1200
#' r_f=c(2,rep(2,K))
#' r_all=sum(r_f)
#' gamma1=0.5*rep(1,r_all)
#' pm=rep(p/K,K)
#' beta_0=c(rep(0,p))
#' # ------------------------#beta in the alternative with signal strength w rises from 0.05 to 0.2--------------
#' # beta_1=c(rep(0.05,2),rep(0,p/K-2),rep(0.05,2),rep(0,p/K-2),rep(0.05,2),rep(0,p/K-2),rep(0.05,2),rep(0,p/K-2))
#' # beta_2=c(rep(0.1,2),rep(0,p/K-2),rep(0.1,2),rep(0,p/K-2),rep(0.1,2),rep(0,p/K-2),rep(0.1,2),rep(0,p/K-2))
#' # beta_3=c(rep(0.15,2),rep(0,p/K-2),rep(0.15,2),rep(0,p/K-2),rep(0.15,2),rep(0,p/K-2),rep(0.15,2),rep(0,p/K-2))
#' # beta_4=c(rep(0.2,2),rep(0,p/K-2),rep(0.2,2),rep(0,p/K-2),rep(0.2,2),rep(0,p/K-2),rep(0.2,2),rep(0,p/K-2))
#' library(Matrix)
#' B_c<-matrix(rnorm(p*r_f[1]),nrow=p)
#' B=B_c
#' B_s<-list()
#'   for(k in 1:K){
#'     B_s[[k]]<-matrix(rnorm(pm[k]*r_f[k+1]),nrow=pm[k])
#'   }
#' block_diag <- as.matrix(bdiag(B_s[[1]], B_s[[2]]))
#' block_diag <- as.matrix(bdiag(block_diag, B_s[[3]]))
#' block_diag <- as.matrix(bdiag(block_diag, B_s[[4]]))
#' B=cbind(B_c,as.matrix(block_diag))
#' F<-matrix(rnorm(n*r_all),nrow=n)
#' V<-matrix(rnorm(n*p),nrow=n)
#' X=F%*%t(B)+V
#' ss=1
#' XX=list()
#' for(k in 1:K){
#'   XX[[k]]=X[,ss:(ss+pm[k]-1)]
#'   ss=ss+pm[k]
#' }
#' Y=F%*%gamma1+rnorm(n,0,1)
#' res <- betatest(XX,Y,K)
#' @export
betatest<-function(X,y,K,rmax=8,method='ic',weight=TRUE, type = "IC3",lambda.seq=NULL,alpha=0.05){


  fit=estimation(X,y,K,rmax=rmax,method=method,weight=weight, type =type,lambda.seq = lambda.seq)
  beta_hat=fit$beta
  hatV=fit$V
  Fhat=fit$Fhat
  p=ncol(hatV)
  n=nrow(hatV)

  #-------------------Use Nodewise regression to estimate Theta, the inverse matrix of the population covariance------
  C<-as.matrix(diag(rep(1,p)))
  T<-c()
  for(j in 1:ncol(hatV)){
    fit_u = glmnet(hatV[,-j], hatV[,j], intercept=FALSE,
                   lambda=cv.glmnet(hatV[,-j], hatV[,j],intercept=FALSE)$lambda.1se)
    beta<-as.vector(fit_u$beta)
    C[j,-j]<--beta
    T<-c(T,1/n*sum((hatV[,j]-hatV[,-j]%*%beta)^2)+fit_u$lambda/2*sum(abs(beta)))
    if (j%%10==0){
      # print(j)
    }
  }
  T1<-diag(1/T)
  Theta<-T1%*%C #estimated Theta
  #------------------Conduct bootstrap to compute the critical values----------

    boostrap<-c()
    boostres=NULL
    boostres2=NULL
    for(j in 1:2500){
      tmp_a=abs(Theta%*%t(hatV)%*%rnorm(n,mean=0,sd=1))
      hs=rep(0,K)
      ss=1
      for(k in 1:K){
        hs[k]=max(tmp_a[ss:(ss+ncol(X[[k]])-1),1])
        ss=ss+ncol(X[[k]])
        }
      boostres=cbind(boostres,hs)#modality
      boostres2=cbind(boostres2,tmp_a)#each
      boostrap<-c(boostrap,max(tmp_a))#all
    }



   XX=X[[1]]
  pm=rep(0,K)
  pm[1]=ncol(X[[1]])
  for(kk in 2:K){
    XX=cbind(XX,X[[kk]])
    pm[kk]=ncol(X[[kk]])
  }
  sd_alter<-rcv(y,XX,n,pm,K,rmax=rmax,method=method,weight=weight,lambda.seq = lambda.seq,type=type)
  d_beta<-beta_hat+1/n*Theta%*%t(hatV)%*%(y-hatV%*%beta_hat-Fhat%*%fit$gamma)


  c_alpha<-1/sqrt(n)*quantile(boostrap,(1-alpha))
  res_all=sum(sqrt(n)*max(abs(d_beta))/sqrt(sd_alter[1])>c_alpha)


  tmp= apply(boostres,1,function(x)  quantile(x,(1-alpha)))
  c_alpha_modality<-1/sqrt(n)*tmp
  ss=1
  res_modality=rep(0,K)
  for(kk in 1:K){
    res_modality[kk]=sum(sqrt(n)*max(abs(d_beta[ss:(ss+pm[kk]-1)]))/sqrt(sd_alter[1])>c_alpha_modality[kk])
    ss=ss+pm[kk]
  }
  res=c(res_all,res_modality)

  tmp= apply(boostres2,1,function(x)  quantile(x,(1-alpha)))
  c_alpha_each<-1/sqrt(n)*tmp
  res_each=rep(0,p)
  for(kk in 1:p){
    res_each[kk]=sum(sqrt(n)*(abs(d_beta[kk]))/sqrt(sd_alter[1])>c_alpha_each[kk])
  }
  res=c(res,res_each)
  names(res)=c('all',paste0('modality',1:K),paste0('variable',1:p))

  p_value_all<-sum(1/sqrt(n)*boostrap>((sqrt(n)*max(abs(d_beta))/sqrt(sd_alter[1]))))/2500 #Compute p-values
  ss=1
  p_value_modality=NULL
  for(kk in 1:K){
    p_value_modality<-c(p_value_modality,sum(1/sqrt(n)*boostres[kk,]>((sqrt(n)*max(abs(d_beta[ss:(ss+pm[kk]-1)]))/sqrt(sd_alter[1]))))/2500) #Compute p-values
    ss=ss+pm[kk]
  }
  pvals=c(p_value_all,p_value_modality)
  p_value_each=NULL
  for(kk in 1:p){
    tmp=sum(1/sqrt(n)*boostres2[kk,]>((sqrt(n)*abs(d_beta[kk])/sqrt(sd_alter[1]))))/2500
    p_value_each<-c(p_value_each,tmp)
  }
  pvals=c(pvals,p_value_each)
  names(pvals)=c('all',paste0('modality',1:K),paste0('variable',1:p))





  return(list(res=res,pvals=pvals))


}




