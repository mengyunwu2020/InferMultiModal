#' Estimation function for the proposed model.
#'
#' @param X A list of the multi-modality data, each element is a data matrix of each group with dimension n*p_k.
#' @param y A response vector.
#' @param K The number of modalities.
#' @param rmax The maximum factor numbers of all modalities.
#' @param r0 The number of common factors. Default is NULL, the algorithm will automatically estimate the number of global factors. If you have prior information about the true number of global factors, you can set it manually.
#' @param r The number of local factors in each group. Default is NULL, the algorithm will automatically estimate the number of local factors. If you have prior information, set it manually as an integer vector of length K (the number of groups).
#' @param method The method used in the APM algorithm, default is gap, it can also be ic.
#' @param weight The weight of each projection matrix in the APM algorithm, default is TRUE. if weight = FALSE, then simply calculate the mean of all projection matrices.
#' @param type The method used in estimating the factor numbers in each modality initially, default is BIC3.
#' @param lambda.seq A user supplied non-negative tuning parameters.
#' @param factor_res A list containing the following components:
#' \describe{
#' \item{Ghat}{The supplied common factors. A numeric matrix or vector.}
#' \item{loading_G}{A list consisting of the common factor loadings. Each element of the list corresponds to one common factor.}
#' \item{Fhat}{The supplied specific factors. A numeric matrix or vector.}
#' \item{loading_F}{A list consisting of the specific factor loadings. Each list element corresponds to one specific factor.}
#' \item{residual}{A list consisting of the residuals.}
#' }
#' @param gamma_MCP The tuning parameter of the MCP penalty. Default is 3 for MCP.
#' @param res_fac A list (or object) returned by \code{estimation_fac},
#'   containing the estimated factor components.
#' @param tun_ncvreg Character string specifying the tuning method for the regression step. Options are \code{"BIC"} or \code{"CV"}.
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
#' r<-0.3
#' n=ceiling(log(p)/(r^2/4^2))
#' dis=NULL
#' rr_hat=NULL
#' F<-matrix(rnorm(n*r_all),nrow=n)
#' V<-matrix(rnorm(n*p),nrow=n)
#' X=F%*%t(B)+V
#' Y=F%*%gamma1+V%*%beta+rnorm(n,mean=0,sd=1)
#' ss=1
#' XX=list()
#' for(k in 1:K){
#'   XX[[k]]=X[,ss:(ss+pm[k]-1)]
#'   ss=ss+pm[k]
#' }
#' res_fac<-estimation_fac(XX,K,method='gap')
#' fit<-estimation(XX,Y,K,res_fac=res_fac)
#' @export
#'
#'
estimation<-function(X,y,K,rmax=10,r0=NULL,r=NULL,method='gap',weight=TRUE, type = "BIC3",lambda.seq=NULL,factor_res=NULL,res_fac,gamma_MCP=3,tun_ncvreg='CV'){


  list2env(res_fac, envir = environment())
  penalty_factors=c(rep(0,rhat),rep(1,p))

  if(tun_ncvreg=='BIC'){
    if(is.null(lambda.seq)){
      cfit <- ncvreg(Re(FV),y,penalty.factor = penalty_factors, penalty="MCP",gamma=gamma_MCP)
    }else{
      cfit <- ncvreg(Re(FV),y,penalty.factor = penalty_factors, penalty="MCP", lambda = lambda.seq,gamma=gamma_MCP)

    }
    bic_values <- BIC(cfit)

    best_idx <- which.min(bic_values)
    cv_res <- cfit$beta[, best_idx]

  }else{
    set.seed(5432)
    if(is.null(lambda.seq)){
      cvfit <- cv.ncvreg(Re(FV), y,penalty.factor = penalty_factors, penalty="MCP",gamma=gamma_MCP)

    }else{
      cvfit <- cv.ncvreg(Re(FV), y,penalty.factor = penalty_factors, penalty="MCP",lambda = lambda.seq,gamma=gamma_MCP)
    }

    fit <- cvfit$fit
    cv_res<-fit$beta[,cvfit$min]
  }


  res<-list(Fhat=Fhat,V=V,FV=FV,gamma=cv_res[2:(rhat+1)],beta=cv_res[(rhat+2):(rhat+p+1)],rr_hat=rr_hat,loading_c=loading_c,loading_s=loading_s,res_node=res_node)
  class(res) <- "estimation"
  return(res)
}

#' @export
estimation_fac<-function(X,K,rmax=10,r0=NULL,r=NULL,method='gap',weight=TRUE, type = "BIC3",factor_res=NULL,nodewise_step=FALSE,boot_rep=5,boot_each=500,alpha=0.05){
  if(is.null(factor_res)){
    factor_res<-APM(X, rmax = rmax,r0=r0,r=r,localfactor = TRUE, method = method,weight=weight, type = type)

  }
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
    rr_hat[1]=ncol(fc)
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
  rhat=dim(FV)[2]
  p=dim(V)[2]
  Fhat=FV
  FV=cbind(FV,V)



  res_node=NULL
  if(nodewise_step){
    group_idx_list <- list()
    start <- 1
    for(k in 1:K){
      group_idx_list[[k]] <- start:(start + ncol(X[[k]]) - 1)
      start <- start + ncol(X[[k]])
    }
    res_node=nodewise_step1(V,boot_rep=boot_rep,boot_each=boot_each,K=K,group_idx_list=group_idx_list,alpha=alpha)
  }

  res<-list(Fhat=Fhat,V=V,FV=FV,rr_hat=rr_hat,loading_c=loading_c,loading_s=loading_s,res_node=res_node,rhat=rhat,p=p)

  return(res)
}

nodewise_step1<-function(hatV,boot_rep=5,boot_each=500,K,group_idx_list,alpha){
  p=ncol(hatV)
  sds <- apply(hatV, 2, sd)
  n=nrow(hatV)
  C <- as.matrix(diag(rep(1, p)))
  T_vec <- c()
  for(j in 1:p){
    fit_u = glmnet(hatV[,-j], hatV[,j], intercept=FALSE, lambda=cv.glmnet(hatV[,-j], hatV[,j], intercept=FALSE)$lambda.1se)
    beta_u <- as.vector(fit_u$beta)

    min_tau_sq <- 0.0001 * mean(sds^2)


    C[j,-j] <- -beta_u
    rsd_tmp= max((1/n*sum((hatV[,j]-hatV[,-j]%*%beta_u)^2)+ fit_u$lambda/2*sum(abs(beta_u))),min_tau_sq)
    T_vec <- c(T_vec,rsd_tmp)
  }
  Theta <- diag(1/T_vec) %*% C

  all_samples <- c(); mod_samples <- matrix(NA, nrow=K, ncol=0); each_samples <- matrix(NA, nrow=p, ncol=0)
  quant_vec <- c(); quants_mod <- NULL; quants_each <- NULL
  theta_tv=Theta %*% t(hatV)


  set.seed(1)

  for(m in 1:boot_rep){
    b_all <- c(); b_mod <- NULL; b_each <- NULL
    for(j in 1:boot_each){
      tmp_a = abs( theta_tv %*% rnorm(n, mean=0, sd=1))
      hs = sapply(1:K, function(k) max(tmp_a[group_idx_list[[k]], 1]))
      b_mod = cbind(b_mod, hs); b_each = cbind(b_each, tmp_a); b_all = c(b_all, max(tmp_a))
    }

    quant_vec <- c(quant_vec, quantile(b_all, 1-alpha))
    quants_mod <- cbind(quants_mod, apply(b_mod, 1, function(x) quantile(x, 1-alpha)))
    quants_each <- cbind(quants_each, apply(b_each, 1, function(x) quantile(x, 1-alpha)))

    all_samples <- c(all_samples, b_all)
    mod_samples <- cbind(mod_samples, b_mod)
    each_samples <- cbind(each_samples, b_each)
  }

  c_alpha <- 1/sqrt(n) * mean(quant_vec)
  c_alpha_mod <- 1/sqrt(n) * rowMeans(quants_mod)
  c_alpha_each <- 1/sqrt(n) * rowMeans(quants_each)
  return(list(c_alpha=c_alpha,c_alpha_mod=c_alpha_mod,c_alpha_each=c_alpha_each,all_samples=all_samples,mod_samples=mod_samples,each_samples=each_samples,theta_tv=theta_tv))
}

