#' Hypothesis Test of the Effects of Idiosyncratic Residuals
#'
#' @param X A list of the multi-modality data, each element is a data matrix of each group with dimension n*p_k.
#' @param y A response vector.
#' @param K The number of modalities.
#' @param fit An object returned by [estimation()].
#' @param rmax The maximum factor numbers of all modalities.
#' @param r0 The number of common factors. Default is NULL, the algorithm will automatically estimate the number of global factors. If you have prior information about the true number of global factors, you can set it manually.
#' @param r The number of local factors in each group. Default is NULL, the algorithm will automatically estimate the number of local factors. If you have prior information, set it manually as an integer vector of length K (the number of groups).
#' @param method The method used in the APM algorithm, default is gap, it can also be ic.
#' @param weight The weight of each projection matrix in the APM algorithm, default is TRUE. if weight = FALSE, then simply calculate the mean of all projection matrices.
#' @param type The method used in estimating the factor numbers in each modality initially, default is BIC3.
#' @param lambda.seq A user supplied non-negative tuning parameters.
#' @param alpha The significance level used in hypothesis testing. The default is 0.05.
#' @param factor_res A list containing the following components:
#' \describe{
#' \item{Ghat}{The supplied common factors. A numeric matrix or vector.}
#' \item{loading_G}{A list consisting of the common factor loadings. Each element of the list corresponds to one common factor.}
#' \item{Fhat}{The supplied specific factors. A numeric matrix or vector.}
#' \item{loading_F}{A list consisting of the specific factor loadings. Each list element corresponds to one specific factor.}
#' \item{residual}{A list consisting of the residuals.}
#' }
#' @param boot_rep Integer. The number of bootstrap repetitions (outer loop).
#' @param boot_each Integer. The number of samples per repetition (inner loop).
#' @param gamma_MCP The tuning parameter of the MCP penalty. Default is 3 for MCP.
#' @param seed_sis An integer specifying the random seed used by the SIS method for noise variance estimation.
#' @param tune_sis Method for tuning the regularization parameter of the penalized likelihood
#' subproblems and the final model in the SIS noise variance estimation step.
#' Options include \code{"bic"}, \code{"ebic"}, \code{"aic"}, and \code{"CV"}.
#' @param tun_ncvreg Character string specifying the tuning method for the regression step. Options are \code{"BIC"} or \code{"CV"}.
#'
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
#' res_fac = estimation_fac(XX, K, nodewise_step=TRUE)
#' fit = estimation(XX,Y, K,res_fac =res_fac)
#' res <- betatest(XX,Y,K,fit=fit)
#' @export
betatest <- function(X, y, K, fit, rmax=10, r0=NULL, r=NULL, method='gap', weight=TRUE, type="BIC3", lambda.seq=NULL, alpha=0.05, factor_res=NULL,boot_rep=5,boot_each=500,gamma_MCP=3,seed_sis=1,tune_sis='ebic',tun_ncvreg='CV'){

  beta_hat = fit$beta; hatV = fit$V; Fhat = fit$Fhat; p = ncol(hatV); n = nrow(hatV)

  list2env(fit$res_node, envir = environment())

  group_idx_list <- list()
  start <- 1
  for(k in 1:K){
    group_idx_list[[k]] <- start:(start + ncol(X[[k]]) - 1)
    start <- start + ncol(X[[k]])
  }




  XX = do.call(cbind, X); pm = sapply(X, ncol)
  sd_alter <- rcv(y, XX, n, pm, K, rmax=rmax, r0=r0, r=r, method=method, weight=weight, lambda.seq=lambda.seq, type=type, factor_res=factor_res,gamma_MCP=gamma_MCP,seed_sis=seed_sis,tune_sis =tune_sis,tun_ncvreg=tun_ncvreg)
  d_beta <- beta_hat + 1/n * theta_tv %*% (y - hatV %*% beta_hat - Fhat %*% fit$gamma)

  res_all = sum(sqrt(n)*max(abs(d_beta))/sqrt(sd_alter[1]) > c_alpha)
  res_modality = sapply(1:K, function(kk) sum(sqrt(n)*max(abs(d_beta[group_idx_list[[kk]]]))/sqrt(sd_alter[1]) > c_alpha_mod[kk]))

  res_each = sapply(1:p, function(kk) sum(sqrt(n)*(abs(d_beta[kk]))/sqrt(sd_alter[1]) > c_alpha_each[kk]))

  p_val_all <- mean(all_samples/sqrt(n) > (sqrt(n)*max(abs(d_beta))/sqrt(sd_alter[1])))
  p_val_mod <- sapply(1:K, function(kk) mean(mod_samples[kk,]/sqrt(n) > (sqrt(n)*max(abs(d_beta[group_idx_list[[kk]]]))/sqrt(sd_alter[1]))))
  p_val_each <- sapply(1:p, function(kk) mean(each_samples[kk,]/sqrt(n) > (sqrt(n)*abs(d_beta[kk])/sqrt(sd_alter[1]))))

  res = c(res_all, res_modality, res_each)
  pvals = c(p_val_all, p_val_mod, p_val_each)
  names(res) = names(pvals) = c('all', paste0('modality',1:K), paste0('variable',1:p))

  return(list(res=res, pvals=pvals,fit = fit))
}

