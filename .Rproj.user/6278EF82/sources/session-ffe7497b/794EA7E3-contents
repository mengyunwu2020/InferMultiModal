
library(MASS)
get_idx <- function(k, X) {
  col_counts <- sapply(X, ncol)
  start <- sum(col_counts[1:(k-1)]) + 1
  end <- sum(col_counts[1:k])
  return(start:end)
}


prediction_joint<-function(X_test,Y_test,K,fit,maxitier=50,eps=10^(-6),Vhat_update=TRUE,Vhat_imput){
  beta_hat1<-as.vector(fit$beta)
  gamma_hat1<-as.vector(fit$gamma)
  names(gamma_hat1)=rep(0:K,fit$rr_hat)
  loading_c=fit$loading_c
  loading_s=fit$loading_s
  loading_c_mat=loading_c[[1]]
  loading_s_mat=loading_s[[1]]
  for(k in 1:(K-1)){
    loading_c_mat=rbind(loading_c_mat,loading_c[[k+1]])
    loading_s_mat=as.matrix(bdiag(loading_s_mat,loading_s[[k+1]]))
  }

  loading_mat=cbind(loading_c_mat,loading_s_mat)
  projection <- ginv(loading_mat)
  Fhat <- X_test %*% t(projection)


  if(Vhat_update){

  Vhat <- X_test - Fhat%*%t(loading_mat)

  }else{
    Vhat=Vhat_imput
  }


  y_predict=Vhat%*%beta_hat1+Fhat%*%gamma_hat1






  Y_c=Y_test

  for(i in 1:maxitier){

    Y_c=Y_c-Vhat%*%beta_hat1
    XY_c=cbind(X_test,Y_c)


    loading_mat_g=rbind(loading_mat,gamma_hat1)
    projection <- ginv(loading_mat_g)
    Fhat_p=Fhat
    Fhat <- XY_c %*% t(projection)


    if(Vhat_update){

      Vhat=  X_test - Fhat%*%t(loading_mat)
    }



    if(sum((Fhat_p-Fhat)^2)<eps) break;


  }


  y_predict_joint=Vhat%*%beta_hat1+Fhat%*%gamma_hat1



  res<-list(y_predict_joint=y_predict_joint, y_predict= y_predict)
  return(res)


}


joint_func<-function(X,y,K,rmax=10,method='gap',weight=TRUE, type = "BIC3",lambda, maxitier=50,eps=10^(-6),fit,r_update=TRUE,Vhat_update=TRUE,Vhat_imput,beta_imput){

   BIC_score=10e+10

  if(Vhat_update){

    Vhat=fit[["FV"]][,-c(1:sum(fit[["rr_hat"]]))]
  }else{
    Vhat=Vhat_imput
  }

  Fhat=fit[["Fhat"]]
  n=nrow(Fhat)
  Fc=Fhat[,1:fit$rr_hat[1]]
  Fs=Fhat[,(1+fit$rr_hat[1]):(sum(fit$rr_hat))]

  if(Vhat_update){
    betahat=fit[["beta"]]
  }else{
    betahat=beta_imput
  }

  gammahat_s=(fit$gamma)[(1+fit$rr_hat[1]):(sum(fit$rr_hat))]
  gammahat_c=(fit$gamma)[1:(fit$rr_hat[1])]

  Fs_l=gammahat_s_l=list()
  tmp=1
  for(k in 1:K){
    Fs_l[[k]]=Fs[,(tmp):(tmp+(fit$rr_hat)[(k+1)]-1)]
    gammahat_s_l[[k]]=gammahat_s[(tmp):(tmp+(fit$rr_hat)[(k+1)]-1)]
    tmp=tmp+(fit$rr_hat)[(k+1)]
  }


  loading_c=fit$loading_c
  loading_s=fit$loading_s
  loading_c_mat=loading_c[[1]]
  loading_s_mat=loading_s[[1]]

  for(k in 1:(K-1)){
    if(length(loading_c)==K){
      loading_c_mat=rbind(loading_c_mat,loading_c[[k+1]])
    }
    if(!(anyNA(loading_s[[k+1]]))){
      loading_s_mat=as.matrix(bdiag(loading_s_mat,loading_s[[k+1]]))
    }
  }


  tmp_c=Fs%*%t(loading_s_mat)
  Y_c=y-Fs%*%gammahat_s
  pm <- sapply(X, ncol)
  XX=do.call(cbind,X)
  p=ncol(XX)
  if(r_update){
    for(i in 1:maxitier){

      Y_c=Y_c-Vhat%*%betahat
      X_c=XX- tmp_c
      XY_c=cbind(X_c,Y_c)

      cov<-XY_c%*%t(XY_c)/n
      ee=eigen(cov)

      eigvec=ee$vectors
      eigval=ee$values


      K_est<-which.min(diff(log(eigval[1:10]))) #estimate factor numbers using eigenvalue ratios.

      Fc_old=Fc
      Fc <- Re(sqrt(n) * eigvec[, 1:K_est, drop = FALSE])

      hatB2<-t(1/n*t(Fc)%*%XY_c) #Estimated Factor Loading
      loading_c_mat=hatB2[1:p,]
      gamma_c_hat=hatB2[(p+1),]

      Y_c=Y_c-Fc%*%gamma_c_hat
      tmp=1
      for(k in 1:K){
        Y_c=Y_c+matrix(Fs_l[[k]],ncol=length(gammahat_s_l[[k]]))%*%gammahat_s_l[[k]]
        XY_c=cbind((X[[k]]-Fc%*%t(loading_c_mat[tmp:(tmp+pm[k]-1),])),Y_c)
        tmp=tmp+pm[k]
        cov<-XY_c%*%t(XY_c)/n

        ee=eigen(cov)

        eigvec=ee$vectors
        eigval=ee$values


        K_est<-which.min(diff(log(eigval[1:10]))) #estimate factor numbers using eigenvalue ratios.

        Fs_l[[k]]<- Re(sqrt(n) * eigvec[, 1:K_est, drop = FALSE])

        hatB2<-t(1/n*t(Fs_l[[k]])%*%XY_c) #Estimated Factor Loading
        loading_s[[k]]=hatB2[1:(pm[k]),]
        gammahat_s_l[[k]]=hatB2[(pm[k]+1),]
        Y_c=Y_c-matrix(Fs_l[[k]],ncol=K_est)%*%gammahat_s_l[[k]]
      }

      loading_s_mat=loading_s[[1]]

      for(k in 1:(K-1)){
        if(!(anyNA(loading_s[[k+1]]))){
          loading_s_mat=as.matrix(bdiag(loading_s_mat,loading_s[[k+1]]))
        }
      }


      tmp_c=NULL
      for(k in 1:K){
        tmp_c=cbind(tmp_c,Fs_l[[k]]%*%t(loading_s[[k]]))
      }


      if(Vhat_update){

        Vhat=XX-Fc%*%t(loading_c_mat)-tmp_c
      }
      Y_c=Y_c+Vhat%*%betahat
      if(Vhat_update){
        fit3=ncvreg(Re(Vhat), Y_c, penalty="MCP",lambda=ifelse(length(lambda)>1,lambda,c(lambda,lambda)))
        betaold=betahat
        betahat <- coef(fit3, which = which.min(BIC(fit3)))[-1]
        BIC_score=min(BIC(fit3))

      }

      Y_c=Y_c+Fc%*%gamma_c_hat
      if(sum((betahat-betaold)^2)<eps) break;


    }
  }else{
    for(i in 1:maxitier){

      Y_c=Y_c-Vhat%*%betahat
      X_c=XX- tmp_c
      XY_c=cbind(X_c,Y_c)

      cov<-XY_c%*%t(XY_c)/n

      ee=eigen(cov)
      eigvec=ee$vectors
      eigval=ee$values

      K_est<-fit$rr_hat[1]
      Fc_old=Fc
      Fc<-Re(sqrt(n)*eigvec[,1:K_est]) #estimated factors

      hatB2<-t(1/n*t(Fc)%*%XY_c) #Estimated Factor Loading
      loading_c_mat=hatB2[1:p,]
      gamma_c_hat=hatB2[(p+1),]

      Y_c=Y_c-Fc%*%gamma_c_hat
      tmp=1
      for(k in 1:K){
        Y_c=Y_c+Fs_l[[k]]%*%gammahat_s_l[[k]]
        XY_c=cbind((X[[k]]-Fc%*%t(loading_c_mat[tmp:(tmp+pm[k]-1),])),Y_c)
        tmp=tmp+pm[k]
        cov<-XY_c%*%t(XY_c)/n
        ee=eigen(cov)

        eigvec=ee$vectors
        eigval=ee$values

        K_est<-fit$rr_hat[(k+1)]

        Fs_l[[k]]<-Re(sqrt(n)*eigvec[,1:K_est]) #estimated factors

        hatB2<-t(1/n*t(Fs_l[[k]])%*%XY_c) #Estimated Factor Loading
        loading_s[[k]]=hatB2[1:(pm[k]),]
        gammahat_s_l[[k]]=hatB2[(pm[k]+1),]
        Y_c=Y_c-matrix(Fs_l[[k]],ncol=K_est)%*%gammahat_s_l[[k]]
      }

      loading_s_mat=loading_s[[1]]

      for(k in 1:(K-1)){
        if(!(anyNA(loading_s[[k+1]]))){
          loading_s_mat=as.matrix(bdiag(loading_s_mat,loading_s[[k+1]]))
        }
      }


      tmp_c=NULL
      for(k in 1:K){
        tmp_c=cbind(tmp_c,Fs_l[[k]]%*%t(loading_s[[k]]))
      }


      if(Vhat_update){

        Vhat=XX-Fc%*%t(loading_c_mat)-tmp_c
      }
      Y_c=Y_c+Vhat%*%betahat
      if(Vhat_update){
        fit3=ncvreg(Re(Vhat), Y_c, penalty="MCP",lambda=ifelse(length(lambda)>1,lambda,c(lambda,lambda)))
        betaold=betahat
        betahat <- coef(fit3, which = which.min(BIC(fit3)))[-1]
        BIC_score=min(BIC(fit3))
      }
      Y_c=Y_c+Fc%*%gamma_c_hat
      if(sum((Fc-Fc_old)^2)<eps) break;


    }
  }




  if(!anyNA(Fc)){
    FV=Fc
  }else{
    FV=NULL
  }
  rr_hat=rep(0,length(Fs_l)+1)
  if(!anyNA(Fc)){
    rr_hat[1]=dim(Fc)[2]
  }else{
    rr_hat[1]=0
  }
  for(i in 1:length(Fs_l)){
    if(!anyNA(Fs_l[[i]])){
      FV=cbind(FV,Fs_l[[i]])
      rr_hat[i+1]=dim(Fs_l[[i]])[2]
    }
  }
  r=dim(FV)[2]
  Fhat=FV
  FV=cbind(FV,Vhat)

  tmp=1
  for(k in 1:K){
    loading_c[[k]]= loading_c_mat[tmp:(tmp+pm[k]-1),]
    tmp=tmp+pm[k]

  }
  res<-list(fit_ours=fit,Fhat=Fhat,V=Vhat,FV=FV,gamma=c(gamma_c_hat,unlist(gammahat_s_l)),beta=betahat,rr_hat=rr_hat,loading_c=loading_c,loading_s=loading_s, BIC_score= BIC_score)

  return(res)


}


betatest_joint<-function(X,y,K,rmax=10,r0=NULL,r=NULL,method='gap',weight=TRUE, type = "BIC3",lambda=c(0.1,0.05),alpha=0.05,factor_res=NULL,fit,r_update=TRUE,boot_rep=5,boot_each=500,lam_in_iteration){

  tmp=1e+100
  if(lam_in_iteration){
    fit=joint_func(X,y,K,rmax=rmax,method=method,weight=weight, type =type,lambda = lambda,fit=fit,r_update=r_update)
  }else{
    for(i in 1:length(lambda)){
      res=joint_func(X,y,K,rmax=rmax,method=method,weight=weight, type =type,lambda = lambda[i],fit=fit,r_update=r_update)
      if(tmp>res$BIC_score){
        fit=res
        tmp=res$BIC_score
      }
    }
  }

  beta_hat = fit$beta; hatV = fit$V; Fhat = fit$Fhat; p = ncol(hatV); n = nrow(hatV)

  C <- as.matrix(diag(rep(1, p)))
  T_vec <- c()
  for(j in 1:p){
    fit_u = glmnet(hatV[,-j], hatV[,j], intercept=FALSE, lambda=cv.glmnet(hatV[,-j], hatV[,j], intercept=FALSE)$lambda.1se)
    beta_u <- as.vector(fit_u$beta)
    C[j,-j] <- -beta_u
    T_vec <- c(T_vec, 1/n*sum((hatV[,j]-hatV[,-j]%*%beta_u)^2) + fit_u$lambda/2*sum(abs(beta_u)))
  }
  Theta <- diag(1/T_vec) %*% C

  all_samples <- c(); mod_samples <- matrix(NA, nrow=K, ncol=0); each_samples <- matrix(NA, nrow=p, ncol=0)
  quant_vec <- c(); quants_mod <- NULL; quants_each <- NULL
  theta_tv=Theta %*% t(hatV)
  for(m in 1:boot_rep){
    b_all <- c(); b_mod <- NULL; b_each <- NULL
    for(j in 1:boot_each){
      tmp_a = abs(theta_tv %*% rnorm(n, mean=0, sd=1))
      hs = sapply(1:K, function(k) max(tmp_a[get_idx(k, X), 1]))
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

  XX = do.call(cbind, X); pm = sapply(X, ncol)

  sd_alter<-rcv4(y,XX,n,pm,K,rmax=rmax,method=method,weight=weight, type =type,lambda=lambda,factor_res=NULL,r_update=r_update,lam_in_iteration=lam_in_iteration)


  d_beta <- beta_hat + 1/n * theta_tv %*% (y - hatV %*% beta_hat - Fhat %*% fit$gamma)

  res_all = sum(sqrt(n)*max(abs(d_beta))/sqrt(sd_alter[1]) > c_alpha)
  res_modality = sapply(1:K, function(kk) sum(sqrt(n)*max(abs(d_beta[get_idx(kk, X)]))/sqrt(sd_alter[1]) > c_alpha_mod[kk]))
  res_each = sapply(1:p, function(kk) sum(sqrt(n)*(abs(d_beta[kk]))/sqrt(sd_alter[1]) > c_alpha_each[kk]))

  p_val_all <- mean(all_samples/sqrt(n) > (sqrt(n)*max(abs(d_beta))/sqrt(sd_alter[1])))
  p_val_mod <- sapply(1:K, function(kk) mean(mod_samples[kk,]/sqrt(n) > (sqrt(n)*max(abs(d_beta[get_idx(kk, X)]))/sqrt(sd_alter[1]))))
  p_val_each <- sapply(1:p, function(kk) mean(each_samples[kk,]/sqrt(n) > (sqrt(n)*abs(d_beta[kk])/sqrt(sd_alter[1]))))

  res = c(res_all, res_modality, res_each)
  pvals = c(p_val_all, p_val_mod, p_val_each)
  names(res) = names(pvals) = c('all', paste0('modality',1:K), paste0('variable',1:p))

  return(list(res=res, pvals=pvals))




}



gtest_joint<-function(X,y,K,ratio1,ratio2,rmax=10,r0=NULL,r=NULL,method='gap',weight=TRUE,type="BIC3",lambda.seq=NULL,r_f=NULL,alpha=0.05,factor_res=NULL,lambda,r_update=TRUE,lam_in_iteration){


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
        tmp <- capture.output(out <-rcv4(y,X_1,n,pm,K,rmax=rmax,method=method,weight=weight, type =type,lambda=lambda,factor_res=NULL,r_update=r_update,lam_in_iteration=lam_in_iteration), file = NULL)
        out
      }
    )),
    silent = TRUE
  ) #estimate variance through function rcv
  if (!inherits(sig_est, "try-error")) {
  var_est1<-sig_est[1] #estimated variance via isis.
  }else{
    var_est1<-1
  }
  output11 <- try(
    suppressMessages(suppressWarnings(
      {
        tmp <- capture.output(
          out <- selectfeature4(y, X_1, n, K, pm, ratio1 = ratio1, ratio2 = ratio2,
                                rf = r_f, rmax = rmax, method = method, weight = weight,
                                lambda = lambda, alpha = alpha, type = type,
                                r_update = r_update,lam_in_iteration=lam_in_iteration),
          file = NULL
        )
        out
      }
    )),
    silent = TRUE
  )

  if (!inherits(output11, "try-error")) {
    output1 <- output11$Q_stat
    p=sum(pm)

    pvals=rep(0,(K+1))

    r_f=output11$rf
    for(kk in 1:(K+1)){
      tmp=pchisq(q = output1[,kk]/var_est1, df = r_f[kk], lower.tail = FALSE)

      maxp=ifelse(length(tmp)==1,tmp,maximump(tmp)$p)

      restmp=c(ACAT(tmp),maxp)


      pvals[kk]=mean(restmp)


    }

    names(pvals)=c('common',paste0('specific',1:K))
    return(pvals)

  }else {

    pvals=rep(NA,(K+1))
    return(pvals)
  }


}


selectfeature4<-function(Y1,X,n,K,pm,ratio1=0.8,ratio2=0.9,rf,rmax=10,r0=NULL,r=NULL,method='gap',weight=TRUE,lambda,alpha=0.05, type = "BIC3",factor_res=NULL,r_update=TRUE,lam_in_iteration){

  time_test= floor((n-floor(n^(ratio1)))/n^ratio2)
  ss=1
  XX=list()
  for(k in 1:K){
    XX[[k]]=X[1:floor(n^(ratio1)),ss:(ss+pm[k]-1)]
    ss=ss+pm[k]
  }
  p=sum(pm)


  if(!(is.null(factor_res))){
    factor_res1=list()
    factor_res1$loading_c=factor_res$loading_G
    factor_res1$loading_s=factor_res$loading_F
    Idx=1:floor(n^(ratio1))
    factor_res1$fc<-(factor_res$Ghat)[Idx,]
    # factor_res1$fs<-(factor_res$Fhat)[Idx,]

    factor_res1$fs=list()
    for(k in 1:K){

      factor_res1$fs[[k]]<-(factor_res$Fhat)[[k]][Idx,]
    }


    factor_res1$v<-(factor_res$residual)[Idx,]

  }else{
    factor_res1=NULL
  }

  tmp=1e+100

  res_fac = estimation_fac(XX, K,rmax=rmax,method=method,weight=weight,type=type,factor_res=NULL,nodewise_step=FALSE)

  fit=estimation(XX,Y1[1:floor(n^(ratio1))],K,rmax=rmax,method=method,weight=weight,type=type,lambda.seq=NULL,res_fac =res_fac)
  if(lam_in_iteration){

    fit1=joint_func(XX,Y1[1:floor(n^(ratio1))],K,rmax=rmax,method=method,weight=weight, type =type,lambda = lambda,fit=fit,r_update=r_update)

  }else{
    for(i in 1:length(lambda)){

    res=joint_func(XX,Y1[1:floor(n^(ratio1))],K,rmax=rmax,method=method,weight=weight, type =type,lambda = lambda[i],fit=fit,r_update=r_update)
    if(tmp>res$BIC_score){
      fit1=res
      tmp=res$BIC_score
    }
  }
  }


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

  if(!(is.null(factor_res))){
    factor_res1=list()
    factor_res1$loading_c=factor_res$loading_G
    factor_res1$loading_s=factor_res$loading_F
    Idx=(floor(n^(ratio1))+1):n
    factor_res1$fc<-(factor_res$Ghat)[Idx,]
    factor_res1$fs=list()
    for(k in 1:K){

      factor_res1$fs[[k]]<-(factor_res$Fhat)[[k]][Idx,]
    }
    factor_res1$v<-(factor_res$residual)[Idx,]
    factor_res=factor_res1

  }else{


    res_fac = estimation_fac(XX, K,rmax=rmax,method=method,weight=weight,type=type,factor_res=NULL,nodewise_step=FALSE)

    fit=estimation(XX,Y1[(floor(n^(ratio1))+1):n],K,rmax=rmax,method=method,weight=weight,type=type,lambda.seq=NULL,res_fac =res_fac)
    if(lam_in_iteration){

      fit1=joint_func(XX,Y1[(floor(n^(ratio1))+1):n],K,rmax=rmax,method=method,weight=weight, type =type,lambda = lambda,fit=fit,r_update=r_update)

    }else{
      tmp=1e+100
      for(i in 1:length(lambda)){

        res=joint_func(XX,Y1[(floor(n^(ratio1))+1):n],K,rmax=rmax,method=method,weight=weight, type =type,lambda = lambda[i],fit=fit,r_update=r_update)
        if(tmp>res$BIC_score){
          fit1=res
          tmp=res$BIC_score
        }
      }
    }
  }

  hatfc<-fit1$Fhat[,1:fit1$rr_hat[1]]#estimated common factors
  fs<-fit1$Fhat[,-c(1:fit1$rr_hat[1])]#estimated specific factors
  V<-fit1$V#Estimated Idiosyncratic component


  rr_hat=fit1$rr_hat




  FV=cbind(fit1$Fhat,V[,select1])

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
rcv4<-function(Y,X,n,pm,K,rmax=10,method='gap',weight=TRUE, type = "BIC3",lambda = NULL,factor_res=NULL,r_update=TRUE,lam_in_iteration){
  p=ncol(X)
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



  if(!(is.null(factor_res))){
    factor_res1=list()
    factor_res1$loading_c=factor_res$loading_G
    factor_res1$loading_s=factor_res$loading_F
    factor_res1$fc<-(factor_res$Ghat)[Idx,]
    factor_res1$fs=list()
    for(k in 1:K){

      factor_res1$fs[[k]]<-(factor_res$Fhat)[[k]][Idx,]
    }
    factor_res1$v<-(factor_res$residual)[Idx,]

    factor_res2=list()
    factor_res2$loading_c=factor_res$loading_G
    factor_res2$loading_s=factor_res$loading_F
    factor_res2$fc<-(factor_res$Ghat)[-Idx,]

    factor_res2$fs=list()
    for(k in 1:K){

      factor_res2$fs[[k]]<-(factor_res$Fhat)[[k]][-Idx,]
    }

    # factor_res2$fs<-(factor_res$Fhat)[-Idx,]
    factor_res2$v<-(factor_res$residual)[-Idx,]
  }else{
    factor_res1=factor_res2=NULL
  }


  tmp=1e+100

  res_fac = estimation_fac(XX0, K,rmax=rmax,method=method,weight=weight,type=type,factor_res=NULL,nodewise_step=FALSE)

  fit=estimation(XX0,Y0,K,rmax=rmax,method=method,weight=weight,type=type,lambda.seq=NULL,res_fac =res_fac)
  if(lam_in_iteration){

    fit1=joint_func(XX0,Y0,K,rmax=rmax,method=method,weight=weight, type =type,lambda = lambda,fit=fit,r_update=r_update)

  }else{
     for(i in 1:length(lambda)){

    res=joint_func(XX0,Y0,K,rmax=rmax,method=method,weight=weight, type =type,lambda = lambda[i],fit=fit,r_update=r_update)
    if(tmp>res$BIC_score){
      fit1=res
      tmp=res$BIC_score
    }
  }
  }


  tmp=1e+100

  res_fac = estimation_fac(XX2, K,rmax=rmax,method=method,weight=weight,type=type,factor_res=NULL,nodewise_step=FALSE)

  fit=estimation(XX2,Y2,K,rmax=rmax,method=method,weight=weight,type=type,lambda.seq=NULL,res_fac =res_fac)
  if(lam_in_iteration){

    fit2=joint_func(XX2,Y2,K,rmax=rmax,method=method,weight=weight, type =type,lambda = lambda,fit=fit,r_update=r_update)

  }else{
    for(i in 1:length(lambda)){

      res=joint_func(XX2,Y2,K,rmax=rmax,method=method,weight=weight, type =type,lambda = lambda[i],fit=fit,r_update=r_update)
      if(tmp>res$BIC_score){
        fit2=res
        tmp=res$BIC_score
      }
    }

  }


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

