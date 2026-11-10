beta_fan<-function(X,n,p,Y){


  X<-X-rowMeans(X) #normalize
  #-----------Use covariance matrix to estimate factors-----------
  cov<-X%*%t(X)/n #covariance matrix
  eigvec<-eigen(cov)$vectors
  eigval<-eigen(cov)$values
  KK<-which.min(diff(log(eigval[1:10]))) #select number of factors
  hatf<-Re(sqrt(n)*eigvec[,1:KK]) #estimated factors
  hatB<-t(1/n*t(hatf)%*%X) #estimate factor loading
  hatU<-X-hatf%*%t(hatB) #estimate idiosyncratic component
  #-------------------Use Nodewise regression to estimate Theta, the invese matrix of the population covariance------
  C<-as.matrix(diag(rep(1,p)))
  T<-c()
  for( j in 1:ncol(hatU)){
    fit_u = glmnet(hatU[,-j], hatU[,j], intercept=FALSE,
                   lambda=cv.glmnet(hatU[,-j], hatU[,j],intercept=FALSE)$lambda.1se)
    beta<-as.vector(fit_u$beta)
    C[j,-j]<--beta
    T<-c(T,1/n*sum((hatU[,j]-hatU[,-j]%*%beta)^2)+fit_u$lambda/2*sum(abs(beta)))

  }
  T1<-diag(1/T)
  Theta<-T1%*%C #estimated Theta
  #------------------Conduct bootstrap to compute the critical values----------
  quant<-c() #threshold
  for (m in 1:5){
    boostrap<-c()
    for(j in 1:200){
      boostrap<-c(boostrap,max(abs(Theta%*%t(hatU)%*%rnorm(n,mean=0,sd=1))))
    }
    quant<-c(quant,quantile(boostrap,0.95))
  }
  c_alpha<-1/sqrt(n)*mean(quant)


  #-------------Estimate the variances----------------
  sd_alter1<-rcv_fan(Y,X,n)
  #-----------Factor adjusted Y in those five different regimes --------
  Y_new1<-Y-1/n*hatf%*%t(hatf)%*%Y #under the null
  #-----------Model fitting in those five different regimes-----------------
  fit1 = glmnet(hatU, Y_new1, intercept=FALSE,
                lambda=cv.glmnet(hatU,Y_new1,intercept=FALSE)$lambda.1se) #under the null
  #----------------Fitted values------------------------------------------
  lambda_beta1<-fit1$lambda #fitted tuning parameters
  beta_hat1<-as.vector(fit1$beta)
  d_beta1<-beta_hat1+1/n*Theta%*%t(hatU)%*%(Y-hatU%*%beta_hat1)

  idx1<-c()
  #-------------------Compute the size of the test or power of the test----------------------
  idx1<-c(idx1,sum(sqrt(n)*max(abs(d_beta1))/sqrt(sd_alter1[1])>c_alpha))

  res=c(idx1)
  return(res)


}


#-------------Helper function: rcv_fan: compute the variance of Y given X-----------
#input Y: response
#X: covariate
#n number of observations
rcv_fan<-function(Y,X,n){
  Idx = 1:(n/2)
  X0 = X[Idx,]; X2 = X[-Idx,]
  Sigma0 = tcrossprod(X0)*2/n; Sigma2 = tcrossprod(X2)*2/n
  eigen0 = eigen(Sigma0); eigen2 = eigen(Sigma2)
  eigvec0 = eigen0$vectors; eigvalue0 = eigen0$values
  eigvec2 = eigen2$vectors; eigvalue2 = eigen2$values
  K0 = which.min(diff(log(eigvalue0[1:100])))
  K2 = which.min(diff(log(eigvalue2[1:100])))
  F.hat0 = eigvec0[,1:K0]*sqrt(n/2); F.hat2 = eigvec2[,1:K2]*sqrt(n/2)
  B.hat0.T = t(F.hat0)%*%X0*2/n; B.hat2.T = t(F.hat2)%*%X2*2/n
  U.hat0 = X0 - F.hat0%*%B.hat0.T; U.hat2 = X2 - F.hat2%*%B.hat2.T
  Y0 = Y[Idx]; Y2 = Y[-Idx]
  tmp0 = tcrossprod(F.hat0); tmp2 = tcrossprod(F.hat2)
  Y0.new = Y0 - tmp0%*%Y0*2/n
  Y2.new = Y2 - tmp2%*%Y2*2/n
  SIS0 = SIS(U.hat0,Y0.new)
  SIS2 = SIS(U.hat2,Y2.new)
  S.hat0.ISIS = SIS0$ix; S.hat0.SIS = SIS0$ix0
  S.hat2.ISIS = SIS2$ix; S.hat2.SIS = SIS2$ix0
  X0=cbind(F.hat0,U.hat0)
  X2=cbind(F.hat2,U.hat2)
  X0.ISIS = X0[,c(1:K0,S.hat2.ISIS+K0)]
  X0.SIS = X0[,c(1:K0,K0+S.hat2.SIS)]
  X2.ISIS = X2[,c(1:K2,K2+S.hat0.ISIS)]; X2.SIS = X2[,c(1:K2,K2+S.hat0.SIS)]
  lm0.ISIS = lm(Y0~X0.ISIS-1)
  lm0.SIS = lm(Y0~X0.SIS-1)
  lm2.ISIS = lm(Y2~X2.ISIS-1)
  lm2.SIS = lm(Y2~X2.SIS-1)
  Q0.ISIS.H1 = sum((resid(lm0.ISIS))^2)
  sigma.hat0.ISIS = Q0.ISIS.H1/(n/2 - K0 - length(S.hat2.ISIS))
  Q0.SIS.H1 = sum((resid(lm0.SIS))^2)
  sigma.hat0.SIS = Q0.SIS.H1/(n/2 - K0 - length(S.hat2.SIS))
  Q2.ISIS.H1 = sum((resid(lm2.ISIS))^2)
  sigma.hat2.ISIS = Q2.ISIS.H1/(n/2 - K2 - length(S.hat0.ISIS))
  Q2.SIS.H1 = sum((resid(lm2.SIS))^2)
  sigma.hat2.SIS = Q2.SIS.H1/(n/2 - K2 - length(S.hat0.SIS))
  sigma.hat.ISIS = mean(c(sigma.hat0.ISIS, sigma.hat2.ISIS))
  sigma.hat.SIS = mean(c(sigma.hat0.SIS, sigma.hat2.SIS))
  c(sigma.hat.ISIS,sigma.hat.SIS)} #output estimated variance using (ISIS method, SIS method), respectively.
#---------------- Main Procedures -------------------

beta_fan9<-function(X,n,p,Y,Y2,Y3,Y4,Y5,Y6,Y7,Y8,Y9){

  X<-X-rowMeans(X) #normalize
  #-----------Use covariance matrix to estimate factors-----------
  cov<-X%*%t(X)/n #covariance matrix
  eigvec<-eigen(cov)$vectors
  eigval<-eigen(cov)$values
  KK<-which.min(diff(log(eigval[1:10]))) #select number of factors
  hatf<-sqrt(n)*eigvec[,1:KK] #estimated factors
  hatB<-t(1/n*t(hatf)%*%X) #estimate factor loading
  hatU<-X-hatf%*%t(hatB) #estimate idiosyncratic component
  #-------------------Use Nodewise regression to estimate Theta, the invese matrix of the population covariance------
  C<-as.matrix(diag(rep(1,p)))
  T<-c()
  for( j in 1:ncol(hatU)){
    fit_u = glmnet(hatU[,-j], hatU[,j], intercept=FALSE,
                   lambda=cv.glmnet(hatU[,-j], hatU[,j],intercept=FALSE)$lambda.1se)
    beta<-as.vector(fit_u$beta)
    C[j,-j]<--beta
    T<-c(T,1/n*sum((hatU[,j]-hatU[,-j]%*%beta)^2)+fit_u$lambda/2*sum(abs(beta)))

  }
  T1<-diag(1/T)
  Theta<-T1%*%C #estimated Theta
  #------------------Conduct bootstrap to compute the critical values----------
  quant<-c() #threshold
  for (m in 1:5){
    boostrap<-c()
    for(j in 1:200){
      boostrap<-c(boostrap,max(abs(Theta%*%t(hatU)%*%rnorm(n,mean=0,sd=1))))
    }
    quant<-c(quant,quantile(boostrap,0.95))
  }
  c_alpha<-1/sqrt(n)*mean(quant)

  sd_alter1<-rcv_fan(Y,X,n) #variances under the alternative
  sd_alter2<-rcv_fan(Y2,X,n)

  sd_alter3<-rcv_fan(Y3,X,n)
  sd_alter4<-rcv_fan(Y4,X,n) #variances under the alternative
  sd_alter5<-rcv_fan(Y5,X,n)
  sd_alter6<-rcv_fan(Y6,X,n) #variances under the alternative
  sd_alter7<-rcv_fan(Y7,X,n)
  sd_alter8<-rcv_fan(Y8,X,n) #variances under the alternative
  sd_alter9<-rcv_fan(Y9,X,n)

  Y_new1<-Y-1/n*hatf%*%t(hatf)%*%Y
  Y_new2<-Y2-1/n*hatf%*%t(hatf)%*%Y2
  Y_new3<-Y3-1/n*hatf%*%t(hatf)%*%Y3

  Y_new4<-Y4-1/n*hatf%*%t(hatf)%*%Y4
  Y_new5<-Y5-1/n*hatf%*%t(hatf)%*%Y5

  Y_new6<-Y6-1/n*hatf%*%t(hatf)%*%Y6
  Y_new7<-Y7-1/n*hatf%*%t(hatf)%*%Y7

  Y_new8<-Y8-1/n*hatf%*%t(hatf)%*%Y8
  Y_new9<-Y9-1/n*hatf%*%t(hatf)%*%Y9

  fit1 = glmnet(hatU, Y_new1, intercept=FALSE,
                lambda=cv.glmnet(hatU,Y_new1,intercept=FALSE)$lambda.1se)

  fit2 = glmnet(hatU, Y_new2, intercept=FALSE,
                lambda=cv.glmnet(hatU,Y_new2,intercept=FALSE)$lambda.1se)
  fit3 = glmnet(hatU, Y_new3, intercept=FALSE,lambda=cv.glmnet(hatU,Y_new2,intercept=FALSE)$lambda.1se)

  fit4 = glmnet(hatU, Y_new4, intercept=FALSE,
                lambda=cv.glmnet(hatU,Y_new4,intercept=FALSE)$lambda.1se)
  fit5 = glmnet(hatU, Y_new5, intercept=FALSE,lambda=cv.glmnet(hatU,Y_new4,intercept=FALSE)$lambda.1se)

  fit6 = glmnet(hatU, Y_new6, intercept=FALSE,
                lambda=cv.glmnet(hatU,Y_new6,intercept=FALSE)$lambda.1se)
  fit7 = glmnet(hatU, Y_new7, intercept=FALSE,lambda=cv.glmnet(hatU,Y_new6,intercept=FALSE)$lambda.1se)

  fit8 = glmnet(hatU, Y_new8, intercept=FALSE,
                lambda=cv.glmnet(hatU,Y_new8,intercept=FALSE)$lambda.1se)
  fit9 = glmnet(hatU, Y_new9, intercept=FALSE,lambda=cv.glmnet(hatU,Y_new8,intercept=FALSE)$lambda.1se)


  lambda_beta1<-fit1$lambda
  lambda_beta2<-fit2$lambda
  lambda_beta3<-fit3$lambda

  lambda_beta4<-fit4$lambda
  lambda_beta5<-fit5$lambda


  lambda_beta6<-fit6$lambda
  lambda_beta7<-fit7$lambda
  lambda_beta9<-fit9$lambda
  lambda_beta8<-fit8$lambda

  beta_hat1<-as.vector(fit1$beta)
  d_beta1<-beta_hat1+1/n*Theta%*%t(hatU)%*%(Y-hatU%*%beta_hat1) #



  beta_hat2<-as.vector(fit2$beta)
  beta_hat3<-as.vector(fit3$beta)
  d_beta2<-beta_hat2+1/n*Theta%*%t(hatU)%*%(Y2-hatU%*%beta_hat2) #debiased beta under the alternative
  d_beta3<-beta_hat3+1/n*Theta%*%t(hatU)%*%(Y3-hatU%*%beta_hat3)


  beta_hat4<-as.vector(fit4$beta)
  beta_hat5<-as.vector(fit5$beta)
  d_beta4<-beta_hat4+1/n*Theta%*%t(hatU)%*%(Y4-hatU%*%beta_hat4) #debiased beta under the alternative
  d_beta5<-beta_hat5+1/n*Theta%*%t(hatU)%*%(Y5-hatU%*%beta_hat5)


  beta_hat6<-as.vector(fit6$beta)
  beta_hat7<-as.vector(fit7$beta)
  d_beta6<-beta_hat6+1/n*Theta%*%t(hatU)%*%(Y6-hatU%*%beta_hat6) #debiased beta under the alternative
  d_beta7<-beta_hat7+1/n*Theta%*%t(hatU)%*%(Y7-hatU%*%beta_hat7)


  beta_hat8<-as.vector(fit8$beta)
  beta_hat9<-as.vector(fit9$beta)
  d_beta8<-beta_hat8+1/n*Theta%*%t(hatU)%*%(Y8-hatU%*%beta_hat8) #debiased beta under the alternative
  d_beta9<-beta_hat9+1/n*Theta%*%t(hatU)%*%(Y9-hatU%*%beta_hat9)




  idx1<-c()
  idx2<-c()
  idx3<-c()
  idx4<-c()
  idx5<-c()
  idx6<-c()
  idx7<-c()
  idx8<-c()
  idx9<-c()


  idx1<-c(idx1,sum(sqrt(n)*max(abs(d_beta1))/sqrt(sd_alter1[1])>c_alpha))
  idx4<-c(idx4,sum(sqrt(n)*max(abs(d_beta4))/sqrt(sd_alter4[1])>c_alpha)) #size of the test under the alternative
  idx5<-c(idx5,sum(sqrt(n)*max(abs(d_beta5))/sqrt(sd_alter5[1])>c_alpha))
  idx2<-c(idx2,sum(sqrt(n)*max(abs(d_beta2))/sqrt(sd_alter2[1])>c_alpha)) #size of the test under the alternative
  idx3<-c(idx3,sum(sqrt(n)*max(abs(d_beta3))/sqrt(sd_alter3[1])>c_alpha))
  idx6<-c(idx6,sum(sqrt(n)*max(abs(d_beta6))/sqrt(sd_alter6[1])>c_alpha)) #size of the test under the alternative
  idx7<-c(idx7,sum(sqrt(n)*max(abs(d_beta7))/sqrt(sd_alter7[1])>c_alpha))
  idx8<-c(idx8,sum(sqrt(n)*max(abs(d_beta8))/sqrt(sd_alter8[1])>c_alpha)) #size of the test under the alternative
  idx9<-c(idx9,sum(sqrt(n)*max(abs(d_beta9))/sqrt(sd_alter9[1])>c_alpha))

  res=c(idx1,idx2,idx3,idx4,idx5,idx6,idx7,idx8,idx9)
  return(res)

}
