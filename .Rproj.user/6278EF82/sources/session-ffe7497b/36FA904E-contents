#' @export
rcv<-function(Y,X,n,pm,K,rmax=10,r0=NULL,r=NULL,method='gap',weight=TRUE, type = "BIC3",lambda.seq = NULL,factor_res=NULL,gamma_MCP=3,seed_sis,tune_sis,tun_ncvreg){
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
    factor_res1$loading_G=factor_res$loading_G
    factor_res1$loading_F=factor_res$loading_F
    factor_res1$Ghat<-(factor_res$Ghat)[Idx,,drop=FALSE]


    factor_res1$Fhat=list()
    factor_res1$residual=list()
    for(k in 1:K){

      factor_res1$Fhat[[k]]<-factor_res$Fhat[[k]][Idx,]

      factor_res1$residual[[k]]<-factor_res$residual[[k]][Idx,]
    }


    factor_res2=list()
    factor_res2$loading_G=factor_res$loading_G
    factor_res2$loading_F=factor_res$loading_F
    factor_res2$Ghat<-(factor_res$Ghat)[-Idx,,drop=FALSE]

    factor_res2$Fhat=list()
    factor_res2$residual=list()
    for(k in 1:K){

      factor_res2$Fhat[[k]]<-factor_res$Fhat[[k]][-Idx,]

      factor_res2$residual[[k]]<-factor_res$residual[[k]][-Idx,]
    }
  }else{
    factor_res1=factor_res2=NULL
  }


  res_fac = estimation_fac(XX0, K, rmax=rmax, r0=r0,r=r,method=method,weight=weight,type = type, factor_res=factor_res1,nodewise_step=FALSE)

  fit1<-estimation(XX0,Y0,K,rmax=rmax,r0=r0,r=r,method=method,weight=weight,lambda.seq = lambda.seq, type = type,factor_res=factor_res1,res_fac=res_fac,gamma_MCP=gamma_MCP,tun_ncvreg=tun_ncvreg)

  res_fac = estimation_fac(XX2, K, rmax=rmax, r0=r0,r=r,method=method,weight=weight,type = type, factor_res=factor_res2,nodewise_step=FALSE)

  fit2<-estimation(XX2,Y2,K,rmax=rmax,r0=r0,r=r,method=method,weight=weight,lambda.seq = lambda.seq, type = type,factor_res=factor_res2,res_fac=res_fac,gamma_MCP=gamma_MCP,tun_ncvreg=tun_ncvreg)

  r_all<-dim(fit1$FV)[2]-p;  r_all2<-dim(fit2$FV)[2]-p
  V.hat0 = fit1$FV[,(r_all+1):(r_all+p)]; V.hat2 = fit2$FV[,(r_all2+1):(r_all2+p)]
  F.hat0 = fit1$FV[,1:(r_all)]; F.hat2 = fit2$FV[,1:(r_all2)]
  Y0.new = Y0 -  F.hat0%*%fit1$gamma; Y2.new = Y2 -  F.hat2%*%fit2$gamma

  SIS0 = SIS(V.hat0,Y0.new,tune =tune_sis,seed=seed_sis); SIS2 = SIS(V.hat2,Y2.new,tune =tune_sis,seed=seed_sis)
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
