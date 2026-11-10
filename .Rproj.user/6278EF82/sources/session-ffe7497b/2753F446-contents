library(GrFA)
library(ncvreg)
library(Matrix)
library(glmnet)
library(SIS)
library(MASS)

library(here)
source(here::here("alternative", "beta_fan.R"))
devtools::load_all(here::here("InferMultiModal"))
source(here::here("alternative","joint_func.R"), echo=TRUE)
source(here::here("alternative","beta_li.R"))

rate_F <- seq(0,55,by=5)
args=c(1200,1,1,1,1)
# Number of variables
p <- as.numeric(args[1])
library(metap)

# Settings for generating factors F
# 1 = independent factors
choose_F=as.numeric(args[2])

# Settings for generating error structure V
# 1 = independent errors
# 2 = banded errors -- Band (0.2)
# 3 = moderate mixing errors -- AR (0.2)
choose_V <- as.numeric(args[3])

# Settings for noise distribution
# 1 = Gaussian noise -- N(0,1)
# 2 = Uniform noise -- Uniform [-sqrt(3),sqrt(3)]
choose_noise <- as.numeric(args[4])



# Settings for loading distribution
# 1 = Gaussian N(0,1)
# 2 = Uniform [-sqrt(3),sqrt(3)]
choose_B <- as.numeric(args[5])



# Display the chosen settings
cat("Running simulation with settings:\n")
cat(" - p =", p, "\n")
cat(" - choose_V =", choose_V, "(1 = independent, 2 = Band (0.2), 3 = AR (0.2)\n")
cat(" - choose_noise =", choose_noise, "(1 = Gaussian, 2 = Uniform)\n")
cat(" - choose_B =", choose_B, "(1 = Gaussian, 2 = Uniform)\n")


n=300 #number of observations

set.seed(1)
K=4
r_f=c(2,rep(2,K)) #number of factors
r_all=sum(r_f)
pm=rep(p/K,K)

#-------------Compute the critical value-------------
c_alpha=rep(0,5)
for(ii in 1:5){

  c_alpha[ii]<-qchisq(0.95,df=r_f[ii])
}

if(choose_V==2){
  ## v2
  top_mat <- diag(1,p)
  for(i in 1:(p-1)){
    top_mat[i, i+1] <- .2
    top_mat[i+1, i] <- .2
  }
}else if(choose_V==3){
  ## v3
  top_mat<-matrix(nrow=p,ncol=p) #covariance matrix
  for (i in c(1:p)){
    for (j in c(1:p)){
      top_mat[i,j]=0.2^abs(i-j)
    }
  }
}



beta=c(rep(1,2),rep(0,pm[1]-2),rep(1,2),rep(0,pm[2]-2),rep(1,2),rep(0,pm[3]-2),rep(1,2),rep(0,pm[4]-2)) #generate true beta



if(choose_B==1){
  B_c<-matrix(rnorm(p*r_f[1],sd=1),nrow=p) #generate common factor loading or using gaussian.
}else{
  B_c<-matrix(runif(p*r_f[1],-sqrt(3),sqrt(3)),nrow=p) #generate common factor loading or using gaussian.
}

B=B_c
B_s<-list()

if(choose_B==1){
  for(k in 1:K){
    B_s[[k]]<-matrix(rnorm(pm[k]*r_f[k+1],sd=1),nrow=pm[k]) #generate common factor loading or using gaussian.
  }
}else{
  for(k in 1:K){
    B_s[[k]]<-matrix(runif(pm[k]*r_f[k+1],-sqrt(3),sqrt(3)),nrow=pm[k])
  }
}
block_diag <- as.matrix(bdiag(B_s[[1]], B_s[[2]]))
block_diag <- as.matrix(bdiag(block_diag, B_s[[3]]))
block_diag <- as.matrix(bdiag(block_diag, B_s[[4]]))

B=cbind(B_c,as.matrix(block_diag))






#-----------Record the size or power of test when $m$ takes 6 values------------#
mat1<-NULL
idx<-NULL
idx2<-NULL
for (m in rate_F){
  gamma1=  c(rep(m/10,r_all))
  isis1<-NULL
  isis2=NULL
  isis3<-NULL
  isis4=NULL
  isis5=NULL
  sigma_isis<-c()
  sigma_sis<-c()


  res=list()
  res1=res2=res3=res4=res5=res6=res7=res_li=NULL
  res1_ours=res_ours=NULL
  for (i in 1:2000){

    set.seed(i)
    if (choose_F==1){
      F<-matrix(rnorm(n*r_all),nrow=n) #generate factor
    }



    if(choose_V==1){
      V<-matrix(rnorm(n*p),nrow=n)
    }else{
      V<-mvrnorm(n, rep(0,p),top_mat) #generate idiosyncratic component
    }

    X=F%*%t(B)+V #design matrix
    if (choose_noise==1){ #generate noise distribution 1:Gaussian; 2: Uniform
      Y2=F%*%gamma1+V%*%beta+rnorm(n,mean=0,sd=1)

      sigma=1
    }else if (choose_noise==2){
      Y2=F%*%gamma1+V%*%beta+runif(n,-sqrt(3),sqrt(3))

      sigma=1
    }

    factor_res1=NULL
    factor_res1$loading_G=B_c
    factor_res1$loading_F=as.matrix(block_diag)
    factor_res1$Ghat=F[,1:2]
    factor_res1$Fhat=list()
    rtmp=3
    for(k in 1:K){

      factor_res1$Fhat[[k]]=F[,(rtmp):(rtmp+1)]
      rtmp=rtmp+2
    }
    factor_res1$residual=list()
    ptmp=1
    for(k in 1:K){

      factor_res1$residual[[k]]=V[,ptmp:(ptmp+pm[k]-1)]
      ptmp=ptmp+pm[k]
    }





    ss=1
    XX=list()
    for(k in 1:K){
      XX[[k]]=X[,ss:(ss+pm[k]-1)]
      ss=ss+pm[k]
    }

    #---------using ours-----------
    gres=gammatest(XX,Y2,K,ratio1=0.9,ratio2=0.6,rmax=10,method='gap',weight=TRUE,lambda.seq=NULL,r_f=r_f,alpha=0.05,type="BIC3",factor_res=NULL)#specific
    gres2=gammatest(XX,Y2,K,ratio1=0.9,ratio2=0.52,rmax=10,method='gap',weight=TRUE,lambda.seq=NULL,r_f=r_f,alpha=0.05,type="BIC3",factor_res=NULL)#common
    res_ours=rbind(res_ours,gres[1:(K+1)]<0.05)
    res1_ours=rbind(res1_ours,gres2[1:(K+1)]<0.05)
    #-----------------Li--------------
    reli=li_estimate_f1(X,XX,K,n,p,Y2)
    f_hat_li=reli$hatf_all
    lam=c(0.3,0.01,0.1,0.05,0.005)
    resfu=rep(0,K)
    ss=1
    for(k in 1:K){
      resfu[k]=li_betam(Y2,X[,-c(ss:(ss+pm[k]-1))],f_hat_li[[k]],lambda2=lam,reli$sigma[1])
      ss=ss+pm[k]
    }
    id_li=resfu


    res_li=rbind(res_li,id_li)

    #---------oracle-based-----------
    gres=gammatest(XX,Y2,K,ratio1=0.9,ratio2=0.6,rmax=10,r0=2,r=rep(2,K),method='gap',weight=TRUE,lambda.seq=NULL,r_f=r_f,alpha=0.05,type="BIC3",factor_res=factor_res1)#specific
    gres2=gammatest(XX,Y2,K,ratio1=0.9,ratio2=0.52,rmax=10,r0=2,r=rep(2,K),method='gap',weight=TRUE,lambda.seq=NULL,r_f=r_f,alpha=0.05,type="BIC3",factor_res=factor_res1)#common
    res=rbind(res,gres[1:(K+1)]<0.05)
    res1=rbind(res1,gres2[1:(K+1)]<0.05)
    #---------semi-oracle-based-----------
    gres=gammatest(XX,Y2,K,ratio1=0.9,ratio2=0.6,rmax=10,r0=2,r=rep(2,K),method='gap',weight=TRUE,lambda.seq=NULL,r_f=r_f,alpha=0.05,type="BIC3",factor_res=NULL)#specific
    gres2=gammatest(XX,Y2,K,ratio1=0.9,ratio2=0.52,rmax=10,r0=2,r=rep(2,K),method='gap',weight=TRUE,lambda.seq=NULL,r_f=r_f,alpha=0.05,type="BIC3",factor_res=NULL)#common
    res2=rbind(res2,gres[1:(K+1)]<0.05)
    res3=rbind(res3,gres2[1:(K+1)]<0.05)

    #---------joint-----------
    gres= gtest_joint(XX,Y2,K,ratio1=0.9,ratio2=0.6,rmax=10,r0=NULL,r=NULL,method='gap',weight=TRUE,type="BIC3",lambda.seq=NULL,r_f=r_f,alpha=0.05,factor_res=NULL,lambda=seq(.1,0.01,length.out=10),r_update = TRUE,lam_in_iteration=TRUE)
    gres2= gtest_joint(XX,Y2,K,ratio1=0.9,ratio2=0.52,rmax=10,r0=NULL,r=NULL,method='gap',weight=TRUE,type="BIC3",lambda.seq=NULL,r_f=r_f,alpha=0.05,factor_res=NULL,lambda=seq(.1,0.01,length.out=10),r_update = TRUE,lam_in_iteration=TRUE)
    res4=rbind(res4,gres<0.05)
    res5=rbind(res5,gres2<0.05)


  }
}




