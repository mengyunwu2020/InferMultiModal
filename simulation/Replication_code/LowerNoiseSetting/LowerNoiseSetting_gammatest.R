
library("SIS")
library("MASS")
library(GrFA)
library(ncvreg)
library(Matrix)
library(glmnet)
library(ACAT)
library(metap)
library(harmonicmeanp)
library(here)
source(here::here("alternative", "beta_fan.R"))
devtools::load_all(here::here("InferMultiModal"))
source(here::here("alternative","joint_func.R"), echo=TRUE)
source(here::here("alternative","beta_li.R"))



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
# 1 = Gaussian noise
choose_noise <- as.numeric(args[4])



# Settings for loading distribution
# 1 = Gaussian
choose_B <- as.numeric(args[5])



# Display the chosen settings
cat("Running simulation with settings:\n")
cat(" - p =", p, "\n")
cat(" - choose_V =", choose_V, "(1 = independent, 2 = Band (0.2), 3 = AR (0.2))\n")
cat(" - choose_noise =", choose_noise, "(1 = Gaussian, 2 = Uniform)\n")
cat(" - choose_B =", choose_B, "(1 = Gaussian, 2 = Uniform)\n")


n=300

set.seed(1)
K=4
r_f=c(2,rep(2,K)) #number of factors
r_all=sum(r_f)
pm=rep(p/K,K)

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




beta=c(rep(1,2),rep(0,pm[1]-2),rep(1,2),rep(0,pm[2]-2),rep(1,2),rep(0,pm[3]-2),rep(1,2),rep(0,pm[4]-2))


if(choose_B==1){
  B_c<-matrix(rnorm(p*r_f[1],sd=1),nrow=p) #generate common factor loading or using gaussian.
}else{
  B_c<-matrix(runif(p*r_f[1],-sqrt(3),sqrt(3)),nrow=p) #generate common factor loading or using gaussian.
}

B=B_c
B_s<-list()


if(choose_B==1){
  for(k in 1:K){
    B_s[[k]]<-matrix(rnorm(pm[k]*r_f[k+1],mean = 0),nrow=pm[k]) #generate weak factor loading or using gaussian.
  }
}else{
  for(k in 1:K){
    B_s[[k]]<-matrix(runif(pm[k]*r_f[k+1],(-sqrt(3)),(sqrt(3))),nrow=pm[k])
  }
}

block_diag <- as.matrix(bdiag(B_s[[1]], B_s[[2]]))
block_diag <- as.matrix(bdiag(block_diag, B_s[[3]]))
block_diag <- as.matrix(bdiag(block_diag, B_s[[4]]))

B=cbind(B_c,as.matrix(block_diag))






ratio1=0.9

mat1<-NULL
idx<-NULL
idx2<-NULL
idx3<-NULL
for (m in seq(0,50,by=5)){
   gamma1=  c(rep(m/10,r_all))
  isis1<-NULL
  isis2=NULL
  isis3<-NULL
  isis4=NULL
  isis5=NULL
  sigma_isis<-c()
  sigma_sis<-c()


  res=NULL
  res1=res2=res3=res4=res5=res6=res7=res8=res9=res10= rrhat1=NULL

  for (i in 1:2000){
    t1=Sys.time()
    set.seed(i)
    if (choose_F==1){
      F<-matrix(rnorm(n*r_all),nrow=n)
    }



    if(choose_V==1){
      V<-matrix(rnorm(n*p),nrow=n)
    }else{
      V<-mvrnorm(n, rep(0,p),top_mat) #generate idiosyncratic component
    }

    X=F%*%t(B)+V #design matrix
    Y2=F%*%gamma1+V%*%beta+rnorm(n,0,.6)

    factor_res1=NULL
    factor_res1$loading_G=B_c
    factor_res1$loading_F=as.matrix(block_diag)
    factor_res1$Ghat=F[,1:r_f[1]]
    factor_res1$Fhat=list()
    rtmp=r_f[1]+1
    for(k in 1:K){

      factor_res1$Fhat[[k]]=F[,(rtmp):(rtmp+r_f[k+1]-1)]
      rtmp=rtmp+r_f[k+1]
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


    #---------joint betatest-----------
    gres= gtest_joint(XX,Y2,K,ratio1=0.85,ratio2=0.5,rmax=10,r0=NULL,r=NULL,method='gap',weight=TRUE,type="BIC3",lambda.seq=NULL,r_f=r_f,alpha=0.05,factor_res=NULL,lambda=seq(.1,0.01,length.out=10),r_update = TRUE,lam_in_iteration=TRUE)
    gres2= gtest_joint(XX,Y2,K,ratio1=0.85,ratio2=0.51,rmax=10,r0=NULL,r=NULL,method='gap',weight=TRUE,type="BIC3",lambda.seq=NULL,r_f=r_f,alpha=0.05,factor_res=NULL,lambda=seq(.1,0.01,length.out=10),r_update = TRUE,lam_in_iteration=TRUE)
    res4=rbind(res4,gres<0.05)
    res5=rbind(res5,gres2<0.05)
    #---------using ours-----------
    gres=gammatest(XX,Y2,K,ratio1=0.85,ratio2=0.5,rmax=10,method='gap',weight=TRUE,lambda.seq=NULL,r_f=r_f,alpha=0.05,type="BIC3",factor_res=NULL)#specific IC3
    gres2=gammatest(XX,Y2,K,ratio1=0.85,ratio2=0.51,rmax=10,method='gap',weight=TRUE,lambda.seq=NULL,r_f=r_f,alpha=0.05,type="BIC3",factor_res=NULL)#common IC3
    res9=rbind(res9,gres[1:(K+1)]<0.05)
    res10=rbind(res10,gres2[1:(K+1)]<0.05)
    rrhat1=rbind(rrhat1,gres[-c(1:(K+1))])
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


    res8=rbind(res8,id_li)

  }

}


