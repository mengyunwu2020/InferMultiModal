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


args=c(1200,1,1,1,1)
# Number of variables
p <- as.numeric(args[1])


# Settings for generating factors F
# 1 = independent factors
choose_F <- as.numeric(args[2])


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
# 2 = Uniform
choose_B <- as.numeric(args[5])




# Display the chosen settings
cat("Running simulation with settings:\n")
cat(" - p =", p, "\n")
cat(" - choose_V =", choose_V, "(1 = independent, 2 = Band (0.2), 3 = AR (0.2)\n")
cat(" - choose_noise =", choose_noise, "(1 = Gaussian, 2 = Uniform)\n")
cat(" - choose_B =", choose_B, "(1 = Gaussian, 2 = Uniform)\n")




library(GrFA)
library(ncvreg)
library(Matrix)
library(glmnet)
library(SIS)
library(MASS)
#-----------Basic Settings----------------
n=300
set.seed(1)
K=4
r_f=c(2,rep(2,K)) #number of factors
r_all=sum(r_f)
gamma1=0.5*rep(1,r_all)
pm=rep(p/K,K)


if(choose_V==2){
  ## v2
  top_mat <- diag(1,p)
  for(i in 1:(p-1)){
    top_mat[i, i+1] <- .2
    top_mat[i+1, i] <- .2
  }
}else if(choose_V==3){
  top_mat<-matrix(nrow=p,ncol=p) #covariance matrix
  for (i in c(1:p)){
    for (j in c(1:p)){
      top_mat[i,j]=0.2^abs(i-j)
    }
  }
}


beta_0=c(rep(0,p)) #beta under the null hypothesis
#------------------------#beta in the alternative with signal strength w rises from 0.05 to 0.35--------------
beta_1=c(rep(0.05,2),rep(0,p/K-2),rep(0.05,2),rep(0,p/K-2),rep(0.05,2),rep(0,p/K-2),rep(0.05,2),rep(0,p/K-2))
# beta_2=c(rep(0.1,2),rep(0,p/K-2),rep(0.1,2),rep(0,p/K-2),rep(0.1,2),rep(0,p/K-2),rep(0.1,2),rep(0,p/K-2))
# beta_3=c(rep(0.15,2),rep(0,p/K-2),rep(0.15,2),rep(0,p/K-2),rep(0.15,2),rep(0,p/K-2),rep(0.15,2),rep(0,p/K-2))
# beta_4=c(rep(0.2,2),rep(0,p/K-2),rep(0.2,2),rep(0,p/K-2),rep(0.2,2),rep(0,p/K-2),rep(0.2,2),rep(0,p/K-2))
# beta_5=c(rep(0.25,2),rep(0,p/K-2),rep(0.25,2),rep(0,p/K-2),rep(0.25,2),rep(0,p/K-2),rep(0.25,2),rep(0,p/K-2))
# beta_6=c(rep(0.3,2),rep(0,p/K-2),rep(0.3,2),rep(0,p/K-2),rep(0.3,2),rep(0,p/K-2),rep(0.3,2),rep(0,p/K-2))
# beta_7=c(rep(0.35,2),rep(0,p/K-2),rep(0.35,2),rep(0,p/K-2),rep(0.35,2),rep(0,p/K-2),rep(0.35,2),rep(0,p/K-2))


if(choose_B==1){
  B_c<-matrix(rnorm(p*r_f[1]),nrow=p) #generate common factor loading or using gaussian.
}else{
  B_c<-matrix(runif(p*r_f[1],-sqrt(3),sqrt(3)),nrow=p) #generate common factor loading or using gaussian.
}
B_c[,1]<-B_c[,1]/(p^(1/10))

B=B_c
B_s<-list()

if(choose_B==1){
  for(k in 1:K){
    B_s[[k]]<-matrix(rnorm(pm[k]*r_f[k+1],mean = 1),nrow=pm[k])  #generate weak factor loading or using gaussian.
    B_s[[k]][,1]<- B_s[[k]][,1]/(p^(1/10))
  }
}else{
  for(k in 1:K){
    B_s[[k]]<-matrix(runif(pm[k]*r_f[k+1],(-sqrt(3)+1),(sqrt(3)+1)),nrow=pm[k])
    B_s[[k]][,1]<- B_s[[k]][,1]/(p^(1/10))
  }
}
block_diag <- as.matrix(bdiag(B_s[[1]], B_s[[2]]))
block_diag <- as.matrix(bdiag(block_diag, B_s[[3]]))
block_diag <- as.matrix(bdiag(block_diag, B_s[[4]]))

B=cbind(B_c,as.matrix(block_diag))


repli=2000
idx3=idx1<-matrix(NA,nrow = 18,ncol = repli+1)
idx4=idx2=matrix(NA,72,repli+1)
idx5=matrix(NA,36,repli+1)
idx6=matrix(NA,9,repli+1)

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


    ss=1
    XX=list()
    for(k in 1:K){
      XX[[k]]=X[,ss:(ss+pm[k]-1)]
      ss=ss+pm[k]
    }


    #---------------- Main Procedures -------------------
    if (choose_noise==1){
      Y=F%*%gamma1+rnorm(n,0,1)
    }else if (choose_noise==2){
      Y=F%*%gamma1+runif(n,-sqrt(3),sqrt(3))
    }


    #------------Response under the alternative-------------
    # Y2=Y+V%*%beta_1
    # Y3=Y+V%*%beta_2
    # Y4=Y+V%*%beta_3
    # Y5=Y+V%*%beta_4
    # Y6=Y+V%*%beta_5
    # Y7=Y+V%*%beta_6
    # Y8=Y+V%*%beta_7

    #---------ours-----------

    res_fac1 = estimation_fac(XX, K, rmax=10, r0=NULL, r=NULL, method='gap', weight=TRUE, type="BIC3",  factor_res=NULL,nodewise_step=TRUE,boot_rep=1,boot_each=2500,alpha=0.05)


    fit_e1 = estimation(XX,Y, K, rmax=10, r0=NULL, r=NULL, method='gap', weight=TRUE, type="BIC3", lambda.seq=NULL, factor_res=NULL,res_fac =res_fac1)
    fit1=betatest(XX,Y,K,fit_e1,rmax=10,method='gap',weight=TRUE, type = "BIC3",lambda.seq=NULL,alpha=0.05,factor_res=NULL,boot_rep = 1,boot_each = 2500)

    idx1[1,i+1]<-fit1$pvals[1]<0.05
    idx2[1:4,i+1]<-fit1$pvals[2:5]<0.05

    #-----------------Fan---------------
    res_fan=beta_fan(X,n,p,Y)
    idx3[1,i+1]<-res_fan[1]
    #-----------------Li0--------------
    reli=li_estimate_f1(X,XX,K,n,p,Y)
    f_hat_li=reli$hatf_all
    lam=c(0.3,0.01,0.1,0.05,0.005)
    res=rep(0,K)
    ss=1
    for(k in 1:K){
      res[k]=li_betam(Y,X[,-c(ss:(ss+pm[k]-1))],f_hat_li[[k]],lambda2=lam,reli$sigma[1])
      ss=ss+pm[k]
    }
    id_li=res
    idx4[1:4,i+1]=id_li

  }


