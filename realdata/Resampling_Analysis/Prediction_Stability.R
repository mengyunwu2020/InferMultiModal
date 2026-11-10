library(GrFA)
library(ncvreg)
library(Matrix)
library(glmnet)
library(SIS)
library(MASS)
library(here)

devtools::load_all(here::here("InferMultiModal"))

source(here::here("alternative","beta_li.R"), echo=TRUE)

file_data = here::here("Data","RealData","ProcessedData","lung_data.RData")


load(here::here("realdata", "Results", "fit_allsample.Rdata"))

load(here::here("realdata", "Results", "res_allsample.Rdata"))


set.seed(1)
K=4
load(file_data)
Y=lung_df[,1]
X=lung_df[,-1]
pm=c(300,300,300,194)
p=sum(pm)
n=nrow(X)
p=sum(pm)


tmp_mean <- colMeans(X)
X_train <- scale(X, center = tmp_mean, scale = FALSE)

factor_res1=NULL
tmp_mean <- mean(Y)
y_train <- Y-tmp_mean

ss=1
XX_train=list()
for(k in 1:K){
  XX_train[[k]]=as.matrix(X_train[,ss:(ss+pm[k]-1)])
  ss=ss+pm[k]
}
repli=100
test_error_mat=matrix(NA,nrow = repli,ncol=4)
for(rate in c(20)){
  betares=NULL
  gammares=NULL
  for (re in 1:repli){


    tt1=Sys.time()
    set.seed(re)

    sub<-sample(1:nrow(X),round(nrow(X)*rate/(rate+1)))

    X_train<-X[sub,]
    X_test<-X[-sub,]

    tmp_mean <- colMeans(X)
    X_train <- scale(X_train, center = tmp_mean, scale = FALSE)
    X_test <- scale(X_test, center = tmp_mean, scale = FALSE)

    factor_res1=NULL
    factor_res1$Ghat <- res_ori$Fhat[sub, 1, drop = FALSE]
    factor_res1$Fhat=list()
    rtmp=2
    for(k in 1:K){
      fing=res_ori$rr_hat[(k+1)]
      factor_res1$Fhat[[k]]=res_ori$Fhat[sub,(rtmp):(rtmp+fing-1)]
      rtmp=rtmp+fing
    }
    factor_res1$residual=list()
    ptmp=1
    for(k in 1:K){
      factor_res1$residual[[k]]=res_ori$V[sub,ptmp:(ptmp+pm[k]-1)]
      ptmp=ptmp+pm[k]
    }

    B_c=res_ori$loading_c[[1]]
    B_k=res_ori$loading_s[[1]]
    for(k in 2:4){
      B_c=rbind(B_c,res_ori$loading_c[[k]])
      B_k=bdiag(B_k,res_ori$loading_s[[k]])
    }

    factor_res1$loading_G=B_c
    factor_res1$loading_F=as.matrix(B_k)

    y_train<-Y[sub]
    y_test<-Y[-sub]

    tmp_mean <- mean(Y)
    y_train <- y_train-tmp_mean
    y_test <- y_test-tmp_mean



    data=X_train
    n=dim(data)[1]
    n_test=dim(X_test)[1]
    ss=1
    XX_train=list()
    for(k in 1:K){
      XX_train[[k]]=as.matrix(X_train[,ss:(ss+pm[k]-1)])
      ss=ss+pm[k]
    }

    cat("the ", re,"th replicate:","\n")


    ##--ours-----
    set.seed(123)
    res_fac = estimation_fac(XX_train, K, rmax=10, method='gap', weight=TRUE, type='BIC3',  factor_res=factor_res1,nodewise_step=TRUE,boot_rep=1,boot_each=2500,alpha=0.05)

    fit_e111 = estimation(XX_train,y_train, K, rmax=10, method='gap', weight=TRUE, type='BIC3', lambda.seq =seq(1,.01,length.out=100), factor_res=factor_res1,res_fac =res_fac,gamma_MCP =3,tun_ncvreg='CV')

    fit_e111$beta=fit_ori$fit$beta
    fit_e111$gamma=fit_ori$fit$gamma
    fit<-betatest(XX_train,y_train,K,fit_e111,rmax=10,method='gap',weight=TRUE, type ='BIC3',alpha=0.05,factor_res=factor_res1,lambda.seq =seq(1,.01,length.out=100),boot_each = 2500,boot_rep = 1,gamma_MCP =3,seed_sis=123,tune_sis='ebic',tun_ncvreg='CV')

    set.seed(123)
    gres=gammatest(XX_train,y_train,K,ratio1=.79,ratio2=.89,rmax=10,method='gap',weight=TRUE,lambda.seq =seq(1,.01,length.out=100),r_f=c(1,2,3,9,3),alpha=0.05,type="BIC3",factor_res=factor_res1,gamma_MCP =3,seed_sis=1,tune_sis ='ebic',tun_ncvreg='CV')#specific


    betares=rbind(betares,fit$res)
    gammares=rbind(gammares,gres)
    # write.csv(betares,file=here::here("realdata","Results","betares20.csv"))
    # write.csv(gammares,file=here::here("realdata","Results","gammares20.csv"))
    res_fac = estimation_fac(XX_train, K, rmax=10, method='gap', weight=TRUE, type='BIC3',  factor_res=NULL,nodewise_step=FALSE,boot_rep=1,boot_each=2500,alpha=0.05)
    fit_cv = estimation(XX_train,y_train, K, rmax=10, method='gap', weight=TRUE, type='BIC3', lambda.seq =NULL, factor_res=NULL,res_fac =res_fac,gamma_MCP =3,tun_ncvreg='CV')
    beta_hat1<-as.vector(fit_cv$beta)
    gamma_hat1<-as.vector(fit_cv$gamma)
    names(gamma_hat1)=rep(0:K,fit_cv$rr_hat)
    r_f=fit_cv$rr_hat
    if(sum(is.na(fit_cv$loading_c))>0){
      loading_c_mat=fit_cv$loading_c
      loading_s=fit_cv$loading_s
      loading_s_mat=Re(loading_s[[1]])
      for(k in 1:(K-1)){
        loading_s_mat=as.matrix(bdiag(loading_s_mat,Re(loading_s[[k+1]])))
      }}else{
        loading_c=fit_cv$loading_c
        loading_s=fit_cv$loading_s
        loading_c_mat=Re(loading_c[[1]])
        loading_s_mat=Re(loading_s[[1]])
        for(k in 1:(K-1)){
          loading_c_mat=rbind(loading_c_mat,Re(loading_c[[k+1]]))
          loading_s_mat=as.matrix(bdiag(loading_s_mat,Re(loading_s[[k+1]])))
        }

      }



    cvfit <- cv.ncvreg(X_train,y_train, penalty="MCP")
    fit <- cvfit$fit
    cv_res<-fit$beta[,cvfit$min]
    beta_hat2<-as.vector(cv_res[-1])


    cov<-X_train%*%t(X_train)/n
    eigvec<-eigen(cov)$vectors
    eigval<-eigen(cov)$values
    K_est<-which.min(diff(log(eigval[1:10])))
    hatf2<-Re(sqrt(n)*eigvec[,1:K_est])
    hatB2<-t(1/n*t(hatf2)%*%X_train)
    hatB2_fan=hatB2
    hatU2<-X_train-hatf2%*%t(hatB2)
    Y_1<-y_train-1/n*hatf2%*%t(hatf2)%*%y_train
    cvfit3 =cv.ncvreg(hatU2,Y_1,penalty="lasso")
    fit3 <- cvfit3$fit
    cv_res3<-fit3$beta[,cvfit3$min]
    beta_hat3=cv_res3[-1]
    gamma_hat3= 1/n*t(hatf2)%*%y_train


    ##
    hatf_all=NULL
    hatU2=NULL
    rr_hat=rep(0,K)
    hatB2_li=NULL
    for(k in 1:K){
      cov<-XX_train[[k]]%*%t(XX_train[[k]])/n
      eigvec<-eigen(cov)$vectors
      eigval<-eigen(cov)$values
      K_est<-which.min(diff(log(eigval[1:10])))
      rr_hat[k]=K_est
      hatf2<-Re(sqrt(n)*eigvec[,1:K_est])
      hatf_all=cbind(hatf_all,hatf2)
      hatB2<-t(1/n*t(hatf2)%*%XX_train[[k]])
      if(k==1){
        hatB2_li=hatB2
      }else{
        # hatB2    <- Re(hatB2)
        hatB2_li=as.matrix(bdiag(hatB2_li,hatB2))
      }
      tmp=XX_train[[k]]-hatf2%*%t(hatB2)
      hatU2<-cbind(hatU2,tmp)
    }
    penalty_factors1=c(rep(0,dim(hatf_all)[2]),rep(1,p))
    cvfit4 <- cv.ncvreg(cbind(hatf_all,hatU2),y_train,penalty.factor = penalty_factors1, penalty="MCP")
    fit4 <- cvfit4$fit
    cv_res4<-fit4$beta[,cvfit4$min]
    beta_hat4=cv_res4[(dim(hatf_all)[2]+2):(dim(hatf_all)[2]+p+1)]
    gamma_hat4=cv_res4[2:(dim(hatf_all)[2]+1)]
    names(gamma_hat4)=rep(1:K,rr_hat)


    if(dim(loading_s_mat)[1]!=p){
      BtB_inv <- solve(t(loading_c_mat) %*% loading_c_mat)
      projection <- BtB_inv %*% t(loading_c_mat)
      F_hat <- X_test %*% t(projection)

      V_hat <- X_test - F_hat%*%t(loading_c_mat)
      tmp=which(is.na(loading_s_mat),arr.ind = TRUE)[1]

      ss=1
      XX_test=list()
      for(k in 1:K){
        XX_test[[k]]=as.matrix(X_test[,ss:(ss+pm[k]-1)])
        ss=ss+pm[k]
      }
      ss=1
      VV=list()
      for(k in 1:K){
        if(is.element(k,tmp)){
          VV[[k]]=as.matrix(V_hat[,ss:(ss+pm[k]-1)])
        }else{
          loading_s_mat2=loading_s[[k]]
          BtB_inv <- solve(t(loading_s_mat2) %*% loading_s_mat2)
          projection <- BtB_inv %*% t(loading_s_mat2)


          F_hat2<-V_hat[,ss:(ss+pm[k]-1)] %*% t(projection)
          F_hat=cbind(F_hat,F_hat2)
          VV[[k]] <- V_hat[,ss:(ss+pm[k]-1)] - F_hat2%*%t(loading_s_mat2)
        }
        ss=ss+pm[k]
      }

      V_hat= do.call(cbind, VV)
    }else{
      if(sum(is.na(loading_c_mat))>0){
        loading_mat= loading_s_mat
        BtB_inv <- solve(t(loading_mat) %*% loading_mat)       # (B^T B)^{-1}
        projection <- BtB_inv %*% t(loading_mat)     # (B^T B)^{-1} B^T
        F_hat <- X_test %*% t(projection)       #   f_t^T
        V_hat <- X_test - F_hat%*%t(loading_mat)
      }else{
        loading_mat=cbind(loading_c_mat,loading_s_mat)
        BtB_inv <- solve(t(loading_mat) %*% loading_mat)       # (B^T B)^{-1}
        projection <- BtB_inv %*% t(loading_mat)     # (B^T B)^{-1} B^T
        F_hat <- X_test %*% t(projection)       #   f_t^T
        V_hat <- X_test - F_hat%*%t(loading_mat)
      }

    }
    y_predict=V_hat%*%beta_hat1+F_hat%*%gamma_hat1
    test_error1=mean((y_test-y_predict)^2)


    y_predict=X_test%*%beta_hat2
    test_error2=mean((y_test-y_predict)^2)



    BtB_inv <- solve(t(hatB2_fan) %*% hatB2_fan)       # (B^T B)^{-1}
    projection <- BtB_inv %*% t(hatB2_fan)     # (B^T B)^{-1} B^T
    F_hat <- X_test %*% t(projection)       #   f_t^T
    V_hat <- X_test - F_hat%*%t(hatB2_fan)
    y_predict=V_hat%*%beta_hat3+F_hat%*%gamma_hat3
    test_error3=mean((y_test-y_predict)^2)



    BtB_inv <- solve(t(hatB2_li) %*% hatB2_li)
    projection <- BtB_inv %*% t(hatB2_li)
    F_hat <- X_test %*% t(projection)
    V_hat <- X_test - F_hat%*%t(hatB2_li)
    y_predict=V_hat%*%beta_hat4+F_hat%*%gamma_hat4
    test_error4=mean((y_test-y_predict)^2)

    test_error_mat[re,]=c(test_error1,test_error2,test_error3,test_error4)
    # write.csv(test_error_mat,file=here::here("realdata","Results","rate20error.csv"))

    tt2=Sys.time()
    print(tt2-tt1)
  }

}

