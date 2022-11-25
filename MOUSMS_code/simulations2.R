###########################Supported functions for the 5 methods we compared:
#1. An intercept will be included in X as.default.
#2. Parameters are denoted as : psi = c(beta, sigmasq)
#3. models is a list with one entry for each model. Each entry is an integer vector
#that specifies the columns of matrix data to be used as a regressor in that model.

library(rlist)
library(MASS)
library(parallel)
source("Nested_glm.R")#nested MCS method
source("NISfuns.R")#Greedy SIS for GAM
source("mcb_cv.R")
library(SIS)
library(ggplot2)

ar1_cor <- function(n, rho) {
  exponent <- abs(matrix(1:n - 1, nrow = n, ncol = n, byrow = TRUE) -
                    (1:n - 1))
  return(rho^exponent)
}
Generate.X <- function(case,n,p, rho, fixed.b=F){
  non.zeros = 3
  b=rep(1, non.zeros)
  if(fixed.b==F){
  z = abs(rnorm(non.zeros)); signs = rbinom(non.zeros,1,0.4)
  zero.entries = which(signs==0); signs[zero.entries] = -1
  f1 = 5*log(n)/n^(1/2); f2 = z; b = (f1 + f2)*signs
  }
  true.beta = rep(0,p)
  true.beta[1:non.zeros] = b
  if(case==0){
    non.zeros = 2
    true.beta = c(1,1,0)
    x1 = rbinom(n, 1, 0.5)
    x2 = rnorm(n,0,1)
    x= cbind(x1 = x1, x2=x2, x2sq = x2^2)
  }
  if(case == 1){#independent features
    x = matrix(rnorm(n*p, mean=0, sd=1), n, p)
  }
  if(case==2){#exponential decay correlation
    Sigma=ar1_cor(p,rho)
    x=mvrnorm(n,mu=rep(0,p), Sigma = Sigma)
  }
  if(case == 3){#undecay correlation
    corrmat = diag(rep(1-rho, p)) + matrix(rho, p, p)
    cholmat = chol(corrmat)
    x = matrix(rnorm(n*p, 0, 1), n, p)
    x = x%*%cholmat
  }
  else if(case == 4){#nonlinear relation
    true.beta = rep(0,p)
    true.beta[1:non.zeros] = 1
    x=matrix(rnorm(n*p),n,p)
    x[,2]=-1/3*x[,1]^3+rnorm(n)
  }
  S = 1:non.zeros
  #Check irrespresentable condition
  C=1/n*t(x)%*%x
  C11=C[S,S]
  C21=C[-S,S]
  s=sign(true.beta[S])
  IR=all(abs(C21 %*% solve(C11) %*% s) -1<0.0001)
  data=list(x=x,n=n,p=p,non.zeros=non.zeros,var_t=S,true.beta=true.beta,IR=IR )
  return(data)
}
#every sublist should have same structure
flatten = function(alist,K){
  newlist=alist[[1]]
  if (K >1){
    for(j in 2:K){
      newlist=Map(rbind,newlist,alist[[j]])
    }
  }
  return(newlist)
}



eva.lop=function(p, newlist, delta,K, n,B,sigmasq_t){
  hat_prob= as.numeric(newlist$hat_prob)
  hat_M_index= as.numeric(newlist$hat_M_index)

  hat_logp = sapply(delta, function(x)logp_delta(hat_prob,x))
  prob = mean(hat_M_index)
  logP = sapply(delta, function(x)logp_delta(prob,x))
  Elogp = colMeans(hat_logp)
  RB = (Elogp - logP)/logP * 100
  sd_hat_logp = apply(hat_logp, 2, sd)
  CV = sd_hat_logp/ abs(Elogp)
  result = data.frame(p,n, B, Elogp, logP, RB, CV,delta,sigmasq_t)
  return(result)
}

eva.cs = function(p, width,cover,sis_select, K, alpha,n,B,sigmasq_t){
  ##result of confidence sets
  covering_prob = colSums(cover)/K
  #length = colSums(2^width)/K
  #ACP=covering_prob/ length
  AW= colSums(width)/K
  C.level= (1-alpha)*100
  SIS_rate= mean(sis_select)
  result = data.frame(p, n, B, C.level, covering_prob, AW, SIS_rate)
  return(result)
}




#beta is a vector of all fixed effects
rep_3=function(case,rho,family, sigmasq_t, B, alpha,delta, n, p, penalty,tune,screen,iter){
  data=Generate.X(case,n,p,rho)
  X= data$x
  var_t= data$var_t
  predy=data$x %*% data$true.beta
  #generate data under Mopt
  Y <- Generate.Y(predy, sigmasq_t,n, family)
  data$y=Y
  predictors=1:p
  sis_select=T
  if(screen==TRUE){
  if(family=="gam"){
    model_ini=greedINIS(data,quant=0.95,gnum = 20)
    predictors= model_ini$ISISind
    X= X[,predictors]
  }else{
  nsis=n-1
  model_ini = SIS(X, Y, family = family,nsis=nsis,iter = iter)
  predictors= model_ini$ix0
  X= X[,predictors]
  }
  if (!all(var_t %in% predictors)){
    sis_select=F
  }
  }
  pn=0
  if (is.numeric(predictors))pn=length(predictors)
  result=list()
  #NMCS
  if (pn>0){
  nested = NMCS(Y, X, family, B, alpha, delta, penalty,tune)
  width =sapply(1: length(alpha), function(x)nested$mcs[[x]]$width)
  cover = sapply(1: length(alpha),function(x)all(
    predictors[nested$mcs[[x]]$LB] %in%c(0,var_t))
    &all(var_t %in% predictors[nested$mcs[[x]]$UB]))
  #LogP
  hat_M_index = as.numeric(all(predictors[nested$hat_M$var_M] %in% c(0,var_t))&
                             all(var_t %in% predictors[nested$hat_M$var_M]))
  hat_prob =  nested$hat_prob
  #MCB
  mcb.width=NULL
  mcb.cover=NULL
  #if (family=="gaussian"){
  #MCB=mcb.cv(X, Y, B, penalty=penalty, tune=tune)$mcbframe
  #mcb.width = sapply(alpha, function(x)MCB$width[which(MCB$bcr>=1-x)[1]])
  #mcb.cover = sapply(alpha,function(x)all(
  #  predictors[as.numeric(MCB$lbm[which(MCB$bcr>=1-x)[1]][[1]])] %in%  c(0,var_t))
  #   &all(var_t  %in% predictors[as.numeric(MCB$ubm[which(MCB$bcr>=1-x)[1]][[1]])]))
  #}
  result=list(width=width,cover=cover,hat_M_index=hat_M_index,
              hat_prob=hat_prob, sis_select=sis_select,
              mcb.width=mcb.width, mcb.cover=mcb.cover,pn=pn,
              CRE=nested$CRE,
              CREb=nested$CREb)
  }else{
    result=NULL####SKIP THIS RUN since no predictors survived!
  }
  return(result)
}

#Include non-zero intercept
simu_3.PL = function(case,family, B,n,p,sigmasq_t, alpha,
                     K, cores,r.names=NULL,
                     rho=0,penalty="scad", tune="bic",delta=0.0001,screen=FALSE,iter=TRUE){
  if (is.null(r.names))r.names=family
  Boot = mclapply(1:K, function(x)rep_3(case,rho, family, sigmasq_t, B,
                                        alpha,delta, n, p, penalty,tune,screen,iter),mc.cores=cores)
  valid=which(sapply(Boot, function(x)!is.null(x)))
  Bootn=Boot[valid]
  k=length(Bootn)
  #save(Bootn,file = "Bootn.RData")
  newlist = flatten(Bootn,k)
  ##result of confidence sets
  eva_cs=list()
  eva_logp=list()
  eva_cs= eva.cs(p, newlist$width,newlist$cover,newlist$sis_select,K,
                 alpha,n,B,sigmasq_t)
  eva_logp=eva.lop(p, newlist, delta,K, n, B, sigmasq_t)
  mcb_cs = rep(NULL,length(alpha))
  if (family=="gaussian"){
  mcb.cs= eva.cs(p, newlist$mcb.width,newlist$mcb.cover,newlist$sis_select,K, alpha,n,B,sigmasq_t)
  eva_cs= cbind(eva_cs, method="MCB",mcb.cs)
  }

  write.table(cbind(case=case,rho=rho,penalty=penalty, eva_cs,iter=iter,k=k,
                    CRE=mean(newlist$CRE),CREb=mean(newlist$CREb)), file =paste(r.names,"_cs.csv",sep=""), sep=",",
              row.names=FALSE, col.names=FALSE,append = TRUE)
  write.table(cbind(case=case,rho=rho,penalty=penalty, eva_logp,iter=iter,k=k), file =paste(r.names,"_logp.csv",sep=""), sep=",",
              row.names=FALSE, col.names=FALSE,append = TRUE)
  return(list(eva_cs, eva_logp))
}


cores=24
#setwd("~/NMAC3")
#------------------------------LogP simulation-------------------
nseq=c(50,100,200,500)
B=500
K=200
alpha=c(0.05,seq(0.1,0.7,by=0.1))
sigmasq_t=1
delta=1e-4
family="gaussian"
r.names= "comp1"
#for (n in nseq){
#p=floor(n^(3/4))
#simu_3.PL(1,family=family,B ,n,p, sigmasq_t,
          #alpha,K, cores,penalty="lasso", tune="bic",delta=delta,r.names=r.names)
#simu_3.PL(1,family=family,B ,n,p, sigmasq_t,
         # alpha, K, cores,penalty="scad", tune="bic",delta=delta,r.names=r.names)
#simu_3.PL(1,family=family,B ,n,p, sigmasq_t,
          #alpha, K, cores,penalty="adlasso", tune="bic",delta=delta,r.names=r.names)
#}


#-----------------NMCS simulations----------------
p=1000
n=200
B=200
K=100
iter=FALSE
tune="bic"
r.names=NULL
simu_3.PL(1,"gaussian", B,n,p, sigmasq_t,
 alpha, K, cores,penalty="adlasso", tune=tune,screen = T,iter=iter,r.names=r.names)
simu_3.PL(2,"gaussian", B,n,p, sigmasq_t,
          alpha, K, cores,rho=0.5,penalty="adlasso", tune=tune,screen = T,iter=iter,r.names=r.names)
simu_3.PL(3,"gaussian", B,n,p, sigmasq_t,
         alpha, K, cores, rho=0.5,penalty="adlasso", tune=tune,screen = T,iter=iter,r.names=r.names)

simu_3.PL(1,"binomial", B,n,p, sigmasq_t,
          alpha, K, cores, penalty = "adlasso",tune = tune,screen = T,iter=iter)
simu_3.PL(2,"binomial", B,n,p, sigmasq_t,
          alpha, K, cores, rho=0.5,penalty = "adlasso",tune = tune,screen = T,iter=iter)
simu_3.PL(3,"binomial", B,n,p, sigmasq_t,
          alpha, K, cores, rho=0.5,penalty = "adlasso",tune = tune,screen = T,iter=iter)
B=400
K=100
simu_3.PL(4,"gam", 200, n,p,3,
          alpha, K, cores, screen = T)



