library(miceadds)

setwd("~/rcode_lmm")
source.all("rpql/R")
source("rpqlpath.R")
source("NMCSlmm.R")
# This is the first simulation example from the paper.
# The grid of fraction values is taken from 0.05 to 1 in steps of 0.025
require(mvtnorm)
require(tsutils)
library(lme4)
library(pracma)
library(nmcs)
library(parallel)

glmm.est= function(fit){
  if (class(fit)=="glm"){
    return(list(fixef=summary(fit)$coefficients[,1], ran.cov= NULL, phi=sigma(fit)))
  }
  vcovs= as.data.frame(VarCorr(fit))$vcov
  start.p=0
  Covr= list()
  ngrps.fit= ngrps(fit)
  dim.r= sapply(ranef(fit), ncol)
  for (i in 1:length(ngrps.fit)){
    q=dim.r[i]
    cov_r = diag(vcovs[start.p+1:q])
    cov_r[lower.tri(cov_r)] = vcovs[start.p+q+1:(q*(q-1)/2)]
    cov_r[upper.tri(cov_r)] = t(cov_r)[upper.tri(cov_r)]
    start.p = start.p+ q+q*(q-1)/2
    Covr[[i]]=cov_r
  }
  names(Covr)=names(ngrps.fit)
  sigmasq = vcovs[start.p+1]
  return(list(fixef=summary(fit)$coefficients[,1], ran.cov= Covr, phi=sigmasq))
}

eva.lop=function(newlist, delta,K, B){
  hat_prob= as.numeric(newlist$hat_prob)
  hat_M_index= as.numeric(newlist$hat_M_index)
  
  hat_logp = sapply(delta, function(x)nmcs:::logp_delta(hat_prob,x))
  prob = mean(hat_M_index)
  logP = sapply(delta, function(x)nmcs:::logp_delta(prob,x))
  Elogp = colMeans(hat_logp)
  RB = (Elogp - logP)/logP * 100
  sd_hat_logp = apply(hat_logp, 2, sd)
  CV = sd_hat_logp/ abs(Elogp)
  result = data.frame(B, delta, Elogp, logP, RB, CV)
  return(result)
}

eva.cs = function(width,cover, K, alpha,B){
  ##result of confidence sets
  covering_prob = colSums(cover)/K
  #length = colSums(2^width)/K
  #ACP=covering_prob/ length
  AW= colSums(width)/K
  C.level= (1-alpha)*100
  #SIS_rate= mean(sis_select)
  result = data.frame( B, C.level, covering_prob, AW)
  return(result)
}


flatten = function(alist,K){
  newlist=alist[[1]]
  if (K >1){
    for(j in 2:K){
      newlist=Map(rbind,newlist,alist[[j]])
    }
  }
  return(newlist)
}

run_func=function(n, m, B, lambda_seq,alpha, family=gaussian(),penalty="adl",tune="bic1",delta=delta){
id = kronecker(1:n, rep(1,m))
true.beta = c(1, 1, rep(0,7))
true.D = matrix(c(9,4.8,0.6,0, 4.8,4,1,0, .6,1,1,0, 0,0,0,0), nrow=4, ncol=4)
p=length(true.beta)+nrow(true.D)
predictors=1:p
var_t=which(abs(c(true.beta,diag(true.D)))>0.001)
X= matrix(runif(m*n*9,-2,2), nrow=m*n)
Z= cbind(1,matrix(runif(m*n*3,-2,2), nrow=m*n))
simy <- gendat.glmm(id = list(cluster = id), X = X, beta =true.beta,
                    Z = list(cluster = Z), D = list(cluster = true.D), phi = 1, family = family)
y=simy$y
if (family$family[1]=="gaussian"){
fit_satlme4 <- lmer(y ~ X-1 + (Z-1 | id))
}else{
fit_satlme4 <- glmer(y ~ X-1 + (Z-1 | id),family=family,nAGQ = 0)
}
fit_sat <- build.start.fit(fit_satlme4, gamma = 2)
pen.weights=fit_sat$pen.weights
#pen.weights$fixed[1]=0
run=NMCS_rpql(y = y, X = X, Z =list(g1=Z), id =list(g1=id),
              family=family, penalty=penalty,tune=tune,delta=delta,
                lambda = lambda_seq,start=fit_sat,B=B,alpha=alpha,
              pen.weights=pen.weights)
width =sapply(1: length(alpha), function(x)run$mcs[[x]]$width)
cover = sapply(1: length(alpha),function(x)all(
  predictors[run$mcs[[x]]$LB] %in%c(0,var_t))
  &all(var_t %in% predictors[run$mcs[[x]]$UB]))
#LogP
hat_M_index = as.numeric(all(predictors[run$hat_M$var_M] %in% c(0,var_t))&
                           all(var_t %in% predictors[run$hat_M$var_M]))
hat_prob =  run$hat_prob
result=list(width=width,cover=cover,hat_M_index=hat_M_index,
            hat_prob=hat_prob)
return(result)
}


simu1 = function(n,m, B, alpha, K, cores,lambda,
                     family=gaussian(),penalty="adl",tune="bic1",delta=0.0001){

  #Boot=lapply(1:K, function(x)run_func(B, lambda, alpha,family=family,penalty=penalty,tune=tune,delta=delta))
  Boot = mclapply(1:K, function(x)run_func(n,m, B, lambda, alpha,family=family,
  penalty=penalty,tune=tune,delta=delta),mc.cores=cores)
  save(Boot,file = "Boot.RData")
  k=length(Boot)
  newlist = flatten(Boot,k)
  ##result of confidence sets
  eva_cs=list()
  eva_logp=list()
  eva_cs= eva.cs( newlist$width,newlist$cover,K,
                 alpha,B)
  eva_logp=eva.lop(newlist, delta,K, B)
  write.table(cbind(familly=family$family[1],n,m, K, eva_cs), file ="cs_wo.csv", sep=",",
              row.names=FALSE, col.names=FALSE,append = TRUE)
  write.table(cbind(familly=family$family[1], n, m, K,eva_logp), file ="logp_wo.csv", sep=",",
              row.names=FALSE, col.names=FALSE,append = TRUE)
  return(list(eva_cs, eva_logp))
}


alpha=c(0.05,0.1)
B=200
K=100
cores=20
lambda_seq <- lseq(1e-6,0.5,length=20)

#simu1(n=50,m=10, B=400,alpha,K=100,cores = cores,lambda = lambda_seq, family=binomial())
#simu1(n=50,m=20, B=400,alpha,K=100,cores = cores,lambda = lambda_seq, family=binomial())
#simu1(n=50,m=30, B=400,alpha,K=100,cores = cores,lambda = lambda_seq, family=binomial())
simu1(n=80,m=30, B=400,alpha,K=100,cores = cores,lambda = lambda_seq, family=binomial())

###gaussian
lambda_seq <- lseq(1e-6,0.5,length=30)
#simu1(n=30,m=5, B=200,alpha,K=K,cores = cores,lambda = lambda_seq)
#simu1(n=60,m=10, B=200,alpha,K=K,cores = cores,lambda = lambda_seq)
