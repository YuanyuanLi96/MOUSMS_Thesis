#Nested MAC Algorithms
##X do not have intercept
library(rlist)
library(MASS)
library(parallel)
library(glmnet)
library(ncvreg)
library(gamsel)
library(emdbook)

loglik_lam = function(X, y, beta, family) {
  K = dim(beta)[2]
  n =length(y)
  link = cbind(1, X) %*% beta
  yrep = repmat(y, 1, K)
  if (family == "gaussian")
    return(n*log(apply((yrep - link)^2, 2, mean))/2)
  if (family == "poisson")
    return(apply(exp(link) - yrep * link, 2, sum))
  if (family == "binomial")
    return(apply(log(1 + exp(link)) - yrep * link, 2, sum))
}


getdf = function(coef.beta,eps=1e-6) {
  df=apply(abs(coef.beta) > eps, 2, sum)+1
  return(df)
}

repmat = function(X, m, n) {
  ## R equivalent of repmat (matlab)
  X = as.matrix(X)
  mx = dim(X)[1]
  nx = dim(X)[2]
  matrix(t(matrix(X, mx, nx * n)), mx * m, nx * n, byrow = TRUE)
}


Generate.Y <- function(predy, sigmasq, n, family) {
  if (family %in% c("gaussian","gam")){
    Y = predy + rnorm(n,0, sd = sqrt(sigmasq))
  }
  if (family=="binomial"){
    prob=1/(1+exp(-predy))
    Y = rbinom(n=n,size=1, prob=prob)
  }
  return(Y)
}

half.order=function(beta,rowind,colind, eps=1e-6, all=TRUE){
  ##length(rowind),length(colind)>=1
  if (length(rowind)==1){return(rowind)}
  if (length(colind)==1){
    return(rowind[order(beta[rowind],decreasing = TRUE)])
  }
  sel_beta=beta[rowind,colind]
  enterall = factor(apply(sel_beta,1,function(x)which(x>eps)[1]))
  #break tie by their value
  variable.order=NULL
  for (l in levels(enterall)){
    vnam= which(enterall==l)
    variable.order= c(variable.order,vnam[order(sel_beta[vnam,as.numeric(l)],
                                                decreasing = TRUE)])
  }
  if(all){variable.order=c(variable.order, which(is.na(enterall)))}
  return(rowind[variable.order])
}

#Order of predictors by lasso solution path
#beta is the norm of coefficients matrix
enter.order = function(beta,index.min,eps=1e-6){
  beta=abs(beta)
  hat_M=which(beta[,index.min]> 1e-6)
  hat_MC=which((1:nrow(beta))%in%hat_M==FALSE)
  variable.order=half.order(beta,hat_M,1:index.min)
  raw_order= half.order(beta,1:nrow(beta),1:index.min,all=FALSE)
  CE=all(variable.order%in% raw_order)&all(raw_order%in% variable.order)
  if (ncol(beta)>=(index.min+1)){
    variable.order=c(variable.order, half.order(beta, hat_MC,(index.min+1):ncol(beta)))
  }else{
    variable.order=c(variable.order,hat_MC)
  }                
  return(list(EO=variable.order,CE=CE))
}
#enter.order(beta,index.min = which(fit1$lambda==fit1$lambda.1se))


####################Model selection using different tuning methods
cf_norm=function(dat, nbasis){
  m=nrow(dat)/nbasis
  norm=NULL
  for (i in 1:m){
    data_i=dat[nbasis*(i-1)+1:nbasis,]
    norm=rbind(norm, apply(data_i,2, function(x)sum(x^2)))
  }
  return(norm)
}
cv.gam<- function(y, x,family,eps=0.01){
  n=length(y)
  #lambda = lseq(50, 0.0001, 100)
  nbase=6
  bases=pseudo.bases(x,degree=nbase,df=3)
  gamsel.cv=cv.gamsel(x,y,bases=bases)
  alpha_norm = (gamsel.cv$gamsel.fit$alphas)^2
  betas_norm = cf_norm(gamsel.cv$gamsel.fit$betas,nbase)
  coef_norm = alpha_norm+betas_norm
  #coef_norm = alpha_norm
  index=gamsel.cv$index.1se
  hat_M=which(coef_norm[,index]>eps)
  predy=predict(gamsel.cv$gamsel.fit,x,index=index,type="response")
  sigmasq_est <- mean((predy-y)^2)
  eoinfo= enter.order(coef_norm,index.min = index)
  var.order=eoinfo$EO
  CE=eoinfo$CE
  len= length(hat_M)
  return(list(len= len, var.order=var.order,var_M=hat_M, predy=predy
              ,sigmasq=sigmasq_est,
              CE=CE))
}


sel.method=function(y, x, family, penalty,tune,eps=1e-6){
  n= length(y)
  if (family=="gam"){
    return(cv.gam(y, x,family))
  }
  if(tune=="cv"){
    if(penalty=="scad"){
      cv.fit = cv.ncvreg(x, y, family = family, penalty = "SCAD",nfolds=5)
      betas= cv.fit$fit$beta[-1,]
      lambda.ind = min(which(cv.fit$cve<cv.fit$cve[ cv.fit$min]+cv.fit$cvse[ cv.fit$min]))
      coef.beta = cv.fit$fit$beta[, lambda.ind]  
      # extract coefficients at a single value of lambda, including the intercept
    }else{
      lasso_init <- cv.glmnet(x,y, family=family,nfolds = 10)
      betas= lasso_init$glmnet$beta
      lambda.ind=lasso_init$lambda.1se
      coef.beta =as.numeric(coef(lasso_init,s=lasso_init$lambda.1se))
      if (penalty=="adlasso"){
        penalty.factor <- abs(coef.beta[-1]+1/sqrt(n))^(-1)
        adalasso <- cv.glmnet(x,y, family=family, penalty.factor=penalty.factor,nfolds = 10)
        betas= adalasso$glmnet$beta
        lambda.ind=adalasso$lambda.1se
        coef.beta =as.numeric(coef(adalasso,s=lambda.ind))
      }
    }
  }else{
    if (tune=="aic")k=2
    if (tune=="bic")k=log(n)
    if(penalty=="scad"){
      reg.fit = ncvreg(x, y, family=family,penalty = "SCAD")
      betas = reg.fit$beta#include intercept
    }else{
      lasso_init <- glmnet(x,y, family=family,nfolds = 10)
      betas= coef(lasso_init)
      #df = lasso_init$df
      if (penalty=="adlasso"){
        dev = loglik_lam(x, y, betas, family=family)
        reg.df = getdf(betas[-1, , drop = FALSE])
        obj = 2*dev + k * reg.df
        lambda.ind = which.min(obj)
        coef.beta = betas[, lambda.ind]
        penalty.factor <- abs(coef.beta[-1]+1/sqrt(n))^(-1)
        adalasso <- glmnet(x,y, family=family, penalty.factor=penalty.factor,nfolds = 10)
        betas= rbind(adalasso$a0,adalasso$beta)
      }
    }
    dev = loglik_lam(x, y, betas, family=family)
    reg.df = getdf(betas[-1, , drop = FALSE])
    obj = 2*dev + k * reg.df
    lambda.ind = which.min(obj)
    coef.beta = betas[, lambda.ind]
    betas=betas[-1,]
  }
  #fit glm
  eoinfo=enter.order(abs(betas),index.min = lambda.ind)
  var.order=eoinfo$EO
  CE=eoinfo$CE
  hat_M = which(abs(coef.beta[-1])> eps)
  len = length(hat_M)
  if (len==0){
    fit=glm(y~1, family=family)
  }else{
    Xi = x[,hat_M]
    fit=glm(y~Xi, family=family)
  }
  sigmasq_est = sum(fit$residuals^2)/n
  coef.final=rep(0,dim(x)[2])
  coef.final[hat_M]=fit$coefficients[-1]
  return(list(len=len, var.order=var.order,var_M=hat_M, beta=coef.final,
              predy=predict(fit,type = "link"), sigmasq=sigmasq_est,
              CE=CE))

}
#prob is an element
logp_delta=function(prob, delta){
  logP=rep(0,length(prob))
  for (i in 1:length(prob)){
    if (isTRUE(all.equal(prob[i],0))){
      logP[i]=log(1-delta)
    }else if (isTRUE(all.equal(prob[i], 1))){
      logP[i]=log(delta)
    }else{
      logP[i]=log(1-prob[i])
    }
  }
  return(logP)
}


##################Algorithms for constructing sequence of model confidence sets
#Calculate the Model confidences bounds and their corresponding coverage rate

trunc_bounds = function(x,p){
  if(x<0)x=0
  if(x>p)x=p
  return(x)
}

#alpha ia a given value
MCS_func <- function(bounds, alpha,order_M, len_M,p,eps=1e-6){
  k = which(bounds[1,]- 1+ alpha > -eps)[1]
  BCP=bounds[1,k]
  w=k-1
  j=bounds[2,k]
  #cat("coverage_prob",result$coverage_prob,"\n")
  lower.loc = trunc_bounds(len_M-w+j,p)
  upper.loc = trunc_bounds(len_M+j,p)
  width = upper.loc- lower.loc
  if(upper.loc==0){
    LB=0
    UB=0
  }else{
    if(lower.loc ==0){
      LB = 0
    }else{
      LB =order_M[1:lower.loc]
    }
    UB = order_M[1:upper.loc]
  }
  return(list(alpha= alpha, BCP=BCP, width = width,LB= LB,UB= UB))
}


coverp=function(var_M, len, order, p, B){
  bool=order %in%var_M
  largest_lb=ifelse(!all(bool==TRUE),which(bool==FALSE)[1]-1,p)#largest lower bound
  smallest_ub=ifelse(!all(bool==FALSE),tail(which(bool==TRUE),1), 0)#smallest upper bound
  dl=max(0,len-largest_lb)
  dr=max(0, smallest_ub-len)
  distance=dl+dr
  CP= lapply(0:(2*p),function(x)rep(0,x+1))
  for (w in distance:(2*p)){
    CP[[w+1]][(dr+1):(w-dl+1)]=1/B
  }
  return(CP)
}



#Complexity: O(B)
Conf <- function(var_M, len, order.matrix, p, B){
  CP=coverp(var_M, len[1],order.matrix[1,],p,B)
  for (i in 2:B){
    CP=Map("+",CP,coverp(var_M, len[i],order.matrix[i,],p,B))
  }
  maxCP=sapply(CP, function(x)c(max(x), which.max(x)-1))
  #return 2*2p matrix, rows: CP, j*
  return(maxCP)
}



##################################Main function
#X must include intercept term in the first column
NMCS =function(Y, X, family, B, alpha, delta, penalty,tune){
  n = nrow(X)
  p = ncol(X)#number of predictors exclude intercept

  #select model using original data
  hat_M = sel.method(y=Y, x=X, family=family, penalty,tune)#a list including \hat{M}, \hat{\psi} and rank
  order_M = hat_M$var.order
  sigmasq_M = hat_M$sigmasq
  len_M =hat_M$len
  var_M = hat_M$var_M
  predy_M =hat_M$predy
  CE=hat_M$CE
  ##bootstrap procedure
  order_boot = matrix(0, nrow = B, ncol = p)
  len_boot = rep(0,B)
  CEb=rep(0, B)
  for (b in 1:B){
    #generate data under hat_M
    y_b <- Generate.Y(predy_M, sigmasq_M, n, family)

    #select model using bootstrap data
    hat_M_boot = sel.method(y=y_b, x=X,family=family, penalty,tune)#a list including \hat{M}, \hat{\psi} and rank
    order_boot[b,] = hat_M_boot$var.order
    len_boot[b] = hat_M_boot$len
    CEb[b]=hat_M_boot$CE
  }

  result =list()
  ##confidence sets
  Conf.I = Conf(var_M, len_boot, order_boot, p, B)
  result$mcs = lapply(alpha, function(x)MCS_func(Conf.I, x,order_M, len_M,p))
  ##prob of \hat{M}_b = \hat{M}
  result$hat_prob =  Conf.I[1,1]
  result$hat_logP = logp_delta(result$hat_prob, delta)
  result$hat_M= hat_M
  result$CE=CE
  result$CEb=mean(CEb)
  return(result)
}

#MAC.nested(Y, X , B, alpha, delta, sel.bic, Conf.I2)


