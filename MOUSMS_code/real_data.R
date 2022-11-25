library(SIS)
source("Nested_glm.R")#nested method
source("mcb_cv.R")
library(ggplot2)
library(knitr)
output_NMCS=function(nmcs.r, alpha, predictors=NULL){
  if(is.null(predictors))predictors=1:length(nmcs.r$hat_M$var.order)
  hat_M= sort(predictors[nmcs.r$hat_M$var_M])
  nmcs.r.var=data.frame(CL=sapply(alpha,function(x)1-x),
                        bcp=sapply(1:length(alpha),function(x)
                          nmcs.r$mcs[[x]]$BCP),
                        width = sapply(1:length(alpha),function(x)
                          nmcs.r$mcs[[x]]$width),
                        lbm=matrix(lapply(1:length(alpha),function(x)
                          sort(predictors[nmcs.r$mcs[[x]]$LB]))),
                        ubm=matrix(lapply(1:length(alpha),function(x)
                          sort(predictors[nmcs.r$mcs[[x]]$UB]))))
  return(list(hat_M=hat_M, MCS.frame=nmcs.r.var))
}

output_MCB=function( MCB,alpha,predictors){
  index= sapply(alpha, function(x)which(MCB$bcr>=1-x)[1])
  MCB.var=data.frame(CL=sapply(alpha,function(x)1-x),
                     bcp=sapply(index,function(x)
                       MCB$bcr[[x]]),
                     width=sapply(index,function(x)
                       MCB$width[[x]]),
                     lbm=matrix(lapply(index,function(x)
                       sort(predictors[as.numeric(MCB$lbm[[x]])]))),
                     ubm=matrix(lapply(index,function(x)
                       sort(predictors[as.numeric(MCB$ubm[[x]])]))))
  hat_M= sort(predictors[as.numeric(MCB$lbm[[1]])])
  return(list(hat_M=hat_M, MCS.frame=MCB.var))
}

#X include intercept as the first conlumn
SIS_NMCS=function(Y,X, family="gaussian", B, alpha, screen=FALSE,
                  delta=0.0001, penalty = "scad", tune="bic", nsis=NULL,iter=TRUE){
  n=nrow(X)
  if(is.null(nsis))nsis= n-1
  #use "ISIS-Lasso" method to screen predictors
  predictors=1:ncol(X)
  if(screen){
    model_ini = SIS(X, Y, family=family,iter=iter, nsis= nsis)
    predictors= model_ini$ix0
    X= X[,predictors]
  }
  #nmcs.r
  start_nmcs.r <- Sys.time()
  nmcs.r = NMCS(Y, X , family, B, alpha, delta, penalty,tune)
  end_nmcs.r <- Sys.time()
  nmcs.r_time=end_nmcs.r-start_nmcs.r
  nmcs.r.out=output_NMCS(nmcs.r,alpha,predictors)
  #LogP
  hat_prob = nmcs.r$hat_prob
  logP= logp_delta(hat_prob, delta)
  ##MCB
  MCB.out=NULL
  MCB_time=0
  start_MCB <- Sys.time()
  MCB=mcb.cv(X, Y, B, penalty = penalty,tune = tune)$mcbframe
  end_MCB = Sys.time()
  MCB.out=output_MCB(MCB,alpha,predictors)
  MCB_time=end_MCB-start_MCB
  return(list(nmcs.r=nmcs.r.out, MCB=MCB.out, predictors=predictors,
              nmcs.r_time=nmcs.r_time,
              MCB_time=MCB_time,logP=logP))
}

