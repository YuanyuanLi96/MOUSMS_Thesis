library(miceadds)
#source.all("~/paper/nmcs/R")
library(Matrix)
#X does not include intercept, since response is centered.
#Z: A list with each element being a model matrix for a set of random effects.
#An intercept column must be included if a random intercept is desired.
#start: An object of class "lmerMod" or "glmerMod", obtained when fitting a (G)LMM
#using the lmer and glmer functions in the lme4 package.
NMCS_rpql =function(y, X, Z, id, family=gaussian(), lambda,B=200, alpha=0.05, delta=1e-4,
               penalty="adl",tune="bic1",start,pen.weights = pen.weights){
  n = nrow(X)
  p = ncol(X)#number of predictors
  qs = sapply(Z, ncol)
  q=sum(qs)

  #select model using original data
  hat_M = rpqlpath(y, X, Z, id, family=family, lambda=lambda, pen.type =penalty,
                   tune=tune, pen.weights = pen.weights, start = start)#a list including \hat{M}, \hat{\psi} and rank
  order_M = hat_M$var.order
  var_M= hat_M$var_M
  len_M= length(var_M)
  X_sel= X[,var_M[var_M<=p]]
  ran_e=hat_M$best.fit$nonzero.ranef[[1]]
  if (length(ran_e)==0){
    fit=glm(y ~ X_sel-1,family=family)
  }else{
  Z_sel = Z[[1]][,hat_M$best.fit$nonzero.ranef[[1]]]
  g1=id[[1]]
  fit= glmer(y ~ X_sel-1 + (Z_sel-1 | g1),family=family,nAGQ = 0)
  }
  fit.est= glmm.est(fit)
  ##bootstrap procedure
  order_boot = matrix(0, nrow = B, ncol = p+q)
  len_boot = rep(0,B)
  for (b in 1:B){
    #generate data under hat_M
    simu <-gendat.glmm(id, X=X_sel, beta=fit.est$fixef, list(Z_sel), D=fit.est$ran.cov, family = family, phi= fit.est$phi)
    y_b=simu$y
    #select model using bootstrap data
    hat_M_boot = rpqlpath(y_b, X, Z, id, family, lambda=lambda, 
                          pen.type =penalty,tune=tune, pen.weights = pen.weights, start = start)
    order_boot[b,] = hat_M_boot$var.order
    len_boot[b] = length(hat_M_boot$var_M)
  }
  
  result =list()
  ##confidence sets
  Conf.I = nmcs:::Conf(var_M, len_boot, order_boot, p+q, B)
  result$mcs = lapply(alpha, function(x)nmcs:::MCS_func(Conf.I, x,order_M, len_M,p+q))
  ##prob of \hat{M}_b = \hat{M}
  result$hat_prob =  Conf.I[1,1]
  result$hat_logP = nmcs:::logp_delta(result$hat_prob, delta)
  result$hat_M= hat_M
  return(result)
}

NMCS_mlasso =function(y, X, Z, id, family=gaussian(), lambda,B=200, alpha=0.05, delta=1e-4){
  n = nrow(X)
  p = ncol(X)#number of predictors exclude intercept
  q=ncol(Z)+1
  
  #select model using original data
  hat_M = Pen.LME(y, X, Z, id, t.fracs = lambda) 
  order_M = hat_M$var.order
  var_M= hat_M$var_M
  len_M= length(var_M)
  X_sel= hat_M$X_sel
  Z_sel = hat_M$Z_sel
  fit= hat_M$final
  
  ##bootstrap procedure
  order_boot = matrix(0, nrow = B, ncol = p+q)
  len_boot = rep(0,B)
  for (b in 1:B){
    #generate data under hat_M
    y_b <-simulate(fit)[,1]
    #select model using bootstrap data
    hat_M_boot = Pen.LME(y_b, X, Z, id, t.fracs = lambda) 
    order_boot[b,] = hat_M_boot$var.order
    len_boot[b] = length(hat_M_boot$var_M)
  }
  
  result =list()
  ##confidence sets
  Conf.I = Conf(var_M, len_boot, order_boot, p+q, B)
  result$mcs = lapply(alpha, function(x)MCS_func(Conf.I, x,order_M, len_M,p+q))
  ##prob of \hat{M}_b = \hat{M}
  result$hat_prob =  Conf.I[1,1]
  result$hat_logP = logp_delta(result$hat_prob, delta)
  result$hat_M= hat_M
  return(result)
}



