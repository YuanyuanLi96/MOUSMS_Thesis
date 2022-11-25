#source.all("~/paper/LMM/rpql/R")
####X does not include intercept, intercept will be fitted automatically and excluded for selection.
rpqlpath <- function(y, X, Z, id, family = gaussian(), hybrid.est=FALSE,trial.size = 1, lambda, pen.type = "lasso", tune="bic1", start = NULL, cov.groups = NULL, 
                    pen.weights = NULL, offset = NULL, intercept = FALSE, save.data = FALSE, 
                    control = list(tol = 1e-6, maxit = 100, trace = FALSE, restarts = 5, scad.a = 3.7, mcp.gamma = 2, lasso.lambda.scale = TRUE, seed = NULL), ...) {
  lambda <- as.matrix(lambda)
  if(!control$lasso.lambda.scale)
    warnings("Scaling factor for the second tuning parameter has been turned off for lasso and adaptive lassp penalties. 
             This may cause issues as the group coefficient penalty on the random effects is on a different scale to the single 
             coefficient penalty on the fixed effects.")
  if(ncol(lambda) == 1) {
    lambda[,1] <- sort(lambda[,1], decreasing =F)
  }
  
  
  control <- rpql:::fillin.control(control) 
  
  ## Do a starting fit for formatting reasons
  #start_fit <- rpql(y = y, X = X, Z = Z, id = id, family = family, lambda = lambda[1,], 
                   # pen.type = "lasso", start = start) 
  
  collect_ics <- matrix(Inf, nrow = nrow(lambda), ncol = 8) ## First column in start_fit$ics is DoF
  colnames(collect_ics) <- c("PQL Likelihood", "df","aic", "bic1","bic2","hic1", "hic2","hic3") 
  
  best_fits <- vector("list", 6); 
  names(best_fits) <- c("aic", "bic1","bic2","hic1", "hic2","hic3")
  #rm(start_fit)
  p = ncol(X)
  qr=sapply(Z, ncol)
  q = sum(qr)
  collect_betas <- matrix(0, nrow = p, ncol = nrow(lambda)) ## First column in start_fit$ics is DoF
  collect_b <- matrix(0, nrow = q, ncol = nrow(lambda)) ## First column in start_fit$ics is DoF
  which_min=rep(1, length(best_fits))
  names(which_min) = names(best_fits)
  for(l1 in 1:nrow(lambda)) {
    message("Onto ", l1)
    if(l1 == 1) 
      prev_fit <- start
    
    new_fit <-rpql(y = y, X = X, Z = Z, id = id, family = family, lambda = lambda[l1,], pen.type = pen.type, pen.weights = pen.weights, start = prev_fit, 
                    cov.groups = cov.groups, hybrid.est = hybrid.est, offset = offset, intercept = intercept, save.data = save.data, control = control)
    if(is.null(new_fit)){
      next
    }
    
    if(is.null(best_fits[[1]])){
      for(k in 1:length(best_fits)) 
        best_fits[[k]] <- new_fit
    }else{
      for(l2 in 1:length(best_fits)) {
        if(new_fit$ics[l2+1] < best_fits[[l2]]$ics[l2+1]){ 
          best_fits[[l2]] <- new_fit
          which_min[l2] = l1
        }
      }
    }
    prev_fit <- new_fit
    collect_ics[l1,] <- c(new_fit$logLik1, new_fit$ics)
    collect_betas[new_fit$nonzero.fixef,l1] = 1
    start.ind=0
    for(e in 1:length(qr)) { 
      collect_b[start.ind + new_fit$nonzero.ranef[[e]],l1]= 1
      start.ind= start.ind+ qr[e]
    }
  }
  allv= rbind(collect_betas,collect_b)
  var.order = enter.order(fliplr(allv))
  out <- list(all.best.fits= best_fits, best.fit = best_fits[[tune]], collect.ics = collect_ics, 
              lambda = lambda, control = control, allv=allv,
              var_M=which(allv[,which_min[tune]]>0.01),var.order=var.order)
  return(out)
}



