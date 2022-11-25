###This R script is modified from source codes of the R package-mcb
###sel.method function is defined in Nested_glm.R to include more flexible tunning methods.

#Transfrom a list of variable name denoted selection results to 0-1 matrix result
f01 <- function(object, full.var, p){

  matrix.01 <- data.frame(matrix(0, length(object), p))
  colnames(matrix.01) <- full.var
  for(i in 1:length(object))
  {
    matrix.01[i, (full.var %in% object[[i]])] <- 1
  }
  return(matrix.01)
}

#Calculate the Model confidences bounds and their corresponding freq rate
CI <- function(var.list, var.matrix, p, B)
{
  full.var <- colnames(var.matrix)
  colsum <- apply(var.matrix, 2, sum)
  order <- order(colsum, decreasing = T)
  freq <- vector(length = p+1)
  freq[1] <- 0
  lower <- vector(mode="list", length = p+1)
  upper <- vector(mode="list", length = p+1)
  for(i in 0:p)
  {
    cap <- vector(length=p+1)
    cap[1] <- 0
    for(j in 0:(p-i))
    {
      if (j==0 & i!=0)
      {
        uppertest <- full.var[order[1:i]]
        for(r in 1:length(var.list))
        {
          if(all(var.list[[r]] %in% uppertest)) cap[j+1]<-cap[j+1]+1
        }
      }else{
        if(j!=0){
          lowtest <- full.var[order[1:j]]
          uppertest <- full.var[order[1:(j+i)]]
          for(r in 1:length(var.list))
          {
            if(all(all(lowtest %in% var.list[[r]]),all(var.list[[r]] %in% uppertest))) cap[j+1]<-cap[j+1]+1
          }
        }
        if (j == 0){
          for (r in 1:length(var.list)){
            if (identical(var.list[[r]], character(0))) cap[j+1] <- cap[j+1] + 1
          }
        }
      }
    }
    freq[i+1] <- max(cap)/B
    maxlocation <- which.max(cap)
    if(maxlocation==1)
    {
      if (i != 0){
        lower[[i+1]] <- 0
        upper[[i+1]] <- full.var[order[1:i]]
      } else if (i == 0){
        lower[[1]] <- 0
        upper[[1]] <- 0
      }
    }else{
      lower[[i+1]] <- full.var[order[1:(maxlocation-1)]]
      upper[[i+1]] <- full.var[order[1:(maxlocation-1+i)]]
    }
  }
  result <- list(freq=freq,lower=lower,upper=upper)
  return(result)
}




# GET the residual bootsrap variable selection results by using adaptive lasso and
# you can give a specified lambda
Boot_Sel<-function(data, r, p, full.var,penalty,tune,family="gaussian",eps=1e-6){
  var.instances <- vector(mode="list",length=r)
  x=data[,1:p]
  y=data[,p+1]
  beta_est= sel.method(y, x, family, penalty,tune)$beta
  constant <- as.matrix(x) %*% beta_est
  res_original <- y - constant
  res_after_center <- res_original - mean(res_original)
  for(j in 1:r) {
    ind=sample(1:nrow(x),nrow(x),replace=T)
    new_response <- constant + res_after_center[ind]
    boot.data <- cbind(x, new_response)
    beta_est_boot= sel.method(boot.data[,p+1], as.matrix(boot.data[,1:p]), family, penalty,tune)$beta
    var.instances[[j]]<-full.var[abs(beta_est_boot)>eps]
  }
  return(var.instances)
}


# GET residual boostrap variable selection models with stepwise BIC
RES.BOOT.CI5 <- function(x, p, r, lmbd, full.var){
  if(lmbd==''){
    var<-vector(mode="list",length=r)
    fit<-regsubsets(y~.,data=x,method="seqrep",nvmax=p)
    if (names(coef(fit,which.min(summary(fit)$bic)))[1] == '(Intercept)'){
      beta <- vector(mode = 'numeric', length = p)
      beta[full.var%in%names(coef(fit,which.min(summary(fit)$bic)))] <- coef(fit,which.min(summary(fit)$bic))[-1]
    } else {
      beta <- vector(mode = 'numeric', length = p)
      beta[full.var%in%names(coef(fit,which.min(summary(fit)$bic)))] <- coef(fit,which.min(summary(fit)$bic))
    }
    res_original <- x$y - as.matrix(x[,1:p]) %*% beta
    res_after_center <- res_original - mean(res_original)
    constant <- as.matrix(x[,1:p]) %*% beta
    for(j in 1:r) {
      ind=sample(1:nrow(x),nrow(x),replace=T)
      new_response <- constant + res_after_center[ind]
      boot.data <- cbind(x[,1:p], new_response)
      colnames(boot.data)[p+1] <- "y"
      tem <- regsubsets(y~.,data=boot.data,method="seqrep",nvmax=p)
      if (names(coef(tem,which.min(summary(tem)$bic)))[1] == '(Intercept)'){
        var[[j]]<-names(coef(tem,which.min(summary(tem)$bic)))[-1]
      } else {
        var[[j]]<-names(coef(tem,which.min(summary(tem)$bic)))
      }
    }
    return(var)
  }else{
    var<-vector(mode="list",length=r)
    fit<-regsubsets(y~.,data=x,method="seqrep",nvmax=p)
    if (names(coef(fit,which.min(abs(summary(fit)$bic - lmbd))))[1] == '(Intercept)'){
      beta <- vector(mode = 'numeric', length = p)
      beta[full.var%in%names(coef(fit,which.min(abs(summary(fit)$bic - lmbd))))] <- coef(fit,which.min(summary(fit)$bic))[-1]
    } else {
      beta <- vector(mode = 'numeric', length = p)
      beta[full.var%in%names(coef(fit,which.min(abs(summary(fit)$bic - lmbd))))] <- coef(fit,which.min(summary(fit)$bic))
    }
    res_original <- x$y - as.matrix(x[,1:p]) %*% beta
    res_after_center <- res_original - mean(res_original)
    constant <- as.matrix(x[,1:p]) %*% beta
    for(j in 1:r) {
      ind=sample(1:nrow(x),nrow(x),replace=T)
      new_response <- constant + res_after_center[ind]
      boot.data <- cbind(x[,1:p], new_response)
      colnames(boot.data)[p+1] <- "y"
      tem <- regsubsets(y~.,data=boot.data,method="seqrep",nvmax=p)
      if (names(coef(tem,abs(which.min(summary(tem)$bic - lmbd))))[1] == '(Intercept)'){
        var[[j]]<-names(coef(tem,abs(which.min(summary(tem)$bic - lmbd))))[-1]
      } else {
        var[[j]]<-names(coef(tem,abs(which.min(summary(tem)$bic - lmbd))))
      }
    }
    return(var)
  }
}


mcb.cv <- function(x, y, B=200, penalty = 'lasso', tune="bic",level=0.95,eps=1e-6){
  # methods = c('adlasso', 'lasso', 'scad')
  n = nrow(x)

  data <- cbind(x,y)
  p <- dim(data)[2]-1
  full.var=1:p
  r = B

  # adaptivelasso
    var_adaptivelasso <- Boot_Sel(data, r, p, full.var,penalty = penalty,
                                  tune=tune)
    var_01_ada_lasso<-f01(var_adaptivelasso, full.var, p)
    result <- CI(var_adaptivelasso, var_01_ada_lasso, p, r)

  all_result <- list()

  # ggplot2 - mucplot
  df <- data.frame(result$freq)
  df$x <- seq(0,1,length.out = p+1)
  df$x <- round(df$x,digits = 3)
  muc <- ggplot(df, aes(x=df$x, y=result$freq)) + geom_line() + labs(x = "w/p") + labs(y = "freq") + labs(colour = "Method") + ylim(0,1) + scale_x_continuous(breaks = df$x)
  all_result$mucplot <- muc

  fit_freq = c(result$freq)[result$freq - level >-eps]
  best_fit = which.min(fit_freq - level) + p+1 - sum(result$freq - level > -eps)
  mcb = list()
  mcb$lbm <- result$lower[[best_fit]]
  mcb$ubm <- result$upper[[best_fit]]
  mcb$bcr <- result$freq[best_fit]
  all_result$mcb <- mcb

  mcbframe <- data.frame(lbm = matrix(result$lower), bcr = result$freq, ubm = matrix(result$upper))
  mcbframe$width <- c(0:p)
  mcbframe <- mcbframe[,c('width','lbm','bcr','ubm')]
  all_result$mcbframe <- mcbframe

  return(all_result)
}


