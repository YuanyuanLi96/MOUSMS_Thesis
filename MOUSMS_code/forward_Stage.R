library(glmnet)
library(lars)
library(nmcs)
library(MASS)
set.seed(123)
p=20
n=200
rho=0.5
corrmat = diag(rep(1-rho, p)) + matrix(rho, p, p)
cholmat = chol(corrmat)
x = matrix(rnorm(n*p, 0, 1), n, p)
x = x%*%cholmat
apply(cor(x),1,sum)-1
non.zeros = 3
b=rep(1, non.zeros)
z = abs(rnorm(non.zeros)); signs = rbinom(non.zeros,1,0.4)
zero.entries = which(signs==0); signs[zero.entries] = -1
f1 = 5*log(n)/n^(1/2); f2 = z; b = (f1 + f2)*signs
true.beta = rep(0,p)
true.beta[1:non.zeros] = b
y=x%*% true.beta+rnorm(n)

####
set.seed(100)
n=200
p=20
x=matrix(rnorm(n*p),n,p)
x[,4]=-2*x[,1]+rnorm(n,sd=0.1)
#plot(x[,3],x[,1])
y=x[,1]-x[,2]+x[,3]+rnorm(n)
x_s=x[,1:3]
x_sc=x[,4:p]

abs(1/n*t(x_sc)%*%x_s%*%solve(1/n*t(x_s)%*% x_s))
####

fit1=cv.glmnet(x,y)
#system.time(glmnet(x,y))
plot(fit1$glmnet.fit) # plot coefficient paths
beta=fit1$glmnet$beta
which(abs(coef(fit1,s=fit1$lambda.1se))[-1]>0)
index.min=which(fit1$lambda==fit1$lambda.1se)
enter.order(fit1$glmnet$beta,index.min)


fit2=lars(x,y,type="lasso")
system.time(lars(x,y,type="lasso"))
plot(fit2)
enter.order(t(fit2$beta))

fit3=lars(x,y,type="lar")
plot(fit3)
enter.order(t(fit3$beta))

fit4=lars(x,y,type="forward.stagewise")
system.time(lars(x,y,type="forward.stagewise"))
plot(fit4)
enter.order(t(fit4$beta))



set.seed(1)
n <- 200; p <- 30
x = matrix(rnorm(n*p), ncol = p)
x[,2]=-1/3*x[,1]^3+rnorm(n)
mu = rowSums(x[,1:3])
y = mu + rnorm(n,sd=sqrt(3))
SNR=var(mu) /3
sqrt(SNR)
nbase=6
bases=pseudo.bases(x,degree=nbase,df=3)
gamsel.cv=cv.gamsel(x,y,bases=bases)
alpha_norm = (gamsel.cv$gamsel.fit$alphas)^2
betas_norm = cf_norm(gamsel.cv$gamsel.fit$betas,nbase)
coef_norm = alpha_norm+betas_norm
#coef_norm = alpha_norm
index=gamsel.cv$index.1se
hat_M=which(coef_norm[,index]>0)
hat_M
gamsel.cv$lambda.1se
summary(gamsel.cv$gamsel.fit)
plot(gamsel.cv$gamsel.fit,newx=x,which=1:10)
hat_M
p=10
rho=0.4
n=100
Sigma=ar1_cor(p,rho)
X=mvrnorm(n,mu=rep(0,p), Sigma = Sigma)
b=sample(c(1,-1),p,replace = TRUE)
B=solve((t(X)%*%X))
all(diag(b)%*% B %*% diag(b) %*% rep(1,p)>0)

b=sample(c(1,-1),n-1,replace = TRUE)
B=solve((t(X)%*%X))
all(diag(b)%*% B %*% diag(b) %*% rep(1,n-1)>0)

fit=adalasso
tLL <- fit$nulldev - deviance(fit)
k <- fit$df
n <- fit$nobs
AICc <- -tLL+2*k+2*k*(k+1)/(n-k-1)
AICc

BIC<-log(n)*k - tLL
plot(BIC)
