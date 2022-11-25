source("real_data.R")
data(ndd1)
y=ndd1[,1]
x=ndd1[,-c(1,2)]
x=apply(x,2,scale)
x=as.matrix(x)
y=as.numeric(y)
n=nrow(x)
alpha=c(0.05, 0.1, 0.2, 0.3)
B=500
iter=FALSE
#r0=SIS_NMCS(y,x, B=B, alpha=alpha, screen = TRUE,penalty = "lasso", tune="bic",iter = iter)
r1=SIS_NMCS(y,x, B=B, alpha=alpha, screen = TRUE,penalty = "adlasso", tune="bic",iter = iter)
#r2=SIS_NMCS(y,x, B=B, alpha=alpha, screen = TRUE,penalty = "scad", tune="bic",iter = iter)
r3=SIS_NMCS(y,x, B=B, alpha=alpha, screen = TRUE,penalty = "adlasso", tune="cv",iter = iter)
save.image(file='ndd1.RData')
#save results
setwd("result")
write.table(apply(cbind(r0$nmcs.r$MCS.frame, r0$MCB$MCS.frame),2,as.character),
            file ="ndd1_result.csv", sep=",",
            row.names=FALSE, col.names=FALSE)
write.table(apply(cbind(r1$nmcs.r$MCS.frame, r1$MCB$MCS.frame),2,as.character),
            file ="ndd1_result.csv", sep=",",
            row.names=FALSE, col.names=FALSE,appen=TRUE)
write.table(apply(cbind(r2$nmcs.r$MCS.frame, r2$MCB$MCS.frame),2,as.character),
            file ="ndd1_result.csv", sep=",",
            row.names=FALSE, col.names=FALSE,appen=TRUE)
write.table(apply(cbind(r3$nmcs.r$MCS.frame, r3$MCB$MCS.frame),2,as.character),
            file ="ndd1_result.csv", sep=",",
            row.names=FALSE, col.names=FALSE,appen=TRUE)


#South African heart disease
SAH= read.table("SAH.txt",
                sep=",",head=T,row.names=1)
y=SAH[,10]
x=SAH[,-10]
glm(chd~., family = "binomial",data=SAH)
xn=apply(x, 2, as.numeric)
xn=apply(xn, 2, scale)
xn[,5]=x[,5]
xn[,5]=as.factor(xn[,5])
head(xn)
B=500
alpha=c(0.05, 0.1, 0.2, 0.3)
result=list()
result[[1]]=SIS_NMCS(y,xn,"binomial", B, alpha, penalty="lasso", tune="aic",iter=FALSE)
result[[2]]=SIS_NMCS(y,xn,"binomial", B, alpha, penalty="lasso", tune="bic",iter=FALSE)
result[[3]]=SIS_NMCS(y,xn,"binomial", B, alpha, penalty="lasso", tune="cv",iter=FALSE)
result[[4]]=SIS_NMCS(y,xn,"binomial", B, alpha, penalty="adlasso", tune="aic",iter=FALSE)
result[[5]]=SIS_NMCS(y,xn,"binomial", B, alpha, penalty="adlasso", tune="bic",iter=FALSE)
result[[6]]=SIS_NMCS(y,xn,"binomial", B, alpha, penalty="adlasso", tune="cv",iter=FALSE)

for (i in 1:6){
  r0=result[[i]]
  write.table(apply(r0$nmcs.r$MCS.frame,2,as.character),
              file ="sah_result.csv", sep=",",
              row.names=FALSE, col.names=FALSE, append = TRUE)
}



