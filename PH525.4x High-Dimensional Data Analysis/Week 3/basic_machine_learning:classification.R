#Conditional Expectation Exercises

n = 1000
y = rbinom(n,1,0.25)
##proportion of ones Pr(Y)
sum(y==1)/length(y)
##expectaion of Y
mean(y)

#1
n = 10000
set.seed(1)
men = rnorm(n,176,7) #height in centimeters
women = rnorm(n,162,7) #height in centimeters
y = c(rep(0,n),rep(1,n))
x = round(c(men,women))
##mix it up
ind = sample(seq(along=y))
y = y[ind]
x = x[ind]

mean(y[x==176])

#2
xs = seq(160,178)
pr =sapply(xs,function(x0) mean(y[x==x0]))
plot(xs,pr)
abline(h=0.5)
abline(v=168)


#Smoothing Exercises

#1
set.seed(5)
N = 250
ind = sample(length(y),N)
Y = y[ind]
X = x[ind]

fit=loess(Y~X)
predict(fit,newdata=data.frame(X=168))

xs = seq(160,178)
Pr =sapply(xs,function(x0) mean(Y[X==x0]))
plot(xs,Pr)
fitted=predict(fit,newdata=data.frame(X=xs))
lines(xs,fitted)

#2
set.seed(5)
B = 1000
N = 250
xs = seq(160,178)
plot(xs,xs,ylim=c(0,1),type="l")
res = replicate(B,{
  ind = sample(length(y),N)
  Y = y[ind]
  X = x[ind]
  fit=loess(Y~X)
  ##optional plots
  fitted=predict(fit,newdata=data.frame(X=xs))
  lines(xs,fitted)
  estimate = predict(fit,newdata=data.frame(X=168))
  return(estimate)
})
popsd(res)


#kNN and Cross Validation Exercises

#1
RNGkind(sample.kind = "Rounding")
library(GSE5859Subset)
data(GSE5859Subset)
y = factor(sampleInfo$group)
X = t(geneExpression)
out = which(geneAnnotation$CHR%in%c("chrX","chrY"))
X = X[,-out]
install.packages('caret')
library(caret)
set.seed(1)
createFolds(y,k=10)

#2
library(class)
library(genefilter)
m=8
k=5
ind = idx[[2]]

pvals = rowttests(t(X[-ind,]),factor(y[-ind]))$p.val
ind2 = order(pvals)[1:m]
predict=knn(X[-ind,ind2],X[ind,ind2],y[-ind],k=k)
sum(predict!=y[ind])

#3
m=8
k=5
result = sapply(idx,function(ind){
  pvals = rowttests(t(X[-ind,]),factor(y[-ind]))$p.val
  ind2 = order(pvals)[1:m]
  predict=knn(X[-ind,ind2],X[ind,ind2],y[-ind],k=k)
  sum(predict!=y[ind])
})
sum(result)/length(y)

#4
ms=2^c(1:11)
ks=seq(1,9,2)
params = expand.grid(k=ks,m=ms)

errors = apply(params,1,function(param){
  k =  param[1]
  m =  param[2]
  result = sapply(idx,function(ind){
    pvals = rowttests(t(X[-ind,]),factor(y[-ind]))$p.val
    ind2 = order(pvals)[1:m]
    predict=knn(X[-ind,ind2],X[ind,ind2],y[-ind],k=k)
    sum(predict!=y[ind])
  })
  sum(result)/length(y)
})
params[which.min(errors),]

#5
pvals = rowttests(t(X),factor(y))$p.val
errors = apply(params,1,function(param){
  k =  param[1]
  m =  param[2]
  result = sapply(idx,function(ind){
    ind2 = order(pvals)[1:m]
    predict=knn(X[-ind,ind2],X[ind,ind2],y[-ind],k=k)
    sum(predict!=y[ind])
  })
  sum(result)/length(y)
})
min(errors)

#6
y = factor(as.numeric(format( sampleInfo$date, "%m")=="06"))

errors = apply(params,1,function(param){
  k =  param[1]
  m =  param[2]
  result = sapply(idx,function(ind){
    pvals = rowttests(t(X[-ind,]),factor(y[-ind]))$p.val
    ind2 = order(pvals)[1:m]
    predict=knn(X[-ind,ind2],X[ind,ind2],y[-ind],k=k)
    sum(predict!=y[ind])
  })
  sum(result)/length(y)
})
min(errors)
