#SVD Exercises

library(tissuesGeneExpression)
data(tissuesGeneExpression)
s = svd(e)
signflips = sample(c(-1,1),ncol(e),replace=TRUE)
signflips
newu= sweep(s$u,2,signflips,FUN="*")
newv= sweep(s$v,2,signflips,FUN="*" )
all.equal( s$u %*% diag(s$d) %*% t(s$v), newu %*% diag(s$d) %*% t(newv))

#1
s = svd(e)
m = rowMeans(e)

cor(s$u[,1],m)

#2
newmeans = rnorm(nrow(e)) ##random values we will add to create new means
newe = e+newmeans ##we change the means
sqrt(crossprod(e[,3]-e[,45]))
sqrt(crossprod(newe[,3]-newe[,45]))
y = e - rowMeans(e)
s = svd(y)
resid = y - s$u %*% diag(s$d) %*% t(s$v)
max(abs(resid))
x=matrix(rep(c(1,2),each=5),5,2)
x
x*c(1:5)
sweep(x,1,1:5,"*")

#4
z = s$d * t(s$v)
sqrt(crossprod(e[,3]-e[,45]))
sqrt(crossprod(y[,3]-y[,45]))
sqrt(crossprod(z[,3]-z[,45]))

abs(sqrt(crossprod(e[,3]-e[,45])) - sqrt(crossprod(z[1:2,3]-z[1:2,45])))

#5
ks = 1:189
realdistance = sqrt(crossprod(e[,3]-e[,45]))
approxdistances = sapply(ks,function(k){
  sqrt(crossprod(z[1:k,3,drop=FALSE]-z[1:k,45,drop=FALSE] )) 
})
percentdiff = 100*abs(approxdistances - realdistance)/realdistance
plot(ks,percentdiff) ##take a look
min(ks[which(percentdiff < 10)])

#6
distances = sqrt(apply(e[,-3]-e[,3],2,crossprod))
approxdistances = sqrt(apply(z[1:2,-3]-z[1:2,3],2,crossprod))
plot(distances,approxdistances) ##take a look
cor(distances,approxdistances,method="spearman")


#MDS Exercises

#1
y = e - rowMeans(e)
s = svd(y)
z = s$d * t(s$v)
library(rafalib)
ftissue = factor(tissue)
mypar(1,1)
plot(z[1,],z[2,],col=as.numeric(ftissue))
legend("topleft",levels(ftissue),col=seq_along(ftissue),pch=1)
d = dist(t(e))
mds = cmdscale(d)

cor(z[1,],mds[,1])

#2
cor(z[2,],mds[,2])

#4
s = svd(geneExpression-rowMeans(geneExpression))
z = s$d * t(s$v)

which.max(cor(sampleInfo$g,t(z)))

#5
s = svd(geneExpression-rowMeans(geneExpression))
z = s$d * t(s$v)
max(cor(sampleInfo$g,t(z)))

#6
which.max(cor(sampleInfo$g,t(z))[-1]) + 1

#7
month = format( sampleInfo$date, "%m")
month = factor( month)

which.max(cor( as.numeric(month), t(z)))
max(cor( as.numeric(month), t(z)))

#8
result = split(s$u[,6],geneAnnotation$CHR)
result = result[ which(names(result)!="chrUn") ]
boxplot(result,range=0)
boxplot(result,range=0,ylim=c(-0.025,0.025))
medians = sapply(result,median)
names(result)[ which.max(abs(medians)) ]
