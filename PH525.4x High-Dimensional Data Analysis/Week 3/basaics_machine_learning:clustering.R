#Hierarchichal Clustering Exercises

RNGkind()

#1
set.seed(1)
m = 10000
n = 24
x = matrix(rnorm(m*n),m,n)
colnames(x)=1:n
hc <- hclust(dist(t(x)))
plot(hc)

#2
set.seed(1)
m = 10000
n = 24
nc = replicate(100,{
  x = matrix(rnorm(m*n),m,n)
  hc = hclust( dist( t(x)))
  length(unique(cutree(hc,h=143)))
})
plot(table(nc))
library(rafalib)
popsd(nc)


#K-means Exercises

#1
library(GSE5859Subset)
data(GSE5859Subset)

mds <- cmdscale(dist(t(geneExpression)))
set.seed(10)
result <- kmeans(t(geneExpression),5)
mypar(1,1)
plot(mds,bg=result$cl,pch=21)
table(sampleInfo$group,result$cluster)
table(sampleInfo$date,result$cluster)
##looks better if we re-order:
table(sampleInfo$date,result$cluster)[,c(4,1,5,3,2)]
