if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install('genefilter')
library(genefilter)

###Error Rates and Procedures Examples

library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename <- "femaleControlsPopulation.csv"
if (!file.exists(filename)) download(url,destfile=filename)

#Monte Carlo simulation:
set.seed(1)
population = unlist( read.csv("femaleControlsPopulation.csv") )
alpha <- 0.05
N <- 12
m <- 10000
pvals <- replicate(m,{
  control = sample(population,N)
  treatment = sample(population,N)
  t.test(treatment,control)$p.value
})
sum(pvals < 0.05) ##This is R/V

#Adding delta when null hyphotheses is false:
alpha <- 0.05
N <- 12
m <- 10000
p0 <- 0.90 ##10% of diets work, 90% don't
m0 <- m*p0
m1 <- m-m0
nullHypothesis <- c( rep(TRUE,m0), rep(FALSE,m1))
delta <- 3

#Generating alternative hypotheses and see whether or not we ignore hte null hypotheses:
set.seed(1)
calls <- sapply(1:m, function(i){
  control <- sample(population,N)
  treatment <- sample(population,N)
  if(!nullHypothesis[i]) treatment <- treatment + delta
  ifelse( t.test(treatment,control)$p.value < alpha, 
          "Called Significant",
          "Not Called Significant")
})

null_hypothesis <- factor( nullHypothesis, levels=c("TRUE","FALSE"))
table(null_hypothesis,calls)

#Demonstrating how V and S change evety time since they are random variables:
B <- 10 ##number of simulations
VandS <- replicate(B,{
  calls <- sapply(1:m, function(i){
    control <- sample(population,N)
    treatment <- sample(population,N)
    if(!nullHypothesis[i]) treatment <- treatment + delta
    t.test(treatment,control)$p.val < alpha
  })
  cat("V =",sum(nullHypothesis & calls), "S =",sum(!nullHypothesis & calls),"\n")
  c(sum(nullHypothesis & calls),sum(!nullHypothesis & calls))
})


###Error Rates and Procedures Exercises

#1
#What does  f(x)  look like?
  #the identity line

#2
B<-1000
minpval <- replicate(B, min(runif(8793,0,1))<0.05)
mean(minpval>=1)

#3
B=10000
cutoffs = 10^seq(-7,-4,0.1) ##we know it has to be small
prob = sapply(cutoffs,function(cutoff){
  minpval =replicate(B, min(runif(8793,0,1))<=cutoff)
  mean(minpval>=1)
})
cutoffs[which.min(abs(prob-0.05))]


###Bonferroni Correction Exercises

#1
alphas <- seq(0,0.25,0.01)
par(mfrow=c(2,2))
for(m in c(2,10,100,1000)){
  plot(alphas,alphas/m - (1-(1-alphas)^(1/m)),type="l")
  abline(h=0,col=2,lty=2)
}

#2
set.seed(1)
B <- 10000
m <- 8793
alpha <- 0.05
pvals <- matrix(runif(B*m,0,1),B,m)
k <- alpha/m
mistakes <- rowSums(pvals<k) 
mean(mistakes>0)

#3
set.seed(1)
B <- 10000
m <- 8793
alpha <- 0.05
pvals <- matrix(runif(B*m,0,1),B,m)
k <- (1-(1-alpha)^(1/m))
mistakes <- rowSums(pvals<k) 
mean(mistakes>0)


###False Discovery Rate

set.seed(1)
population = unlist( read.csv("femaleControlsPopulation.csv") )
alpha <- 0.05
N <- 12
m <- 10000
pvals <- replicate(m,{
  control = sample(population,N)
  treatment = sample(population,N)
  t.test(treatment,control)$p.value
})
alpha <- 0.05
N <- 12
m <- 10000
p0 <- 0.90 ##10% of diets work, 90% don't
m0 <- m*p0
m1 <- m-m0
nullHypothesis <- c( rep(TRUE,m0), rep(FALSE,m1))
delta <- 3

##Vectorised approach:
controls <- matrix(sample(population, N*m, replace=TRUE),nrow=m)

##matrix with control data (rows are tests, columns are mice)
treatments <-  matrix(sample(population, N*m, replace=TRUE),nrow=m)

##add effect to 10% of them
treatments[which(!nullHypothesis),]<-treatments[which(!nullHypothesis),]+delta

##combine to form one matrix
dat <- cbind(controls,treatments)
g <- factor( c(rep(0,N),rep(1,N)) )
pvals <- rowttests(dat,g)$p.value < alpha
sum(pvals  <= alpha/m)

#FDR:
calls <- pvals <= alpha
R <- sum(calls)
V <- sum(nullHypothesis & calls)
Q <- ifelse(R>0,V/R,0)

#Doing 10000 tests 1000 times:
set.seed(1)
g <- factor( c(rep(0,N),rep(1,N)) )
B <- 1000 ##number of simulations
Qs <- replicate(B,{
  ##matrix with control data (rows are tests, columns are mice)
  controls <- matrix(sample(population, N*m, replace=TRUE),nrow=m)
  
  ##matrix with control data (rows are tests, columns are mice)
  treatments <-  matrix(sample(population, N*m, replace=TRUE),nrow=m)
  
  ##add effect to 10% of them
  treatments[which(!nullHypothesis),]<-treatments[which(!nullHypothesis),]+delta
  
  ##combine to form one matrix
  dat <- cbind(controls,treatments)
  
  calls <- rowttests(dat,g)$p.value < alpha
  R=sum(calls)
  Q=ifelse(R>0,sum(nullHypothesis & calls)/R,0)
  return(Q)
})

fdr <- p.adjust(pvals, method="fdr")
mypar(1,1)
plot(pvals,fdr,log="xy")

#Defining procedure that has a false discovery rate 0.05:
alpha <- 0.05
B <- 1000 ##number of simulations. We should increase for more precision
res <- replicate(B,{
  controls <- matrix(sample(population, N*m, replace=TRUE),nrow=m)
  treatments <-  matrix(sample(population, N*m, replace=TRUE),nrow=m)
  treatments[which(!nullHypothesis),]<-treatments[which(!nullHypothesis),]+delta
  dat <- cbind(controls,treatments)
  pvals <- rowttests(dat,g)$p.value 
  ##then the FDR
  calls <- p.adjust(pvals,method="fdr") < alpha
  R=sum(calls)
  Q=ifelse(R>0,sum(nullHypothesis & calls)/R,0)
  return(c(R,Q))
})
Qs <- res[2,]
mypar(1,1)
hist(Qs) ##Q is a random variable, this is its distribution
FDR=mean(Qs)
print(FDR)

###FDR Exercises

library(devtools)
library(rafalib)
install_github("genomicsclass/GSE5859Subset")
library(GSE5859Subset)
data(GSE5859Subset)

#1
g <- factor(sampleInfo$group)
pvals <- rowttests(geneExpression,g)$p.value
sum(pvals<0.05)

#2
k <- 0.05/length(pvals)
sum(pvals<k)

#3
fdr <- p.adjust(pvals,method="fdr")
sum(fdr<0.05)

#4
library(qvalue)
res <- qvalue(pvals)
qvals <- res$qvalues
sum(qvals<0.05)

#5
qvalue(pvals)$pi0

#6
plot(qvalue(pvals)$qvalue/p.adjust(pvals,method="fdr"))
abline(h=qvalue(pvals)$pi0,col=2)

hist(pvals,breaks=seq(0,1,len=21))
expectedfreq <- length(pvals)/20 #per bin
abline(h=expectedfreq*qvalue(pvals)$pi0,col=2,lty=2)

#7
n <- 24
m <- 8793
mat <- matrix(rnorm(n*m),m,n)
delta <- 2
positives <- 500
mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta

set.seed(1)
positives <- 500
g <- factor(rep(c(0,1),each=12))
result <- replicate(B,{
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals = rowttests(mat,g)$p.val
  ##Bonferroni
  FP1 <- sum(pvals[-(1:positives)]<=0.05/m)  
  FP1
})
mean(result/(m-positives))

#8
set.seed(1)
result <- replicate(B,{
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals = rowttests(mat,g)$p.val
  ##Bonferroni
  FP1 <- sum(pvals[-(1:positives)]<=0.05/m)  
  FN1 <- sum(pvals[1:positives]>0.05/m)
  c(FP1,FN1)
})

mean(result[2,]/(positives))

#9
result <- replicate(B,{
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals = rowttests(mat,g)$p.val
  ##Bonferroni
  FP1 <- sum(pvals[-(1:positives)]<=0.05/m)  
  FN1 <- sum(pvals[1:positives]>0.05/m)
  #p.adjust q-value
  qvals1 = p.adjust(pvals,method="fdr")
  FP2 <- sum(qvals1[-(1:positives)]<=0.05)
  c(FP1,FN1,FP2)
})

mean(result[3,]/(m-positives))

#10
result <- replicate(B,{
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals = rowttests(mat,g)$p.val
  ##Bonferroni
  FP1 <- sum(pvals[-(1:positives)]<=0.05/m)  
  FN1 <- sum(pvals[1:positives]>0.05/m)
  #p.adjust q-value
  qvals1 = p.adjust(pvals,method="fdr")
  FP2 <- sum(qvals1[-(1:positives)]<=0.05)
  FN2 <- sum(qvals1[1:positives]>0.05)
  c(FP1,FN1,FP2,FN2)
})

mean(result[4,]/(positives))

#11
result <- replicate(B,{
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals = rowttests(mat,g)$p.val
  ##Bonferroni
  FP1 <- sum(pvals[-(1:positives)]<=0.05/m)  
  FN1 <- sum(pvals[1:positives]>0.05/m)
  #p.adjust q-value
  qvals1 = p.adjust(pvals,method="fdr")
  FP2 <- sum(qvals1[-(1:positives)]<=0.05)
  FN2 <- sum(qvals1[1:positives]>0.05)
  #qvalue q-value
  qvals2 = qvalue(pvals)$qvalue
  FP3 <- sum(qvals2[-(1:positives)]<=0.05)
  c(FP1,FN1,FP2,FN2,FP3)
})

mean(result[5,]/(m-positives))

#12
result <- replicate(B,{
  mat <- matrix(rnorm(n*m),m,n)
  mat[1:positives,1:(n/2)] <- mat[1:positives,1:(n/2)]+delta
  pvals = rowttests(mat,g)$p.val
  ##Bonferroni
  FP1 <- sum(pvals[-(1:positives)]<=0.05/m)  
  FN1 <- sum(pvals[1:positives]>0.05/m)
  #p.adjust q-value
  qvals1 = p.adjust(pvals,method="fdr")
  FP2 <- sum(qvals1[-(1:positives)]<=0.05)
  FN2 <- sum(qvals1[1:positives]>0.05)
  #qvalue q-value
  qvals2 = qvalue(pvals)$qvalue
  FP3 <- sum(qvals2[-(1:positives)]<=0.05)
  FN3 <- sum(qvals2[1:positives]>0.05)  
  c(FP1,FN1,FP2,FN2,FP3,FN3)
})

mean(result[6,]/(positives))