###Poisson Example from RNA-seq

##Simulated gene expression data:
N = 10000    # number of genes
lambdas = 2^seq(1,16, len=N)    # these are the true abundances of genes
y = rpois(N, lambdas)    # note that the null hypothesis is true for all genes
x = rpois(N, lambdas)
ind = which(y>0 & x>0)    # make sure no 0s due to ratio and log

library(rafalib)
mypar()
splot(log2(lambdas), log2(y/x), subset=ind)

##Real gene expression data:
library(parathyroidSE)
data(parathyroidGenesSE)
se = parathyroidGenesSE
x = assay(se)[,23]
y = assay(se)[,24]

ind = which(y>0 & x>0)    # make sure no 0s due to ratio and log
splot((log2(x)+log2(y))/2, log2(x/y), subset=ind)


###Statistical Models Exercises

#1
dbinom(2,4,0.49)

#2
dbinom(4,10,0.49)

#3
1-pbinom(10,20,0.4)

#4
1 - dbinom(0, 189000000, 1/175223510)

#5
1 - pbinom(1, 189000000, 1/175223510)

#6
pbinom(9,20,0.4)-pbinom(7,20,0.4)

#7
b <- (9 - 20*.4)/sqrt(20*.4*.6)
a <- (7 - 20*.4)/sqrt(20*.4*.6)
pnorm(b)-pnorm(a)

#8
exact <- pbinom(450,1000,0.4)-pbinom(350,1000,0.4)
b <- (450 - 1000*.4)/sqrt(1000*.4*.6)
a <- (350 - 1000*.4)/sqrt(1000*.4*.6)
approx <- pnorm(b)-pnorm(a)
abs(exact-approx)

#9
Ns <- c(1000000)
ps <- c(0.001)
k = 1:10
exact = dbinom(k,Ns,ps)
plot(k,exact, main = "exact", type = "l")

a <- (k+0.5 - Ns*ps)/sqrt(Ns*ps*(1-ps))
b <- (k-0.5 - Ns*ps)/sqrt(Ns*ps*(1-ps))
approx = pnorm(a) - pnorm(b)
plot(k, approx, main="approx", type = "l")

#10
N <- 189000000
p <- 1/175223510
dbinom(2,N,p)
a <- (2+0.5 - N*p)/sqrt(N*p*(1-p))
b <- (2-0.5 - N*p)/sqrt(N*p*(1-p))
pnorm(a) - pnorm(b)
dpois(2,N*p)

1 - ppois(1, N*p)


###Maximum Likelihood Estimate (MLE)
datadir="http://www.biostat.jhsph.edu/bstcourse/bio751/data"
x=read.csv(file.path(datadir,"hcmv.csv"))[,2]
breaks=seq(0,4000*round(max(x)/4000),4000)
tmp=cut(x,breaks)
counts=table(tmp)

l<-function(lambda){
  ls <- dpois(counts,lambda,log=TRUE)
  return(sum(ls))
}  
lambdas<-seq(3,7,len=100)
ls <- exp(sapply(lambdas,l))
plot(lambdas,ls,type="l")
mle=optimize(l,c(0,10),maximum=TRUE)
abline(v=mle$maximum)


###MLE Exercises

#1
library(devtools)
install_github("genomicsclass/dagdata")
library(dagdata)
data(hcmv)
library(rafalib)
mypar()
plot(locations,rep(1,length(locations)),ylab="",yaxt="n")

breaks=seq(0,4000*round(max(locations)/4000),4000)
tmp=cut(locations,breaks)
counts=as.numeric(table(tmp))
hist(counts)

probs <- dpois(counts,4)
likelihood <- prod(probs)
likelihood

#more convenient way of counting MLE:
logprobs <- dpois(counts,4,log=TRUE)
loglikelihood <- sum(logprobs)
loglikelihood

func <- function(lambda, counts){
  logprobs <- dpois(counts,lambda,log=TRUE)
  loglikelihood <- sum(logprobs)
  loglikelihood
}
lambdas <- seq(0,15,len=300)
mles <- sapply(lambdas, func, counts)
plot(lambdas, mles)

lambdas[which.max(mles)]

#2
breaks=seq(0,4000*round(max(locations)/4000),4000)
tmp=cut(locations,breaks)
counts=as.numeric(table(tmp))
binLocation=(breaks[-1]+breaks[-length(breaks)])/2
plot(binLocation,counts,type="l",xlab=)

binLocation[which.max(counts)]

#3
max(counts)

#4
lambda = mean(counts[ - which.max(counts) ])
pval = 1 - ppois(13,lambda)
pval

#5
#From the question above, we obtain a p-value smaller than 0.001 for a count of 14. Why is it problematic to report this p-value as strong evidence of a location that is different?
  #We selected the highest region out of 57 and need to adjust for multiple testing.

#6
0.05/57

#7
ps <- (seq(along=counts) - 0.5)/length(counts)
lambda <- mean( counts[ -which.max(counts)])
poisq <- qpois(ps,lambda)
qqplot(poisq,counts)
abline(0,1)


###Models for Variance Exercises

install_github("genomicsclass/tissuesGeneExpression")
library(tissuesGeneExpression)
data("tissuesGeneExpression")
library(genefilter)
y = e[,which(tissue=="endometrium")]

#1
library(genefilter)
s2 <- rowVars(y)
library(rafalib)
mypar(1,2)
qqnorm(s2)
qqline(s2)
##To see the square root transformation does not help much:
qqnorm(sqrt(s2))
qqline(sqrt(s2))
#Which statement is true? (pick one)
  #The normal distribution is not a useful approximation here: the left tail is over estimated and the right tail is underestimated.

#2
BiocManager::install('limma')
library(limma)
estimates <- fitFDist(s2,14)

#3
ps <- (seq(along=s2)-0.5)/length(s2)
theoretical<- qf(ps,14,estimates$df2)*estimates$scale 
LIM <- sqrt( range(c(theoretical,s2)) )
mypar(1,2)
qqplot(sqrt( theoretical ), sqrt( s2 ),ylim=LIM,xlim=LIM)
abline(0,1)
##close up excluding the upper 5%
K <- sqrt( quantile(s2,0.95) )
qqplot( sqrt( theoretical ), sqrt( s2 ),ylim=c(0,K),xlim=c(0,K))
abline(0,1)
