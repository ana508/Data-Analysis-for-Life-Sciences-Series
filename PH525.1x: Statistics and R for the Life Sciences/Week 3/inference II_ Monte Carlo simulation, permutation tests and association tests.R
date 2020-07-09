library(dplyr)
library(rafalib)
library(downloader)

###Monte Carlo Simulation

set.seed(1)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/mice_pheno.csv"
filename <- "mice_pheno.csv"
if (!file.exists(filename)) download(url,destfile=filename)
dat <- read.csv(filename)
controlPopulation <- dat %>% filter(Sex == "F" & Diet == "chow") %>% select(Bodyweight) %>% unlist

#We will build a function that automatically generates a t-statistic under the null hypothesis for a sample size of `n`.

ttestgenerator <- function(n) {
  #note that here we have a false "high fat" group where we actually
  #sample from the chow or control population. 
  #This is because we are modeling the null.
  cases <- sample(controlPopulation,n)
  controls <- sample(controlPopulation,n)
  tstat <- (mean(cases)-mean(controls)) / 
    sqrt( var(cases)/n + var(controls)/n ) 
  return(tstat)
}
ttests <- replicate(1000, ttestgenerator(10))

#We can see how distribution works with histogram
hist(ttests)

#We can see how normal approximation works qq plot
ttests <- replicate(1000, ttestgenerator(10))
qqnorm(ttests)
abline(0,1)

#when sample size is not big enough to follow a normal distribution, we can generate t distribution since it works begtter than a central limit theorem
#we will see much better approximation if we compare 2 conditions
mypar(1,2)
ttests <- replicate(1000, ttestgenerator(3))
qqnorm(ttests)
abline(0,1)
ps <- (seq(0,999)+0.5)/1000
qqplot(qt(ps,df=2*3-2),ttests,xlim=c(-6,6),ylim=c(-6,6))
abline(0,1)

#The t-distribution is a much better approximation in this case, but it is still not perfect. This is due to the fact that the original data is not that well approximated by the normal distribution.
mypar()
qqnorm(controlPopulation)
qqline(controlPopulation)


#The technique we used to motivate random variables and the null
#distribution was a type of Monte Carlo simulation. We had access to
#population data and generated samples at random. In practice, we do
#not have access to the entire population. When we want to
#use Monte Carlo simulations in practice, it is much more typical to
#assume a parametric distribution and generate a population from
#this, which is called a _parametric simulation_. This means that we take
#parameters estimated from the real data (here the mean and the standard
#deviation), and plug these into a model (here the normal simulation. 

#For the case of weights, we could use our knowledge that mice typically weigh 24 grams with a SD of about 3.5 grams, and that the distribution is approximately normal, to generate population data
controls<- rnorm(5000, mean=24, sd=3.5) 

#After we generate the data, we can then repeat the exercise above. We no longer have to use the `sample` function since we can re-generate random normal numbers. The `ttestgenerator` function therefore can be written as follows
ttestgenerator <- function(n, mean=24, sd=3.5) {
  cases <- rnorm(n,mean,sd)
  controls <- rnorm(n,mean,sd)
  tstat <- (mean(cases)-mean(controls)) / 
    sqrt( var(cases)/n + var(controls)/n ) 
  return(tstat)
}


###Monte Carlo Exercises

#1
set.seed(1)
x <- rnorm(5)
t <- sqrt(5)*mean(x)/sd(x)

#2
set.seed(1)
N <- 5
B<- 1000
tstats <- replicate(B,{
  X <- rnorm(N)
  sqrt(N)*mean(X)/sd(X)
})
mean(tstats>2)

#3
Ns<-seq(5,30,5)
B <- 1000
mypar(3,2)
LIM <- c(-4.5,4.5)
for(N in Ns){
  ts <- replicate(B, {
    X <- rnorm(N)
    sqrt(N)*mean(X)/sd(X)
  })
  ps <- seq(1/(B+1),1-1/(B+1),len=B)
  qqplot(qt(ps,df=N-1),ts,main=N,
         xlab="Theoretical",ylab="Observed",
         xlim=LIM, ylim=LIM)
  abline(0,1)
} 
#For which sample sizes does the approximation best work?
  #The approximations are spot on for all sample sizes.

#4
Ns<-seq(5,30,5)
B <- 1000
mypar(3,2)
LIM <- c(-4.5,4.5)
for(N in Ns){
  ts <- replicate(B,{
    x <- rnorm(N)
    y <- rnorm(N)
    t.test(x,y, var.equal = TRUE)$stat
  })
  ps <- seq(1/(B+1),1-1/(B+1),len=B)
  qqplot(qt(ps,df=2*N-2),ts,main=N,
         xlab="Theoretical",ylab="Observed",
         xlim=LIM, ylim=LIM)
  abline(0,1)
}
#For which sample sizes does the approximation best work?
  #The approximations are spot on for all sample sizes.

#5
mypar()
set.seed(1)
N <- 15
B <- 10000
tstats <- replicate(B,{
  X <- sample(c(-1,1), N, replace=TRUE)
  sqrt(N)*mean(X)/sd(X)
})
ps=seq(1/(B+1), 1-1/(B+1), len=B) 
qqplot(qt(ps,N-1), tstats, xlim=range(tstats))
abline(0,1)
#The population data is not normal thus the theory does not apply.
#We check with a Monte Carlo simulation. The qqplot shows a large tail. 
#Note that there is a small but positive chance that all the X are the same.
##In this case the denominator is 0 and the t-statistics is not defined

#6
set.seed(1)
N <- 1000
B <- 10000
tstats <- replicate(B,{
  X <-  sample(c(-1,1), N, replace=TRUE)
  sqrt(N)*mean(X)/sd(X)
})
qqnorm(tstats)
abline(0,1)
#With N=1000, CLT kicks in and the t-statistic is approximated with normal 0,1
##Furthermore, t-distribution with df=999 and normal are practically the same.

#7
set.seed(1)
Ns <- seq(5,45,5)
mypar(3,3)
for(N in Ns){
  medians <- replicate(10000, median ( rnorm(N) ) )
  title <- paste("N=",N,", avg=",round( mean(medians), 2) , ", sd*sqrt(N)=", round( sd(medians)*sqrt(N),2) )
  qqnorm(medians, main = title )
  qqline(medians)
}
##there is an asymptotic result that says SD is sqrt(N*4*dnorm(0)^2)

###Permutations Exercises
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile=filename)
babies <- read.table("babies.txt", header=TRUE)
bwt.nonsmoke <- babies %>% filter(smoke==0) %>% select(bwt) %>% unlist 
bwt.smoke <- babies %>% filter(smoke==1) %>% select(bwt) %>% unlist

#1
N <- 10
set.seed(1)
nonsmokers <- sample(bwt.nonsmoke , N)
smokers <- sample(bwt.smoke , N)
obs <- mean(smokers) - mean(nonsmokers)

dat <- c(smokers,nonsmokers)
shuffle <- sample( dat )
smokersstar <- shuffle[1:N]
nonsmokersstar <- shuffle[(N+1):(2*N)]
mean(smokersstar)-mean(nonsmokersstar)

set.seed(1)
null <- replicate(1000, {
  shuffle <- sample( dat )
  smokersstar <- shuffle[1:N]
  nonsmokersstar <- shuffle[(N+1):(2*N)]
  mean(smokersstar)-mean(nonsmokersstar)
})
( sum( abs(null) >= abs(obs)) +1 ) / ( length(null)+1 ) 
##we add the 1s to avoid p-values=0 but we also accept:
( sum( abs(null) >= abs(obs)) ) / ( length(null) )

#2
set.seed(1)
obs <- median(smokers) - median(nonsmokers)
null <- replicate(1000, {
  shuffle <- sample( dat )
  smokersstar <- shuffle[1:N]
  nonsmokersstar <- shuffle[(N+1):(2*N)]
  median(smokersstar)-median(nonsmokersstar)
})
( sum( abs(null) >= abs(obs)) +1 ) / ( length(null)+1 )
## As before we add 1 to avoid p-value of 0 but we also accept
( sum( abs(null) >= abs(obs)) ) / ( length(null) )

###Association Tests Exercises

d = read.csv("assoctest.csv")

#1
tab <- table(d$allele, d$case)
chisq.test(tab)$statistic

#2
fisher.test(tab)$p.value