library(downloader) 
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename <- basename(url)
download(url, destfile=filename)
x <- unlist( read.csv(filename) )

# make averages5
set.seed(1)
n <- 1000
averages5 <- vector("numeric",n)
for(i in 1:n){
  X <- sample(x,5)
  averages5[i] <- mean(X)
}

# make averages50
set.seed(1)
n <- 1000
averages50 <- vector("numeric",n)
for(i in 1:n){
  X <- sample(x,50)
  averages50[i] <- mean(X)
}


#Normal Distribution Exercises #1
par(mfrow=c(1,2))
#They both look roughly normal, but with a sample size of 50 the spread is smaller.
hist(averages5, main = 'Averages of 5 mice samples')
hist(averages50, main = 'Averages of 50 mice samples')

#Normal Distribution Exercises #2
mean( averages50 < 25 & averages50 > 23)

#Normal Distribution Exercises #3
pnorm( (25-23.9) / 0.43)  - pnorm( (23-23.9) / 0.43) 

#For these exercises, we will be using the following dataset:
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/mice_pheno.csv"
filename <- basename(url)
download(url, destfile=filename)
dat <- read.csv(filename) 

#We will remove the lines that contain missing values:
dat <- na.omit( dat )

#Population, Samples, and Estimates Exercises #1
library(dplyr)
head(dat)
x <- dat %>% filter(Sex =='M' & Diet == 'chow') %>% select(Bodyweight) %>% unlist
mean(x)

#Population, Samples, and Estimates Exercises #2
install.packages('rafalib')
library(rafalib)
popsd(x)

#Population, Samples, and Estimates Exercises #3
set.seed(1)
X <- sample(x, 25)
mean(X)

#Population, Samples, and Estimates Exercises #4
y <- dat %>% filter(Sex =='M' & Diet == 'hf') %>% select(Bodyweight) %>% unlist
mean(y)

#Population, Samples, and Estimates Exercises #5
popsd(y)

#Population, Samples, and Estimates Exercises #6
set.seed(1)
Y <- sample(y, 25)
mean(Y)

#Population, Samples, and Estimates Exercises #7
abs((mean(y)-mean(x)) - (mean(Y)-mean(X)))

#Population, Samples, and Estimates Exercises #8
x <- dat %>% filter(Sex =='F' & Diet == 'chow') %>% select(Bodyweight) %>% unlist
set.seed(2)
X <- sample(x, 25)
y <- dat %>% filter(Sex =='F' & Diet == 'hf') %>% select(Bodyweight) %>% unlist
set.seed(2)
Y <- sample(y, 25)
abs((mean(y)-mean(x)) - (mean(Y)-mean(X)))

#Population, Samples, and Estimates Exercises #9
#For the females, our sample estimates were closer to the population difference than with males. What is a possible explanation for this?
#The population variance of the females is smaller than that of the males; thus, the sample variable has less variability.
popsd(x)
popsd(y)

###Central Limit Theorem (CLT) Exercises

url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/mice_pheno.csv"
filename <- basename(url)
download(url, destfile=filename)
dat <- na.omit( read.csv(filename) )
#Central Limit Theorem Exercises #1
pnorm(1)-pnorm(-1)

#Central Limit Theorem Exercises #2
pnorm(2)-pnorm(-2)

#Central Limit Theorem Exercises #3
pnorm(3)-pnorm(-3)

#Central Limit Theorem Exercises #4
y <- dat %>% filter(Sex=="M" & Diet=="chow") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
mean( abs(z) <=1 )

#Central Limit Theorem Exercises #5
mean( abs(z) <=2)

#Central Limit Theorem Exercises #6
mean( abs(z) <=3)

#Central Limit Theorem Exercises #7
qqnorm(z)
abline(0,1)

#Central Limit Theorem Exercises #8
mypar(2,2)
y <- filter(dat, Sex=="M" & Diet=="chow") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="F" & Diet=="chow") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="M" & Diet=="hf") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)
y <- filter(dat, Sex=="F" & Diet=="hf") %>% select(Bodyweight) %>% unlist
z <- ( y - mean(y) ) / popsd(y)
qqnorm(z);abline(0,1)

#Central Limit Theorem Exercises #9
y <- dat %>% filter(Sex=="M" & Diet=="chow") %>% select(Bodyweight) %>% unlist
set.seed(1)
avgs <- replicate(10000, mean( sample(y, 25)))
hist(avgs)
qqnorm(avgs)
qqline(avgs)
mean(avgs)
#Central Limit Theorem Exercises #10
popsd(avgs)

###CLT in Practice

dir <- 'https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/'
filename <- 'femaleControlsPopulation.csv'
#this is entire population of control mice
population <- read.csv(paste0(dir,filename)) %>% unlist

n <- 10000
nulls <- vector('numeric', n)
for (i in 1:n){
  control <- sample(population, 12)
  treatment <- sample(population, 12)
  nulls[i] <- mean(treatment) - mean(control)
}

hist(nulls, main = 'Null distribution of random variable')

mypar()
qqnorm(nulls)
qqline(nulls)

###T-test

dat <- read.csv('femaleMiceWeights.csv')

control <- dat %>% filter(Diet == 'chow') %>% select(Bodyweight) %>% unlist
treatment <- dat %>% filter(Diet == 'hf') %>% select(Bodyweight) %>% unlist

N <- length(treatment)
obs <- mean(treatment) - mean(control)

#estimating the standard error
se <- sqrt(var(treatment) / N + var(control) / N)

#t statistics
tstat <- obs / se

#central limit theorem tells us that the null distribution 
#is approximated by  a normal distribution which allows us
#to calculate p value without having data of entire population

#pnorm function tells us what what proportion of normally distributed
#are lower than t statistics value
#in this way we get p value 
#we multiplied number by 2 because we need 2 tails
2 * (1 - pnorm(tstat))

#to check how good is the above approximation
nulls <- vector('numeric', n)
for (i in 1:n){
  control <- sample(population, N)
  treatment <- sample(population, N)
  nulls[i] <- (mean(treatment) - mean(control)) / se
}

mypar()
qqnorm(nulls)
abline(0,1)
qqline(nulls)

### T-test in practice

#calculating mean difference and standard error by t.test() function
ttest <- t.test(treatment, control)
#here p value is different since we used t distribution approximation
#t distribution has bigger tail hence smaller p value
ttest

#we can check with qq plot whether or not t distribution is off
qqnorm(control)
qqline(control)

###CLT and t-distribution in Practice Exercises

url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleMiceWeights.csv"
filename <- "femaleMiceWeights.csv"
if(!file.exists("femaleMiceWeights.csv")) download(url,destfile=filename)
dat <- read.csv(filename)

#1
set.seed(1)
n <- 100
p <- 1/6
x <- replicate(10000,{
  z <- sample(1:6,n,replace=TRUE)
  (mean(z==6) - p) / sqrt(p*(1-p)/n)
}) 
qqnorm(x)
abline(0,1) #confirm it's well approximated with normal distribution
mean(abs(x) > 2)

#2

ps <- c(0.5,0.5,0.01,0.01)
ns <- c(5,30,30,100)
mypar(4,2)
for(i in 1:4){
  p <- ps[i]
  sides <- 1/p
  n <- ns[i]
  zs <- replicate(10000,{
    x <- sample(1:sides,n,replace=TRUE)
    (mean(x==1) - p) / sqrt(p*(1-p)/n)
  }) 
  hist(zs,nclass=7)
  qqnorm(zs)
  abline(0,1)
}

#correct one is p=0.5 and n=30

#3
X <- dat %>% filter(Diet=="chow") %>% select(Bodyweight) %>% unlist
Y <- dat %>% filter(Diet=="hf") %>% select(Bodyweight) %>% unlist
mean(X)

#4

#X¯  follows a normal distribution with mean  μX  and standard deviation  σX12√  where  σX  is the population standard deviation.

#5

#mean is 0

#6

sd(X)

#7

2 * ( 1-pnorm(2/sd(X) * sqrt(12) ) )

#8

sqrt(sd(X) ^2/12 + sd(Y) ^2/12)

#9

(mean(Y) - mean(X)) / sqrt(sd(Y) ^2/12 + sd(X) ^2/12)

#10

#If we apply the CLT, what is the distribution of this t-statistic?
  #Normal with mean 0 and standard deviation 1.

#11

Z <- ( mean(Y) - mean(X) ) / sqrt( var(X)/12 + var(Y)/12)
2*( 1-pnorm(Z)) 

#12
t.test(Y, X)
#13

#With the CLT distribution, we obtained a p-value smaller than 0.05 and with the t-distribution, one that is larger. They can't both be right. What best describes the difference?
#These are two different assumptions. The t-distribution accounts for the variability introduced by the estimation of the standard error and thus, under the null, large values are more probable under the null distribution.