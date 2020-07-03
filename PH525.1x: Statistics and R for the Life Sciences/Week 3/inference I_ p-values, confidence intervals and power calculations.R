###T-test Exercises

library(downloader)
library(dplyr)

url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/babies.txt"
filename <- basename(url)
download(url, destfile=filename)
babies <- read.table("babies.txt", header=TRUE)
head(babies)

bwt.nonsmoke <- babies %>% filter(smoke==0) %>% select(bwt) %>% unlist 
bwt.smoke <- babies %>% filter(smoke==1) %>% select(bwt) %>% unlist

library(rafalib)
mean(bwt.nonsmoke)-mean(bwt.smoke)
popsd(bwt.nonsmoke)
popsd(bwt.smoke)

#1
set.seed(1)
N <- 25
dat.ns <- sample(bwt.nonsmoke, N)
dat.s <- sample(bwt.smoke, N)
obs  <- mean(dat.ns) - mean(dat.s)
se <- sqrt(var(dat.ns) / N + var(dat.s) / N)
tval <- abs(obs / se)

#2 

pval <- 1-(pnorm(abs(tval))-pnorm(-abs(tval)))

#3

#Because of the symmetry of the standard normal distribution, there is a simpler way to calculate the probability that a t-value under the null could have a larger absolute value than tval. Choose the simplified calculation from the following:
  #2*pnorm(-abs(tval))

###Confidence Intervals

set.seed(1)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename <- basename(url)
download(url, destfile=filename)
chowPopulation <- unlist( read.csv(filename) ) %>% unlist

mu_chow <- mean(chowPopulation)
print(mu_chow)

N <- 30

chow <- sample(chowPopulation, N)
mean(chow)

#we assume the sample average is normally distributed

se <- sd(chow) / sqrt(N)
se 

#t statistics is 95% of the time between -Q and Q because random variable 0-1 is between -Q and Q 95% of the time
Q <- qnorm(1- 0.05 / 2)
(mean(chow) - mean(mu_chow)) / se

# if we clear out this equation so that the truth is in the middle 
#we get the following equation
#this interval will be on top of the mean of population 95 % of the time. 
interval <- c(mean(chow)-Q*se, mean(chow)+Q*se)

#we can check if it's true of false
interval[1] < mu_chow & interval[2] > mu_chow

B <- 250
mypar()
plot(mean(chowPopulation) + c(-7,7), c(1,1), type = 'n', xlab = 'weigth', ylab = 'interval', ylim = c(1, B))
abline(v=mean(chowPopulation))
for(i in 1:B){
  chow <- sample(chowPopulation, N)
  se <- sd(chow) / sqrt(N)
  interval <- c(mean(chow)-Q*se, mean(chow)+Q*se)
  covered <- mean(chowPopulation) <= interval[2]
  color <- ifelse(covered, 1, 2)
  lines(interval, c(i,i), col= color)
}
#green ones show how 95% of the time interval fall on top of mean of population and red ones show how they do not fall 5% of the time

#we will check if ctl works when sample size is much smaller

N <- 5 
B <- 250
mypar()
plot(mean(chowPopulation) + c(-7,7), c(1,1), type = 'n', xlab = 'weigth', ylab = 'interval', ylim = c(1, B))
abline(v=mean(chowPopulation))
for(i in 1:B){
  chow <- sample(chowPopulation, N)
  se <- sd(chow) / sqrt(N)
  interval <- c(mean(chow)-Q*se, mean(chow)+Q*se)
  covered <- mean(chowPopulation) <= interval[2]
  color <- ifelse(covered, 1, 2)
  lines(interval, c(i,i), col= color)
}
#we can see percentage of red lines is higher than 5% because we did not make intervals big enough
#Q that we used previously is not right anymore since we found it based on ctl


# now check what happens when we find Q with t distribution of approximation

B <- 250
mypar()
plot(mean(chowPopulation) + c(-7,7), c(1,1), type = 'n', xlab = 'weigth', ylab = 'interval', ylim = c(1, B))
abline(v=mean(chowPopulation))
N <- 5
Q <- qt(1-0.05/2, df=4)
for(i in 1:B){
  chow <- sample(chowPopulation, N)
  se <- sd(chow) / sqrt(N)
  interval <- c(mean(chow)-Q*se, mean(chow)+Q*se)
  covered <- mean(chowPopulation) <= interval[2]
  color <- ifelse(covered, 1, 2)
  lines(interval, c(i,i), col= color)
}
#as we increased size of Q, in case of smaller sample we could still get 5% red this time
#when you cannot construct confidence interval based on ctl then you can use t distribution

###Confidence Intervals Exercises

#1

set.seed(1)
N <- 25
dat.ns <- sample(bwt.nonsmoke, N)
dat.s <- sample(bwt.smoke, N)
se <- sqrt( sd( dat.ns)^2/N + sd( dat.s)^2/N )
qt(0.995,48)*se

#2

#Which of the following sentences about a Type I error is not true?
  #From the original data alone, you can tell whether you have made a Type I error.

#3

set.seed(1)
N <- 5
dat.ns <- sample(bwt.nonsmoke, N)
dat.s <- sample(bwt.smoke, N)
t.test(dat.ns, dat.s)$p.value
##Note that if you defone dat.s before dat.ns, you get a different answer
##due to sampling randomness
##tolerance is set to accept both

###Power Calculations

dat <- read.csv("mice_pheno.csv")
controlPopulation <- dat %>% filter(Sex == "F" & Diet == "chow") %>%  
  select(Bodyweight) %>% unlist
hfPopulation <- dat %>% filter(Sex == "F" & Diet == "hf") %>%  
  select(Bodyweight) %>% unlist
mu_hf <- mean(hfPopulation)
mu_control <- mean(controlPopulation)
#we check if we should reject or not the null hypothesis
print(mu_hf - mu_control)
#difference exists hence we can should perform test to reject it
print((mu_hf - mu_control)/mu_control * 100) #9% is important scientifically

#first try to see if we can reject the null hypothesis
set.seed(1)
N <- 5
hf <- sample(hfPopulation,N)
control <- sample(controlPopulation,N)
t.test(hf,control)$p.value #there's not enough evidence to reject it. We have Type II error. We have difference but cannot reject it. This mistake happens because it's a random data.

N <- 12
#Alpha is the character we usually use to name the cutoff at which we reject the hypothesis if the p-value is smaller than .05
alpha <- 0.05
B <- 2000

#function to check whether or not to reject the null hypothesis

reject <- function(N, alpha = 0.05){
  hf <- sample(hfPopulation,N) 
  control <- sample(controlPopulation,N)
  pval <- t.test(hf,control)$p.value
  pval < alpha
}

#if it return TRUE, we reject, if FALSE- we reject.
reject(12)

#see results 2000 times
rejections <- replicate(B,reject(N))

mean(rejections) * 100 #shows we rejected 21% the null hypothesis

#So this is power. It's the probability of rejecting the hypothesis when the alternative is true.

#Let's see how power improves with N. We will use the function `sapply`, which applies a function to each of the elements of a vector. We want to repeat the above for the following sample size

Ns <- seq(5, 50, 5)
power <- sapply(Ns,function(N){
  rejections <- replicate(B, reject(N))
  mean(rejections)
})

#make plot for demonstration
plot(Ns, power, type="b") #For each of the three simulations, the above code returns the proportion of times we reject. Not surprisingly power increases with N

###Power Calculations Exercises

#1

#The p-value is larger than 0.05 so using the typical cut-off, we would not reject. This is a type II error. Which of the following is *not* a way to decrease this type of error?
  #Find a population for which the null is not true.

#2
set.seed(1)
N <- 5
rejects <- replicate(10000,{
  dat.ns <- sample(bwt.nonsmoke , N)
  dat.s <- sample(bwt.smoke , N)
  t.test(dat.s, dat.ns)$p.value < 0.05
})
mean(rejects)

#3

Ns <- seq(30, 120, 30)
power <- sapply(Ns,function(N){
  rejections <- replicate(B, reject(N))
  mean(rejections)
})
Ns[ which.min( abs( power - .8) ) ] 

#4
Ns=c(30,60,90,120)
res <- sapply(Ns, function(N){
  set.seed(1)
  rejects <- replicate(10000,{
    dat.ns <- sample(bwt.nonsmoke , N)
    dat.s <- sample(bwt.smoke , N)
    t.test(dat.s, dat.ns)$p.value < 0.01
  })
  mean(rejects)
})
Ns[ which.min( abs( res - .8) ) ] 