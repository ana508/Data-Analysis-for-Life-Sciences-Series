install.packages('UsingR')

#selecting father.son data set from UsingR library
library(UsingR)
x <- father.son$fheight

#how many individuals
length(x)

#round randomly chosen heigths to the nearest 10th of an inch
round(sample(x, 20), 1)

#creating a histogram
hist(x, breaks = seq(floor(min(x)), ceiling(max(x))), main = 'Height histogram',xlab = 'Height in inches')

#the empirical commutative distribution function
xs <- seq(floor(min(x)), ceiling(max(x)), 0.1)
plot(xs, ecdf(x)(xs), type = 'l', main = 'Height plot', xlab = 'Height in inches', ylab = 'F(x)')


#normal approximation
mean(x> 70)
#we use 1 since we look at right part of distribution
1 - pnorm(70, mean(x), sd(x))
#if we look at left part of distribution
mean(x < 70)
pnorm(70, mean(x), sd(x))

#Quantile-Quantile plot, or in short, Q-Q plot
ps <- seq(0.01, 0.99, 0.1)
#percentile
qs <- quantile(x, ps)
#percentile for normal distribution
normalqs <- qnorm(ps, mean(x), sd(x))
plot(normalqs, qs, xlab = 'Normal percentiles', ylab = 'Height')
abline(0, 1) #identity line

#easier function
qqnorm(x)
qqline(x)

#QQ-plot Exercises #1 and #2
load("skew.RData")
dim(dat)

par(mfrow=c(3,3))
for (i in 1:9) {
  x <- dat[,i]
  qqnorm(x, main = paste(i, "-th column", ssep = ''))
  qqline(x)
}

#examining positive(4) and negative(9) skews in histograms
hist(dat[,4])
hist(dat[,9])


#histogram of not normally distributed data
hist(exec.pay)
qqnorm(exec.pay)
qqline(exec.pay)

#creating a boxplot
boxplot(exec.pay, ylab = '10,000s of dollars', ylim = c(0, 400))
