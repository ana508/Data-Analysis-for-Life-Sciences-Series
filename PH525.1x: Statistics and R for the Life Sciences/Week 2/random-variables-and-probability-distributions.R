library(dplyr)
dat <- read.csv("femaleMiceWeights.csv")

control <- dat %>% filter(Diet == 'chow') %>% select(Bodyweight) %>% unlist
treatment <- dat %>% filter(Diet == 'hf') %>% select(Bodyweight) %>% unlist

#while comparing results of those two values we won't
#always get higher weight in treatment since article said
#that mice with high diet weighted more
#this conclusion was made by finding average of those 2 values
#then after comparing them each other it is obvious that 'hf' Diet made
#mice weigh more

#they are in statistics called random variable
obs <- mean(treatment) - mean(control)

library(downloader)
dir <- 'https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/'
filename <- 'femaleControlsPopulation.csv'
#this is entire population of control mice
population <- read.csv(paste0(dir,filename)) %>% unlist

#results will demonstrate that those 3 give different numbers
mean(sample(population, 12))
mean(control)
mean(treatment)

#Random Variables Exercises #1
mean(population)

#Random Variables Exercises #2
set.seed(1)
abs(mean(sample(population,5)) - mean(population))

#Random Variables Exercises #3
set.seed(5)
abs(mean(sample(population,5)) - mean(population))

#Random Variables Exercises #4

#Why are the answers from 2 and 3 different?
#Because the average of the samples is a random variable.

###null hypothesis
#there is no connection between this 2 values
#because both of them are chosen randomly for no specific characteristic
#difference between averages will be always various
control <- sample(population, 12)
treatment <- sample(population, 12)
mean(control) - mean(treatment)

#finding differences 10000 times and recording all of them
n <- 10000
nulls <- vector('numeric', n)
for (i in 1:n){
  control <- sample(population, 12)
  treatment <- sample(population, 12)
  nulls[i] <- mean(treatment) - mean(control)
}
#the biggest value
max(nulls)

#building histogram
hist(nulls, main = 'Null distribution of random variable')

### p-value
#total number of times null distribution is bigger than observation of an experiment divided by n
sum(nulls > obs) / n

#simpler way of above operation
mean(nulls > obs)

#Null Distributions Exercises #1
set.seed(1)
n <- 1000
averages <- vector('numeric', n)
for (i in 1:n){
  averages[i] <- mean(sample(population, 5))
}
sum(abs(averages - mean(population)) > 1) / n

#Null Distributions Exercises #2
set.seed(1)
n <- 10000
averages <- vector('numeric', n)
for (i in 1:n){
  averages[i] <- mean(sample(population, 5))
}
sum(abs(averages - mean(population)) > 1) / n

#Null Distributions Exercises #3
set.seed(1)
n <- 1000
averages <- vector('numeric', n)
for (i in 1:n){
  averages[i] <- mean(sample(population, 50))
}
sum(abs(averages - mean(population)) > 1) / n


install.packages("gapminder")
library(gapminder)
data(gapminder)
head(gapminder)

#Probability Distributions Exercises #1
mean(gapminder$lifeExp[gapminder$year==1952] <= 40)

#Probability Distributions Exercises #2
mean(gapminder$lifeExp[gapminder$year==1952] <= 60) - mean(gapminder$lifeExp[gapminder$year==1952] <= 40)

#Probability Distributions Exercises #3

#plot the proportions of countries with life expectancy q for a range of different years
prop <- function(q) {
  mean(x <= q)
}
x <- gapminder$lifeExp[gapminder$year==1952]
prop(40)

#Now let's build a range of qs that we can apply the function to:
qs = seq(from=min(x), to=max(x), length=20)

#Use sapply() to apply the prop function to each element of qs:
props <- sapply(qs, prop)

plot(qs, props, main = 'homemade empirical cdf')

#We could also have written this in one line, by defining the prop function inside of sapply() but without naming it:
props <- sapply(qs, function(q) mean(x <= q))

#Plotting with the pre-built function
plot(ecdf(x), main='Pre-built empirical cdf')
