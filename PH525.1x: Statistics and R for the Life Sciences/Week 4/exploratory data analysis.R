library(rafalib)
library(dplyr)
library(downloader)

###Scatter plot

set.seed(1)
data(father.son,package="UsingR")
x <- father.son$fheight
y <- father.son$sheight
mypar()
plot(x, y, xlab = "Father's height in inches", ylab = "Son's height in inches", main = paste("correlation = ",signif(cor(x, y), 2)))
#if we summarized data here with mean and standard deviation we would not be able to see how father's and son's heights both go up. This particular characteristic is best summarized by correlation.
#To calculate correlation we standardize the units and then multiply them together and take the mean of that. That's going to give us 0.5.

boxplot(split(y, round(x)))
print(mean(y[round(x) == 72]))

#Doing the same thin but by standardizing the data
x <- (x-mean(x)) / sd(x)
y <- (y-mean(y)) / sd(y)
means <- tapply(y, round(x*4)/ 4, mean)
fatherheights <- as.numeric(names(means))
plot(fatherheights, means, ylab = "average of strata of son's heights", ylim =range(fatherheights))
abline(0, cor(x, y))

#for data that is normally distributed and 2D, 2 means, 2 standard deviations and correlation are the best way to characterize given data.

#To see what happens when data is uncorrelated
set.seed(1)
a <- rnorm(100); a[1] <- 25
b <- rnorm(100); b[1] <- 26
plot(a, b, main = paste("correlation = ", signif(cor(a, b), 2)))


###Scatter plot Exercises

data(nym.2002, package="UsingR")

#1
males <- nym.2002 %>% filter(gender == "Male")
females <- nym.2002 %>% filter(gender == "Female")
cor(males$age,males$time)

#2
cor(females$age, females$time)

#3
mypar(2,2)
plot(females$age, females$time)
plot(males$age, males$time)
group <- floor(females$age/5) * 5
boxplot(females$time~group)
group <- floor(males$age/5) * 5
boxplot(males$time~group)
#After examining the data, what is a more reasonable conclusion?
  #Finish times are constant up through around 50-60, then they get slower.

###Symmetry of Log Ratios Exercises

time <- sort(nym.2002$time)

#1
min(time) / median(time)

#2
max(time) / median(time)


###Plots to Avoid Exercises

#1
#When is it appropriate to use pie charts or donut charts?
  #Never.

#2
#The use of pseudo-3D plots in the literature mostly adds:
  #Confusion.