###Expressing Design Formula Exercises

#day: A B C
#condition: --------------
#treated | 2 2 2
#control | 2 2 2

#1
#Given the factors we have defined above, and not defining any new ones, which of the following R formula will produce a design matrix (model matrix) that let's us analyze the effect of condition, controlling for the different days:
  #~ day + condition

#Explanation
#Using the ~ and the names for the two variables we want in the model will produce a design matrix controlling for all levels of day and all levels of condition, so ~ day + condition. We do not use the levels A,B,C etc in the design formula.


###Linear Models in Practice I

url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleMiceWeights.csv"
filename <- "femaleMiceWeights.csv"
library(downloader)
if (!file.exists(filename)) download(url, filename)
set.seed(1)

dat <- read.csv("femaleMiceWeights.csv") ##previously downloaded
stripchart(dat$Bodyweight ~ dat$Diet, vertical=TRUE, method="jitter",
           main="Bodyweight over Diet")

#A linear model with one variable
levels(dat$Diet)
#model.matrix() formula chooses non reference level of factor to base its matrix variables
X <- model.matrix(~ Diet, data=dat)
head(X)
#if you want to change levels of factor for model.matrix, you can use relevel() function and set reference desired level of factor


###Linear Models in Practice II

#Running the linear model

fit <- lm(Bodyweight ~ Diet, data=dat)
summary(fit)
coefs <- coef(fit)

#To prove the lm() function is doing the same
Y <- dat$Bodyweight
X <- model.matrix(~ Diet, data=dat)
solve(t(X) %*% X) %*% t(X) %*% Y

#We can get coefficient without lm() function
s <- split(dat$Bodyweight, dat$Diet)
mean(s[["chow"]])
mean(s[["hf"]]) - mean(s[["chow"]])

#T-test would give us the same results
ttest <- t.test(s[["hf"]], s[["chow"]], var.equal=TRUE)
summary(fit)$coefficients[2,3]
ttest$statistic


###Linear Models in Practice Exercises

#1
nx <- 5
ny <- 7
X <- cbind(rep(1,nx + ny),rep(c(0,1),c(nx, ny)))
(t(X) %*% X)[1,1]

#2
t(X) %*% X