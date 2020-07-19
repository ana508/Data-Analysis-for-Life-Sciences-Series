###Inference Review Exercises

g <- 9.8 ## meters per second
h0 <- 56.67
v0 <- 0
n <- 25
tt <- seq(0,3.4,len=n) ##time in secs, t is a base function
y <- h0 + v0 *tt  - 0.5* g*tt^2 + rnorm(n,sd=1)

#To obtain the LSE in R we could write:
  
X <- cbind(1,tt,tt^2)
A <- solve(crossprod(X))%*%t(X)

#1
#Given how we have defined A, which of the following is the LSE of g, the acceleration due to gravity?
  #-2 * (A %*% y) [3]

#2
set.seed(1)
betahat <- replicate(100000,{
  y = 56.67 - 0.5*g*tt^2 + rnorm(n,sd=1)
  betahats = -2*A%*%y
  return(betahats[3])
})
sd(betahat)

#3
library(UsingR)
x <- father.son$fheight
y <- father.son$sheight
n <- length(y)
N <-  50
set.seed(1)
betahat <- replicate(10000,{
  index = sample(n,N)
  sampledat = father.son[index,]
  x = sampledat$fheight
  y = sampledat$sheight
  lm(y~x)$coef[2]
})
sd(betahat)

#4
#Which of the following is closest to the covariance between father heights and son heights?
  #4
Y=father.son$fheight
X=father.son$sheight
mean( (Y - mean(Y))*(X-mean(X) ) )
#Note we know it can't be 0 because only independent variables have covariance 0. And we know it can't be negative because these variables are positively correlated thus (Y - mean(Y)) and (X-mean(X) ) tend to have the same sign


###Technical note on variance

#The standard approach to writing linear models either assume the X are fixed or that we are 
#conditioning on them. Thus  Xβ  has no variance as the X is considered fixed. This is why we 
#write  Var(Yi)=Var(ei)=σ2 . This can cause confusion in practice because if you, for example, compute the following:
x <-  father.son$fheight
beta <-  c(34,0.5)
var(beta[1]+beta[2]*x)

#it is nowhere near 0. This is an example in which we have to be careful in distinguishing code from math. The function 
#var is simply computing the variance of the list we feed it, while the mathematical use of var is considering only quantities 
#that are random variables. In the R code above, x is not fixed at all: we are letting it vary but when we write  Var(Yi)=σ2  we 
#are imposing, mathematically, x to be fixed. Similarly if we use R to compute the variance of Y in our object dropping example we 
#obtain something very different than  σ2=1  (the known variance):
y <- h0 + v0*tt  - 0.5*g*tt^2 + rnorm(n,sd=1)
var(y)

#Again, this is because we are not fixing tt.

###Standard Errors Exercises

x = father.son$fheight
y = father.son$sheight
n = length(y)
N = 50
set.seed(1)
index = sample(n,N)
sampledat = father.son[index,]
x = sampledat$fheight
y = sampledat$sheight
betahat = lm(y~x)$coef

#1
fit = lm(y ~ x)
fit$fitted.values
SSR <- sum((y - fit$fitted.values)^2)
sigma2 <- SSR / 48

#2
X <- cbind(1, x)
solve(t(X) %*% X)[1,1]

#3
#Now we are one step away from the standard error of beta-hat. 
#Take the diagonals from the  (XTX)−1  matrix above, using the diag() function. 
#Now multiply our estimate of  σ2  and the diagonals of this matrix. This is the 
#estimated variance of beta-hat, so take the square root of this. You should end up 
#with two numbers, the standard error for the intercept and the standard error for the slope.
(sqrt(diag(solve(t(X) %*% X)) * sigma2))[2]