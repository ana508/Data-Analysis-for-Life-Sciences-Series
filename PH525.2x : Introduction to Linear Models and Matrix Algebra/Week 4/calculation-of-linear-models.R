###Collinearity

Sex <- c(0,0,0,0,1,1,1,1)
A <-   c(1,1,0,0,0,0,0,0)
B <-   c(0,0,1,1,0,0,0,0)
C <-   c(0,0,0,0,1,1,0,0)
D <-   c(0,0,0,0,0,0,1,1)
X <- model.matrix(~Sex+A+B+C+D-1)
cat("ncol=",ncol(X),"rank=", qr(X)$rank,"\n")

#Rank is the number of columns that are independent from others.
#If you design your experiment correctly, every column should be independent of others.

#The above code needs to be balanced, brecasue of that we will not be able to estimate the effect of sex:
Sex <- c(0,1,0,1,0,1,0,1)
A <-   c(1,1,0,0,0,0,0,0)
B <-   c(0,0,1,1,0,0,0,0)
C <-   c(0,0,0,0,1,1,0,0)
D <-   c(0,0,0,0,0,0,1,1)
X <- model.matrix(~Sex+A+B+C+D-1)
cat("ncol=",ncol(X),"rank=", qr(X)$rank,"\n")


###Collinearity Exercises

#1
m <- matrix(c(1,1,1,1,0,0,1,1,0,1,0,1,0,0,0,1),4,4)
qr(m)$rank
#Which of the above design matrices does NOT have the problem of collinearity?
  #E

#2
sex <- factor(rep(c("female","male"),each=4))
trt <- factor(c("A","A","B","B","C","C","D","D"))
X <- model.matrix( ~ sex + trt)
qr(X)$rank
Y <- 1:8
makeYstar <- function(a,b) Y - X[,2] * a - X[,5] * b

fitTheRest <- function(a,b) {
  Ystar <- makeYstar(a,b)
  Xrest <- X[,-c(2,5)]
  betarest <- solve(t(Xrest) %*% Xrest) %*% t(Xrest) %*% Ystar
  residuals <- Ystar - Xrest %*% betarest
  sum(residuals^2)
}

fitTheRest(1,2)

#3
expand.grid(1:3,1:3)
betas <- expand.grid(-2:8,-2:8)
rss <- apply(betas,1,function(x) fitTheRest(x[1],x[2]))

themin <- min(rss)
betas[which(rss==themin),]
#Which of the following pairs of values minimizes the RSS?
  #All of the above. There is no single minimum.

library(rafalib)
## plot the pairs what are minimum
plot(betas[which(rss==themin),])


###QR Factorization

#Demonstration of numerically unstable solve() function in some circumstances:
n <- 50;M <- 500
x <- seq(1,M,len=n)
X <- cbind(1,x,x^2,x^3)
colnames(X) <- c("Intercept","x","x2","x3")
beta <- matrix(c(1,1,1,1),4,1)
set.seed(1)
y <- X%*%beta+rnorm(n,sd=1)
solve(crossprod(X)) %*% crossprod(X,y)
#We get an error. To see the reason behind that, we see the actual cross product:
options(digits=4)
log10(crossprod(X))
#The problem is that values of columns are way too different.

#Finding LSE with QR:
QR <- qr(X)
Q <- qr.Q( QR )
R <- qr.R( QR )
(betahat <- backsolve(R, crossprod(Q,y) ) )

(betahat <- solve.qr(QR, y)) #built-in fuction for doing the above code faster

#Fitter values:
mypar(1,1)
plot(x,y)
fitted <- tcrossprod(Q)%*%y
lines(x,fitted,col=2)

#Standard errors:
df <- length(y) - QR$rank
sigma2 <- sum((y-fitted)^2)/df
varbeta <- sigma2*chol2inv(qr.R(QR))
SE <- sqrt(diag(varbeta))
cbind(betahat,SE)

#We can get the same result:
summary(lm(y~0+X))$coef


###QR Exercises

url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/spider_wolff_gorb_2013.csv"
filename <- "spider_wolff_gorb_2013.csv"
library(downloader)
if (!file.exists(filename)) download(url, filename)
spider <- read.csv(filename, skip=1)
fit <- lm(friction ~ type + leg, data=spider)
betahat <- coef(fit)
Y <- matrix(spider$friction, ncol=1)
X <- model.matrix(~ type + leg, data=spider)

#1
Q <- qr.Q(qr(X))
Q[1]

#2
R <- qr.R(qr(X))
R[1]

#3
head(crossprod(Q, Y))[1]