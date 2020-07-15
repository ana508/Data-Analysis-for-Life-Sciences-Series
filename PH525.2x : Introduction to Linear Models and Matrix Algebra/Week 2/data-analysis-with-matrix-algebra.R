###Examples

#Finding of average
data(father.son,package="UsingR")
y <- father.son$sheight
print(mean(y))
N <- length(y)
Y<- matrix(y,N,1)
A <- matrix(1,N,1)
barY <- t(A)%*%Y / N
print(barY)

#Finding of sample variance
barY <- crossprod(A,Y) / N
print(barY)
r <- y - barY
crossprod(r)/N
#This above code is equivalent to 
var(y) * (N-1) / N

#Minimizing RSS
x <- father.son$fheight
y <- father.son$sheight
X <- cbind(1,x)
betahat <- solve( t(X) %*% X ) %*% t(X) %*% y
###or
betahat <- solve( crossprod(X) ) %*% crossprod( X, y )

###Matrix Algebra in Practice I

library(rafalib)
set.seed(1)
g <- 9.8 ##meters per second
n <- 25
tt <- seq(0,3.4,len=n) ##time in secs, note: we use tt because t is a base function
d <- 56.67  - 0.5*g*tt^2 + rnorm(n,sd=1) ##meters
mypar()
plot(tt,d,ylab="Distance in meters",xlab="Time in seconds")
lines(tt,d,col=2)

rss <- function(Beta0,Beta1,Beta2){
  r <- y - (Beta0 + Beta1*tt + Beta2*tt^2)
  sum(r^2)
}
Beta2s <- seq(-10,0,len=100)
RSS <- sapply(Beta2s, rss, Beta0=65,Beta1=0)
lines(Beta2s,RSS, type = 'l', col = 2)

###Matrix Algebra in Practice II

X <- cbind(1, tt, tt^2)
head(X)
Beta <- matrix(c(55,0,5), 3, 1)

#residual
y <- y[1:25]
r <- y - X %*% Beta
RSS <- t(r) %*% r
rss(55,0,5)

#same function for above code
RSS <- crossprod(r)
RSS

betahat <- solve(t(X) %*%X) %*% t(X) %*% y
#it's better to write this way:
betahat <- solve(crossprod(X)) %*% crossprod(X,y)

#QR decomposition
QR <- qr(X)
Q <- qr.Q(QR)
R <- qr.R(QR)
backsolve(R,crossprod(Q,y)) #more stable version than solve function

###Matrix Algebra Examples Exercises

X <- matrix(c(1,1,1,1,0,0,1,1),nrow=4)
rownames(X) <- c("a","a","b","b")
beta <- c(5, 2)

a <- X %*% beta

#1
a[1:2,]

#2
a[3:4,]


X <- matrix(c(1,1,1,1,1,1,0,0,1,1,0,0,0,0,0,0,1,1),nrow=6)
rownames(X) <- c("a","a","b","b","c","c")
beta <- c(10,3,-3)

#3
a <- X %*% beta
a[3:4,]

#4
a[5:6,]