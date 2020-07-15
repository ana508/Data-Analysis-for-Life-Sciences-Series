###Matrix Notation Exercises

#1
X = matrix(1:1000,100,10)
X[25,3]

#2
m <-1:10 
x <- cbind(x=m,x1=2*m, x2=3*m, x3=4*m, x4=5*m)
sum(x[7,])

#3
#Which of the following creates a matrix with multiples of 3 in the third column?
  #matrix(1:60,20,3,byrow=TRUE)

###Matrix Operation Exercises

#1
#Which of the following is NOT equivalent to X?
  #X %*% matrix(1,ncol(X) )
#The first transposes the transpose, so we end up with our original X. The third is multiplying each element by 1, and the fourth is multiplying X by the identity. The second is not even guaranteed to have the same dimensions as X.

#2
X = matrix(c(3,2,1,5,4,2,-1,0,-5,2,5,0,1,-1,-5,1),4,4)
y = c(10,5,7,4)
sol = solve(X,y)
sol[ 3 ]

#3
a <- matrix(1:12, nrow=4)
b <- matrix(1:15, nrow=3)
n <- a %*%b
n[3,2]

#4
sum(a[3,] * b[,2])