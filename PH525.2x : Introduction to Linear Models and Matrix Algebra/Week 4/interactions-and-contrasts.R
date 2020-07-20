###Interactions and Contrasts I

library(downloader)
url <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/spider_wolff_gorb_2013.csv"
filename <- "spider_wolff_gorb_2013.csv"
if (!file.exists(filename)) download(url, filename)
spider <- read.csv("spider_wolff_gorb_2013.csv", skip=1)
head(spider)

#Recreating the boxplot from the article
boxplot(spider$friction ~ spider$type * spider$leg, 
        col=c("grey90","grey40"), las=2, 
        main="Comparison of friction coefficients of different leg pairs")

#Creating for demonstration a linear model with one variable
#Subsetting of spider data is necessary in order to have only one variable in this case L1:
spider.sub <- spider[spider$leg == "L1",]
fit <- lm(friction ~ type, data=spider.sub) #y will be friction and x- type
summary(fit)
(coefs <- coef(fit))

#Getting coefficients with a different method
s <- split(spider.sub$friction, spider.sub$type)
mean(s[["pull"]])
mean(s[["push"]]) - mean(s[["pull"]])

#Designing the matrix used in lm() function
X <- model.matrix(~ type, data=spider.sub)
colnames(X)
head(X)
tail(X)

#Making plot of the matrix by putting a black block for the 1's and a white block for the 0's:
library(rafalib)
imagemat(X, main="Model matrix for linear model with one variable")

#Examining the estimated coefficients:
set.seed(1) #same jitter in stripchart
stripchart(split(spider.sub$friction, spider.sub$type), 
           vertical=TRUE, pch=1, method="jitter", las=2, xlim=c(0,3), ylim=c(0,2))
a <- -0.25
lgth <- .1
library(RColorBrewer)
cols <- brewer.pal(3,"Dark2")
abline(h=0)
arrows(1+a,0,1+a,coefs[1],lwd=3,col=cols[1],length=lgth)
abline(h=coefs[1],col=cols[1])
arrows(2+a,coefs[1],2+a,coefs[1]+coefs[2],lwd=3,col=cols[2],length=lgth)
abline(h=coefs[1]+coefs[2],col=cols[2])
legend("right",names(coefs),fill=cols,cex=.75,bg="white")


###Interactions and Contrasts II

#A linear model with 2 variables

#Modeling matrix with type and leg:
X <- model.matrix(~ type + leg, data=spider)
colnames(X)
head(X)
imagemat(X, main="Model matrix for linear model with two factors")

#Fitting a linear model
fitTL <- lm(friction ~ type + leg, data=spider)
summary(fitTL)
(coefs <- coef(fitTL))

#Examining the estimated coefficients
spider$group <- factor(paste0(spider$leg, spider$type))
stripchart(split(spider$friction, spider$group), 
           vertical=TRUE, pch=1, method="jitter", las=2, xlim=c(0,11), ylim=c(0,2))
cols <- brewer.pal(5,"Dark2")
abline(h=0)
arrows(1+a,0,1+a,coefs[1],lwd=3,col=cols[1],length=lgth)
abline(h=coefs[1],col=cols[1])
arrows(3+a,coefs[1],3+a,coefs[1]+coefs[3],lwd=3,col=cols[3],length=lgth)
arrows(5+a,coefs[1],5+a,coefs[1]+coefs[4],lwd=3,col=cols[4],length=lgth)
arrows(7+a,coefs[1],7+a,coefs[1]+coefs[5],lwd=3,col=cols[5],length=lgth)
arrows(2+a,coefs[1],2+a,coefs[1]+coefs[2],lwd=3,col=cols[2],length=lgth)
segments(3+a,coefs[1]+coefs[3],4+a,coefs[1]+coefs[3],lwd=3,col=cols[3])
arrows(4+a,coefs[1]+coefs[3],4+a,coefs[1]+coefs[3]+coefs[2],lwd=3,col=cols[2],length=lgth)
segments(5+a,coefs[1]+coefs[4],6+a,coefs[1]+coefs[4],lwd=3,col=cols[4])
arrows(6+a,coefs[1]+coefs[4],6+a,coefs[1]+coefs[4]+coefs[2],lwd=3,col=cols[2],length=lgth)
segments(7+a,coefs[1]+coefs[5],8+a,coefs[1]+coefs[5],lwd=3,col=cols[5])
arrows(8+a,coefs[1]+coefs[5],8+a,coefs[1]+coefs[5]+coefs[2],lwd=3,col=cols[2],length=lgth)
legend("right",names(coefs),fill=cols,cex=.75,bg="white")
#The coefficient here are not anymore equal to mean().

#Contrasting/Comparing coefficients:
install.packages("contrast")
library(contrast) #Available from CRAN
L3vsL2 <- contrast(fitTL,list(leg="L3",type="pull"),list(leg="L2",type="pull"))
L3vsL2
coefs[4] - coefs[3]
(cT <- L3vsL2$X)
cT %*% coefs

#The reason changing type does not change the contrast is because it leads to addition of the `typepush` effect on both sides of the difference, which cancels out
L3vsL2.equiv <- contrast(fitTL,list(leg="L3",type="push"),list(leg="L2",type="push"))
L3vsL2.equiv$X


###Contrasts Exercises

species <- factor(c("A","A","B","B"))
condition <- factor(c("control","treated","control","treated"))
model.matrix(~ species + condition)

#1
y <- rnorm(4)
fit <- lm(y ~ species + condition)
contrast(fit, list(species="B",condition="control"), list(species="A",condition="treated"))$X

#2
L4vsL2 <- contrast(fitTL,list(leg="L4",type="pull"),list(leg="L2",type="pull"))
L4vsL2

#3
X <- model.matrix(~ type + leg, data=spider)
(Sigma <- sum(fitTL$residuals^2)/(nrow(X) - ncol(X)) * solve(t(X) %*% X))
C <- matrix(c(0,0,-1,0,1),1,5)

Sigma[3,5]


###Interactions and Contrasts III

X <- model.matrix(~ type + leg + type:leg, data=spider) # ':' means interaction
X <- model.matrix(~ type * leg, data=spider) #equivalent of the above formula
colnames(X)
head(X)
imagemat(X, main="Model matrix for linear model with interactions")

fitX <- lm(friction ~ type + leg + type:leg, data=spider)
summary(fitX)
coefs <- coef(fitX)

#Examining the estimated coefficients:
stripchart(split(spider$friction, spider$group), 
           vertical=TRUE, pch=1, method="jitter", las=2, xlim=c(0,11), ylim=c(0,2))
cols <- brewer.pal(8,"Dark2")
abline(h=0)
arrows(1+a,0,1+a,coefs[1],lwd=3,col=cols[1],length=lgth)
abline(h=coefs[1],col=cols[1])
arrows(2+a,coefs[1],2+a,coefs[1]+coefs[2],lwd=3,col=cols[2],length=lgth)
arrows(3+a,coefs[1],3+a,coefs[1]+coefs[3],lwd=3,col=cols[3],length=lgth)
arrows(5+a,coefs[1],5+a,coefs[1]+coefs[4],lwd=3,col=cols[4],length=lgth)
arrows(7+a,coefs[1],7+a,coefs[1]+coefs[5],lwd=3,col=cols[5],length=lgth)
#now the interactions:
segments(3+a,coefs[1]+coefs[3],4+a,coefs[1]+coefs[3],lwd=3,col=cols[3])
arrows(4+a,coefs[1]+coefs[3],4+a,coefs[1]+coefs[3]+coefs[2],lwd=3,col=cols[2],length=lgth)
arrows(4+a,coefs[1]+coefs[2]+coefs[3],4+a,coefs[1]+coefs[2]+coefs[3]+coefs[6],lwd=3,col=cols[6],length=lgth)
segments(5+a,coefs[1]+coefs[4],6+a,coefs[1]+coefs[4],lwd=3,col=cols[4])
arrows(6+a,coefs[1]+coefs[4],6+a,coefs[1]+coefs[4]+coefs[2],lwd=3,col=cols[2],length=lgth)
arrows(6+a,coefs[1]+coefs[4]+coefs[2],6+a,coefs[1]+coefs[4]+coefs[2]+coefs[7],lwd=3,col=cols[7],length=lgth)
segments(7+a,coefs[1]+coefs[5],8+a,coefs[1]+coefs[5],lwd=3,col=cols[5])
arrows(8+a,coefs[1]+coefs[5],8+a,coefs[1]+coefs[5]+coefs[2],lwd=3,col=cols[2],length=lgth)
arrows(8+a,coefs[1]+coefs[5]+coefs[2],8+a,coefs[1]+coefs[5]+coefs[2]+coefs[8],lwd=3,col=cols[8],length=lgth)
legend("right",names(coefs),fill=cols,cex=.75,bg="white")


###Interactions and Contrasts IV


L2push.vs.pull <- contrast(fitX,
                           list(leg="L2", type = "push"), 
                           list(leg="L2", type = "pull"))
L2push.vs.pull
coefs[2] + coefs[6]

#More complicated contrast. Difference of differences:
library(multcomp) ##Available from CRAN
C <- matrix(c(0,0,0,0,0,-1,1,0), 1)
L3vsL2interaction <- glht(fitX, linfct=C)
summary(L3vsL2interaction)
coefs[7] - coefs[6]

# Analysis of Variance:
anova(fitX)


###Interactions Exercises

spider$log2friction <- log2(spider$friction)
boxplot(log2friction ~ type*leg, data=spider)

#1
fit <- lm(log2friction ~ type * leg, data=spider)
summary(fit)

#2
anova(fit)

#3
contrast(fit, list(type="pull",leg="L2"), list(type="pull",leg="L1"))
coef(fit)["legL2"]

#4
contrast(fit, list(type="push",leg="L2"), list(type="push",leg="L1"))
coef(fit)["legL2"] + coef(fit)["typepush:legL2"]

#5
N <- 40
p <- 4
group <- factor(rep(1:p,each=N/p))
X <- model.matrix(~ group)
Y <- rnorm(N,mean=42,7)
mu0 <- mean(Y)
initial.ss <- sum((Y - mu0)^2)
s <- split(Y, group)
after.group.ss <- sum(sapply(s, function(x) sum((x - mean(x))^2)))
(group.ss <- initial.ss - after.group.ss)
group.ms <- group.ss / (p - 1)
after.group.ms <- after.group.ss / (N - p)
f.value <- group.ms / after.group.ms

set.seed(1)
Fs <- replicate(1000, {
  Y = rnorm(N,mean=42,7)
  mu0 = mean(Y)
  initial.ss = sum((Y - mu0)^2)
  s = split(Y, group)
  after.group.ss = sum(sapply(s, function(x) sum((x - mean(x))^2)))
  (group.ss = initial.ss - after.group.ss)
  group.ms = group.ss / (p - 1)
  after.group.ms = after.group.ss / (N - p)
  f.value = group.ms / after.group.ms
  return(f.value)
})
mean(Fs)

hist(Fs, col="grey", border="white", breaks=50, freq=FALSE)
xs <- seq(from=0,to=6,length=100)
lines(xs, df(xs, df1 = p - 1, df2 = N - p), col="red")


###Interactions and Contrasts V

spider$group <- factor(paste0(spider$leg, spider$type))
X <- model.matrix(~ 0 + group, data=spider)
colnames(X)
head(X)
imagemat(X, main="Model matrix for linear model with group variable")

#Fitting linear model
fitG <- lm(friction ~ 0 + group, data=spider)
summary(fitG)
coefs <- coef(fitG)

#Examining the estimated coefficients:
stripchart(split(spider$friction, spider$group), 
           vertical=TRUE, pch=1, method="jitter", las=2, xlim=c(0,11), ylim=c(0,2))
cols <- brewer.pal(8,"Dark2")
abline(h=0)
for (i in 1:8) {
  arrows(i+a,0,i+a,coefs[i],lwd=3,col=cols[i],length=lgth)
}
legend("right",names(coefs),fill=cols,cex=.75,bg="white")

#Differences of differences when there is no intercept
C <- matrix(c(0,0,1,-1,-1,1,0,0), 1)
groupL3vsL2interaction <- glht(fitG, linfct=C)
summary(groupL3vsL2interaction)
names(coefs)
(coefs[6] - coefs[5]) - (coefs[4] - coefs[3])