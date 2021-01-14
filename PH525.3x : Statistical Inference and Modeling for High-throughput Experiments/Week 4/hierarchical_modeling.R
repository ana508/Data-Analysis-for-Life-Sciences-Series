#Bayes' Rule in Practice Exercises

tmpfile <- tempfile()
tmpdir <- tempdir()
download.file("http://seanlahman.com/files/database/lahman-csv_2014-02-14.zip",tmpfile)
##this shows us files
filenames <- unzip(tmpfile,list=TRUE)
players <- read.csv(unzip(tmpfile,files="Batting.csv",exdir=tmpdir),as.is=TRUE)
unlink(tmpdir)
file.remove(tmpfile)

#1
filter(players,yearID==2012) %>% mutate(AVG=H/AB) %>% filter(AB>=500) %>% select(AVG)

#2
mean(unlist(filter(players,yearID==2012&2013&2014) %>% mutate(AVG=H/AB) %>% filter(AB>=500) %>% select(AVG)))

#3
sd(unlist(filter(players,yearID==2012&2013&2014) %>% mutate(AVG=H/AB) %>% filter(AB>=500) %>% select(AVG)))


#Hierarchical Models in Practice Exercises

# install
BiocManager::install("SpikeInSubset")
library(Biobase)
library(SpikeInSubset)
data(rma95)
y <- exprs(rma95)
pData(rma95)
g <- factor(rep(0:1,each=3))
spike <- rownames(y) %in% colnames(pData(rma95))

#1
library(genefilter)
mean(!spike[(rowttests(y,g)$p.val) < 0.01])

#2
sds <- rowSds(y[,g==0])
index <- paste0( as.numeric(spike), as.numeric(rowttests(y,g)$p.val<0.01))
index <- factor(index,levels=c("11","01","00","10"),labels=c("TP","FP","TN","FN"))
boxplot(split(sds,index))

#3
BiocManager::install("limma")
library(limma)
fit <- lmFit(y, design=model.matrix(~ g))
colnames(coef(fit))
fit <- eBayes(fit)

sampleSD = fit$sigma
posteriorSD = sqrt(fit$s2.post)

LIM <- range( c(posteriorSD,sampleSD))
plot(sampleSD, posteriorSD,ylim=LIM,xlim=LIM)
abline(0,1)
abline(v=sqrt(fit$s2.prior))

#4
fit = lmFit(y, design=model.matrix(~ g))
fit = eBayes(fit)
##second coefficient relates to diffences between group
pvals = fit$p.value[,2] 
mean( !spike[pvals < 0.01] )

#Plots Exercises

#1
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SpikeInSubset")
library(SpikeInSubset)
data(mas133)
e <- exprs(mas133)
plot(e[,1],e[,2],main=paste0("corr=",signif(cor(e[,1],e[,2]),3)),cex=0.5)
k <- 3000
b <- 1000 #a buffer
polygon(c(-b,k,k,-b),c(-b,-b,k,k),col="red",density=0,border="red")

mean(e[,1]<k & e[,2]<k)

#2
plot(log2(e[,1]),log2(e[,2]),main=paste0("corr=",signif(cor(log2(e[,1]),log2(e[,2])),2)),cex=0.5)
k <- log2(3000)
b <- log2(0.5)
polygon(c(b,k,k,b),c(b,b,k,k),col="red",density=0,border="red")

#3
e <- log2(exprs(mas133))
plot((e[,1]+e[,2])/2,e[,2]-e[,1],cex=0.5)

sd(e[,2]-e[,1])

#4
sum(abs(e[,2]-e[,1])>1)
