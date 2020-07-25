###An Example of High Throughput Data

install.packages('devtools')
library(devtools)
install_github("genomicsclass/GSE5859Subset")

library(GSE5859Subset)
data(GSE5859Subset) ##this loads the three tables
dim(geneExpression)

#To see if geneExpression matrix columns match rows of sample info matrix:
identical(colnames(geneExpression), sampleInfo$filename)

#If we wanted to know what group each column of geneExpression matrix comes from:
sampleInfo$group #We see first 12 columns are cases and last 12- controls.

#We can get information/features about genes from geneAnnotaion matrix:
head(geneAnnotation)


###R Refresher Exercises

#1
sum(sampleInfo$date=="2005-06-27")

#2
sum(na.omit(geneAnnotation$CHR) == 'chrY')
#or
sum(geneAnnotation$CHR=="chrY",na.rm=TRUE)

#3
i <- which(geneAnnotation$SYMBOL=="ARPC1A")
j <- which(sampleInfo$date=="2005-06-10")
geneExpression[i,j]

#4
median(apply(geneExpression, 2, median))

#5
func <- function(e, group){
  return(t.test( e[group==1], e[group==0])$p.value)
}
g <- factor(sampleInfo$group)
min(apply(geneExpression, 1, func, group=g))


###The challenge of multiple testing

g <- sampleInfo$group

#Taking 25th gene:
e <- geneExpression[25,]

#Comparing control and case:
library(rafalib)
mypar(1,2)

qqnorm(e[g==1])
qqline(e[g==1])

qqnorm(e[g==0])
qqline(e[g==0])

#Making t test:
t.test(e[g==1], e[g==0])$p.value #it's not smaller than 0.05, that makes is insignificant

#Writing function for generating p value with t test:
mytest <- function(x){
  t.test(x[g==1], x[g==0])$p.value
}
mytest(geneExpression[25,])

#Applying above function to all genes in geneExpression matrix:
pvalues <- apply(geneExpression, 1, mytest)

#Finding how many p values of genes are lower than 0.05:
sum(pvalues <= 0.05)

#Monte Carlo simulation:
set.seed(1)
m <- nrow(geneExpression)
n <- ncol(geneExpression)
randomData <- matrix(rnorm(n*m),m,n)
nullpvals <- apply(randomData,1,mytest)
sum(nullpvals<0.05)


###p-values Are Random Variables

#Producing random variables with simulation:
filename <- "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
set.seed(1)
population <- unlist( read.csv(filename) )
N <- 12
B <- 10000
pvals <- replicate(B,{
  control = sample(population,N)
  treatment = sample(population,N)
  t.test(treatment,control)$p.val 
})
hist(pvals)


###Inference in Practice Exercises

set.seed(1)
library(downloader)
url = "https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/femaleControlsPopulation.csv"
filename = "femaleControlsPopulation.csv"
if (!file.exists(filename)) download(url,destfile=filename)
population = read.csv(filename)
pvals <- replicate(1000,{
  control = sample(population[,1],12)
  treatment = sample(population[,1],12)
  t.test(treatment,control)$p.val
})
head(pvals)
hist(pvals)

#1
sum(pvals <= 0.05) / length(pvals)

#2
mean(pvals<= 0.01)

#3
cases = rnorm(10,30,2)
controls = rnorm(10,30,2)
t.test(cases,controls)$p.value

set.seed(100)
pvals <- replicate(20,{
  cases = rnorm(10,30,2)
  controls = rnorm(10,30,2)
  t.test(cases,controls)$p.val
})

sum(pvals<=0.05)

#4
set.seed(100)
B = 1000
plessthan = replicate(B,{
  pvals = replicate(20,{
    cases = rnorm(10,30,2)
    controls = rnorm(10,30,2)
    t.test(cases,controls)$p.value
  })
  sum(pvals<=0.05)
})
mean(plessthan)

#5
mean(plessthan>0)