#Distance Exercises

#1
library(devtools)
install_github("genomicsclass/tissuesGeneExpression")
library(tissuesGeneExpression)
data(tissuesGeneExpression)
head(e)
head(tissue)

sum(tissue == 'hippocampus')

#2
sqrt(crossprod(e[,3]-e[,45]))

#3
sqrt(crossprod(e["210486_at",]-e["200805_at",]))

#4
nrow(e)^2

#5
d = dist(t(e))
length(d)
