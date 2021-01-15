#Confounding Exercises

library(dagdata)
data(admissions)
print( admissions )

#1
index = which(admissions$Gender==0)
accepted= sum(admissions$Number[index] * admissions$Percent[index]/100)
applied = sum(admissions$Number[index])
accepted/applied

#2
index = admissions$Gender==1
men = admissions[index,]
women = admissions[!index,]
menYes = sum(men$Number*men$Percent/100)
menNo = sum(men$Number*(1-men$Percent/100))
womenYes = sum(women$Number*women$Percent/100)
womenNo = sum(women$Number*(1-women$Percent/100))
tab = matrix(c(menYes,womenYes,menNo,womenNo),2,2)
chisq.test(tab)$p.value

#3
major = admissions[1:6,1]
men = admissions[1:6,]
women =admissions[7:12,]
H = (men$Number*men$Percent/100 + women$Number*women$Percent/100) / (men$Number+women$Number)
major[which.min(H)]

#4
major = admissions[1:6,1]
men = admissions[1:6,]
women =admissions[7:12,]
H = (men$Number*men$Percent/100 + women$Number*women$Percent/100) / (men$Number+women$Number)
min(H)

#5
cor(H,men$Number)

#6
cor(H,women$Number)


#Confounding in Genomics Exercises

library(Biobase)
library(devtools)
install_github("genomicsclass/GSE5859")

library(GSE5859Subset)

data(GSE5859Subset)
geneExpression = exprs(e)
sampleInfo = pData(e)
#1

year = format(sampleInfo$date,"%y")
length( unique(year) )
tab = table(year,sampleInfo$ethnicity)
print(tab)
x = rowSums( tab != 0)
sum( x >= 2)

#2
month.year = format(sampleInfo$date,"%m%y")
tab = table(month.year,sampleInfo$ethnicity)
print(tab)
x = rowSums( tab != 0)
mean( x >= 2)

#3
year = factor( format(sampleInfo$date,"%y") )
index = which(year%in% c("02","03") & sampleInfo$ethnicity=="CEU")
year = droplevels(year[index])
pval = rowttests(geneExpression[ ,index], year)$p.value
qval = qvalue(pval)
sum(qval$qvalue < 0.05) 

print( qvalue(pval)$pi0 ) 

#4
year = factor( format(sampleInfo$date,"%y") )
index = which(year%in% c("03","04")  & sampleInfo$ethnicity=="CEU")
year = droplevels(year[index])
pval = rowttests(geneExpression[ ,index], year)$p.value
qval = qvalue(pval)
sum(qval$qvalue < 0.05)

#5
index = which(sampleInfo$ethnicity%in% c("CEU","ASN"))
g = droplevels(sampleInfo$ethnicity[index])
pval = rowttests(geneExpression[ ,index], g)$p.value
qval = qvalue(pval)
sum(qval$qvalue < 0.05)

#6
year = factor( format(sampleInfo$date,"%y") )
index = which(sampleInfo$ethnicity%in% c("CEU","ASN") & year=="05")
g = droplevels(sampleInfo$ethnicity[index])
pval = rowttests(geneExpression[ ,index], g)$p.value
qval = qvalue(pval)
sum(qval$qvalue < 0.05)

#7
year = factor( format(sampleInfo$date,"%y") )
index1 = which(sampleInfo$ethnicity=="ASN" & year=="05")
set.seed(3)
index2 = sample( which(sampleInfo$ethnicity == "CEU" & year=="02"), 3)
index = c( index1, index2)
g = droplevels(sampleInfo$ethnicity[index])
pval = rowttests(geneExpression[ ,index], g)$p.value
qval = qvalue(pval)
sum(qval$qvalue < 0.05)


#Modeling Batch Effects Exercises

#1
sex = sampleInfo$group
month = factor( format(sampleInfo$date,"%m"))
table( sampleInfo$group, month)

library(qvalue)
library(genefilter)
pvals = rowttests(geneExpression,factor(sampleInfo$g))$p.value
qvals = qvalue(pvals)$qvalues
sum(qvals<0.1)

#2
index = geneAnnotation$CHR[qvals<0.1]%in%c("chrX","chrY")
mean(index)

#3
index = which(qvals<0.1 & !geneAnnotation$CHR%in%c("chrX","chrY"))
month = factor( format(sampleInfo$date,"%m"))
pvals = rowttests(geneExpression[index,],month)$p.value
mean(pvals<0.05)

#5
X = model.matrix(~sex+month)
i = 234
y = geneExpression[i,]
fit = lm(y~X-1)
summary(fit)$coef

pvals = sapply(1:nrow(geneExpression),function(i){
  y = geneExpression[i,]
  fit = lm(y~X-1) 
  summary(fit)$coef[2,4]
})
qvals = qvalue(pvals)$qvalue
sum(qvals<0.1)

#6
index = geneAnnotation$CHR[qvals<0.1]%in%c("chrX","chrY")
mean(index)

#7
pvals = sapply(1:nrow(geneExpression),function(i){
  y = geneExpression[i,]
  fit = lm(y~X)
  summary(fit)$coef[3,4]
})
qvals = qvalue(pvals)$qvalue
sum(qvals<0.1)


#Factor Analysis Exercises

#1
y = geneExpression - rowMeans(geneExpression)

sex = sampleInfo$group
mypar(1,2)
cors = cor(y)
image(cors)
o = order(sampleInfo$date)
image(cors[o,o])

#3
pcs = svd(y)$v[,1:2]
o = order(sampleInfo$date)
cols = as.numeric(month)[o]
mypar(2,1)
for(i in 1:2){
  plot(pcs[o,i],col=cols,xaxt="n",xlab="")
  label = gsub("2005-","",sampleInfo$date[o])
  axis(1,1:ncol(y),label,las=2)
}

#4
s = svd(y)
varexplained = s$d^2/ sum(s$d^2)
plot(varexplained)
sum(varexplained>0.10)

#5
month = factor( format(sampleInfo$date,"%m"))
cors = cor( as.numeric(month),s$v)
plot(t(cors))
which.max(abs(cors))

max(abs(cors))

#6
sex = sampleInfo$group
cors = cor( sex,s$v)
plot(t(cors))
which.max(abs(cors))

max(abs(cors))

#7
X <- model.matrix(~sex+s$v[,1:2])

pvals = sapply(1:nrow(geneExpression),function(i){
  y = geneExpression[i,]
  fit = lm(y~X-1)
  summary(fit)$coef[2,4]
})
qvals = qvalue(pvals)$qvalue
sum(qvals<0.1) 

index = geneAnnotation$CHR[qvals<0.1]%in%c("chrX","chrY")
mean(index)


#SVA Exercises

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("sva")
library(sva)

#1
s <- svd(geneExpression-rowMeans(geneExpression))
cor(sampleInfo$group,s$v[,1])

sex = sampleInfo$group
mod = model.matrix(~sex)
svafit = sva(geneExpression,mod)
head(svafit$sv)

for(i in 1:ncol(svafit$sv)){
  print( cor(s$v[,i],svafit$sv[,i]) )
}

X= model.matrix(~sex+svafit$sv)
pvals = sapply(1:nrow(geneExpression),function(i){
  y = geneExpression[i,]
  fit = lm(y~X-1)
  summary(fit)$coef[2,4]
})
qvals = qvalue(pvals)$qvalue
sum(qvals<0.1)

#2
index = geneAnnotation$CHR[qvals<0.1]%in%c("chrX","chrY")
print(mean(index))
res = sapply(1:nrow(geneExpression),function(i){
  y = geneExpression[i,]
  fit = lm(y~X-1)
  summary(fit)$coef[2,c(1,4)]
})
qvals = qvalue(res[2,])$qvalue
pcutoff = max( res[2,qvals < .1] )
library(rafalib)
mypar(1,1)
plot(res[1,],-log10(res[2,]),xlab="M",ylab="log10 p-value")
ind = which(geneAnnotation$CHR=="chrY")
points(res[1,ind],-log10(res[2,ind]),col=1,pch=16)
ind = which(geneAnnotation$CHR=="chrX")
points(res[1,ind],-log10(res[2,ind]),col=2,pch=16)

abline(h=-log10(pcutoff))
legend("bottomleft",c("chrX","chrY"),col=c(2,1),pch=16)
