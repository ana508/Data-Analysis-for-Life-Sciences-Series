library(rafalib)
library(downloader)
library(dplyr)

###Median, MAD, and Spearman Correlation Exercises

data(ChickWeight)
head(ChickWeight)
plot( ChickWeight$Time, ChickWeight$weight, col=ChickWeight$Diet)

#The meaning of this line is: reshape the data from long to wide, where the columns Chick and Diet are the ID's and the column Time indicates different observations for each ID. Now examine the head of this dataset:
chick <- reshape(ChickWeight, idvar=c("Chick","Diet"), timevar="Time",
                direction="wide")

head(chick)

#We also want to remove any chicks that have missing observations at any time points (NA for "not available"). The following line of code identifies these rows and then removes them:
chick <- na.omit(chick)

#1
mean(c(chick$weight.4, 3000))/mean(chick$weight.4)

#2
median(c(chick$weight.4, 3000))/median(chick$weight.4)

#3
sd(c(chick$weight.4, 3000))/sd(chick$weight.4)

#4
mad(c(chick$weight.4, 3000))/mad(chick$weight.4)

#5
# The Pearson correlation between x and y is given in R by cor(x,y). The Spearman correlation is given by cor(x,y,method="spearman").
cor(c(chick$weight.4,3000),c(chick$weight.21,3000)) / cor(chick$weight.4,chick$weight.21)

###Mann-Whitney-Wilcoxon Test Exercises

#1
x <- chick$weight.4[chick$Diet == 1]
y <- chick$weight.4[chick$Diet == 4]
t.test(c(x,200),y)$p.value

#2
wilcox.test(c(x, 200), y, exact=FALSE)$p.value

#3
mypar(1,3)
boxplot(x,y)
boxplot(x,y+10)
boxplot(x,y+100)

t.test(x,y+10)$statistic - t.test(x,y+100)$statistic

#4

#Examine the Wilcoxon test statistic for x and y+10 and for x and y+100. 
#Because the Wilcoxon works on ranks, once the two groups show complete separation, 
#that is all points from group y are above all points from group x, the statistic will 
#not change, regardless of how large the difference grows. Likewise, the p-value has a 
#minimum value, regardless of how far apart the groups are. This means that the Wilcoxon 
#test can be considered less powerful than the t-test in certain contexts. In fact, for small 
#sample sizes, the p-value can't be very small, even when the difference is very large.

wilcox.test(c(1,2,3),c(4,5,6))$p.value

#5
wilcox.test(c(1,2,3),c(400,500,600))$p.value