dat <- read.csv("femaleMiceWeights.csv")

#Getting Started Exercises #1
head(dat)

#Getting Started Exercises #2
dat[12, 2]

#Getting Started Exercises #3
dat$Bodyweight[11]

#Getting Started Exercises #4
length(dat$Bodyweight)

#Getting Started Exercises #5
mean(dat$Bodyweight[dat$Diet=='hf'])

#Getting Started Exercises #6
set.seed(1)
dat$Bodyweight[sample(13:24,1)]

library(dplyr)

controls <- filter(dat, Diet == 'chow')
View(controls)

controls <- select(controls, Bodyweight)
unlist(controls)

controls <- dat %>% filter(Diet=='chow') %>% select(Bodyweight) %>% unlist
controls
mean(controls)


library(downloader)
url="https://raw.githubusercontent.com/genomicsclass/dagdata/master/inst/extdata/msleep_ggplot2.csv"
filename <- basename(url)
download(url,filename)

#dplyr Exercises #1
p <- read.csv('msleep_ggplot2.csv')
class(p)

#dplyr Exercises #2
head(p)
nrow(filter(p, order=='Primates'))

#dplyr Exercises #3
class(filter(p, order=='Primates'))

#dplyr Exercises #4
p %>% filter(order=='Primates') %>% select(sleep_total) %>% class

#dplyr Exercises #5
mean(p %>% filter(order=='Primates') %>% select(sleep_total) %>% unlist)

#dplyr Exercises #6
p %>% filter(order=='Primates') %>% summarize(mean(sleep_total))
